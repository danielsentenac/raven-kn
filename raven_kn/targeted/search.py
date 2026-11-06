"""Light targeted search for a compact-binary signal at a known sky position.

Per sub-segment of the on-source window (default 1024 s):

1. strain from GWOSC, local files, or any callable; resampled, high-passed;
   optional gating of loud glitches (``pycbc.strain.detect_loud_glitches``);
2. PSD by median Welch over the sub-segment's data;
3. matched filter of each template (non-spinning equal-mass TaylorF2 over a
   chirp-mass range) in each detector;
4. each detector's |SNR| shifted by its time delay from the geocentre for the
   known (RA, Dec) at that time, max-pooled in short bins (absorbs the
   rounding of the delays), and combined as rho_net = sqrt(sum_d rho_d^2);
5. the loudest rho_net over templates and the sub-segment's on-source time;
6. background: cyclic time slides of the other detectors (steps larger than
   the light travel time) give the loudest rho_net of noise-only trials of
   the same length; p-value = (n_louder + 1)/(n_slides + 1).

Sub-segment p-values are combined with a trials factor
(:func:`raven_kn.treatment.far.combine_subsegments`).

Limitations (v0): incoherent network statistic, no chi-square / signal
consistency test, aligned zero-spin equal-mass templates. Use PyGRB for
production; this is for fast assessment and validation.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Optional, Sequence

import numpy as np

from ..treatment.far import combine_subsegments, targeted_pvalue

DEFAULT_IFOS = ("H1", "L1", "V1")


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

@dataclass
class TargetedSearchConfig:
    ra_deg: float
    dec_deg: float
    gps_start: float
    gps_end: float
    ifos: Sequence[str] = DEFAULT_IFOS
    mchirp_min: float = 1.0
    mchirp_max: float = 1.5
    f_low: float = 30.0
    sample_rate: int = 2048
    subsegment_s: float = 1024.0
    pad_s: float = 8.0
    min_match: float = 0.97
    n_slides: int = 200
    slide_step_s: float = 3.0
    pool_s: float = 1.0 / 128
    gate_threshold: Optional[float] = 50.0
    max_templates: int = 20000
    label: str = "targeted"

    @classmethod
    def from_request(cls, path: str | Path, **overrides) -> "TargetedSearchConfig":
        """From ``<name>_targeted_search_request.json`` written by ``raven_kn.search.assess``."""
        r = json.loads(Path(path).read_text())
        kw = dict(ra_deg=r["ra_deg"], dec_deg=r["dec_deg"], gps_start=r["gps_start"],
                  gps_end=r["gps_end"], label=r.get("candidate", "targeted"))
        if "mchirp_range" in r:
            kw["mchirp_min"], kw["mchirp_max"] = r["mchirp_range"]
        kw.update(overrides)
        return cls(**kw)


# --------------------------------------------------------------------------
# Template bank
# --------------------------------------------------------------------------

def _mass_equal(mchirp: float) -> float:
    return mchirp * 2 ** (1 / 5)          # m1 = m2 for chirp mass mchirp


def template_duration(mchirp: float, f_low: float) -> float:
    """Leading-order inspiral time from f_low (s)."""
    G_c3 = 4.925490947e-6                  # G Msun / c^3 in s
    return 5 / 256 * (G_c3 * mchirp) ** (-5 / 3) * (math.pi * f_low) ** (-8 / 3)


def _chirp_metric(mchirp: float, f_low: float, f_high: float, psd_func) -> float:
    """Fisher metric g on ln(Mc) for the leading-order TaylorF2 phase,
    psi(f) = 3/128 (pi G Mc f / c^3)^(-5/3), after maximising over the
    coalescence time and phase (projecting out 1 and f). Mismatch for a
    step d ln Mc is ~ g d^2 / 2."""
    G_c3 = 4.925490947e-6
    f = np.linspace(f_low, f_high, 20000)
    w = f ** (-7 / 3) / psd_func(f)
    w /= np.trapezoid(w, f)
    dpsi = -5 / 3 * 3 / 128 * (np.pi * G_c3 * mchirp * f) ** (-5 / 3)
    basis = np.vstack([np.ones_like(f), 2 * np.pi * f])          # d psi / d phi_c, d psi / d t_c
    gram = np.array([[np.trapezoid(w * a * b, f) for b in basis] for a in basis])
    proj = np.array([np.trapezoid(w * dpsi * a, f) for a in basis])
    return float(np.trapezoid(w * dpsi ** 2, f) - proj @ np.linalg.solve(gram, proj))


def _aligo_psd():
    import pycbc.psd
    df = 0.125
    flen = int(4096 / df)
    psd = pycbc.psd.aLIGOZeroDetHighPower(flen, df, 9.0).numpy()
    fr = np.arange(flen) * df

    def f(x):
        return np.interp(x, fr, psd)
    return f


def _direct_match(mc_a: float, mc_b: float, f_low: float, sample_rate: int) -> float:
    import pycbc.psd
    from pycbc.filter import match
    from pycbc.waveform import get_fd_waveform
    dur = template_duration(min(mc_a, mc_b), f_low)
    df = 1.0 / 2 ** int(math.ceil(math.log2(dur + 16)))
    flen = int(sample_rate / 2 / df) + 1
    psd = pycbc.psd.aLIGOZeroDetHighPower(flen, df, f_low)

    def wf(mc):
        m = _mass_equal(mc)
        hp, _ = get_fd_waveform(approximant="TaylorF2", mass1=m, mass2=m, delta_f=df, f_lower=f_low)
        hp.resize(flen)
        return hp
    return float(match(wf(mc_a), wf(mc_b), psd=psd, low_frequency_cutoff=f_low)[0])


def build_bank(mchirp_min: float, mchirp_max: float, f_low: float = 30.0,
               min_match: float = 0.97, sample_rate: int = 2048,
               max_templates: int = 20000, calibrate: bool = True) -> list[float]:
    """Chirp masses of equal-mass, non-spinning TaylorF2 templates spaced with
    the analytic chirp-mass metric (aLIGO design PSD) so that a signal
    half-way between two neighbours keeps a match ``min_match`` (minimal
    match). The leading-order metric underestimates the mismatch at these
    step sizes; the step factor is calibrated once, at ``mchirp_min``, by
    bisection on direct TaylorF2 matches (it is nearly independent of Mc).
    Upper frequency: ISCO of the equal-mass system, capped at Nyquist."""
    if mchirp_min >= mchirp_max:
        return [mchirp_min]
    psd = _aligo_psd()

    def step(mc, c):
        f_isco = 4400.0 / (2 * _mass_equal(mc))
        g = _chirp_metric(mc, f_low, min(f_isco, sample_rate / 2), psd)
        return c * 2 * math.sqrt(2 * (1 - min_match) / g)

    c = 1.0
    if calibrate:
        lo, hi = 0.3, 1.5
        for _ in range(10):
            mid = 0.5 * (lo + hi)
            d = step(mchirp_min, mid)
            if _direct_match(mchirp_min, mchirp_min * math.exp(d / 2), f_low, sample_rate) >= min_match:
                lo = mid
            else:
                hi = mid
        c = lo
    bank = [mchirp_min]
    while bank[-1] < mchirp_max and len(bank) < max_templates:
        bank.append(bank[-1] * math.exp(step(bank[-1], c)))
    if len(bank) > 1 and bank[-1] > mchirp_max:
        bank[-1] = mchirp_max
    return bank


# --------------------------------------------------------------------------
# Data access
# --------------------------------------------------------------------------

StrainFetcher = Callable[[str, float, float], "object"]   # (ifo, start, end) -> pycbc TimeSeries


def gwosc_fetcher(sample_rate: int = 2048, ca_bundle: Optional[str] = None) -> StrainFetcher:
    """Open data from GWOSC via gwpy. ``ca_bundle``: CA file for systems
    where Python's certificate store fails against gwosc.org."""
    if ca_bundle:
        os.environ.setdefault("REQUESTS_CA_BUNDLE", ca_bundle)
        os.environ.setdefault("SSL_CERT_FILE", ca_bundle)

    def fetch(ifo, start, end):
        from gwpy.timeseries import TimeSeries
        ts = TimeSeries.fetch_open_data(ifo, start, end, cache=True)
        return ts.resample(sample_rate).to_pycbc()
    return fetch


def file_fetcher(files: dict[str, str | Path], sample_rate: int = 2048) -> StrainFetcher:
    """Local GWOSC hdf5 / gwf files, one per detector."""
    def fetch(ifo, start, end):
        from gwpy.timeseries import TimeSeries
        f = str(files[ifo])
        if f.endswith((".hdf5", ".h5")):
            ts = TimeSeries.read(f, format="hdf5.gwosc")
        else:
            ts = TimeSeries.read(f, f"{ifo}:GWOSC-4KHZ_R1_STRAIN")
        return ts.crop(start, end).resample(sample_rate).to_pycbc()
    return fetch


# --------------------------------------------------------------------------
# Search
# --------------------------------------------------------------------------

@dataclass
class SubsegmentResult:
    gps_start: float
    gps_end: float
    ifos: list[str]
    loudest_rho_net: float
    loudest_gps: float
    loudest_mchirp: float
    loudest_rho_ifo: dict[str, float]
    background: list[float] = field(default_factory=list)
    p_value: float = 1.0
    gates: dict[str, list[float]] = field(default_factory=dict)
    notes: str = ""


@dataclass
class TargetedSearchResult:
    config: dict
    bank_size: int
    subsegments: list[SubsegmentResult]
    p_value: float
    loudest: dict

    def to_json(self, path: str | Path) -> None:
        d = asdict(self)
        for s in d["subsegments"]:
            s["background"] = [round(x, 3) for x in s["background"]]
        Path(path).write_text(json.dumps(d, indent=2))


def _delays(ifos, ra_deg, dec_deg, gps):
    import lal
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    t = lal.LIGOTimeGPS(float(gps))
    return {i: lal.TimeDelayFromEarthCenter(lal.cached_detector_by_prefix[i].location, ra, dec, t)
            for i in ifos}


def _condition(ts, f_low, gate_threshold):
    from pycbc.filter import highpass
    from pycbc.strain import detect_loud_glitches, gate_data
    ts = highpass(ts, max(f_low - 10, 10.0))
    gates = []
    if gate_threshold:
        try:
            gates = list(detect_loud_glitches(ts, threshold=gate_threshold,
                                              low_freq_cutoff=f_low))
        except Exception:
            gates = []
        if gates:
            ts = gate_data(ts, [(float(g), 0.5, 0.25) for g in gates])
    return ts, [float(g) for g in gates]


def _pool(x, n):
    m = len(x) // n
    return x[: m * n].reshape(m, n).max(axis=1)


def search_subsegment(cfg: TargetedSearchConfig, bank: Sequence[float], fetch: StrainFetcher,
                      start: float, end: float) -> SubsegmentResult:
    import pycbc.psd
    from pycbc.filter import matched_filter
    from pycbc.waveform import get_fd_waveform

    tmax = template_duration(min(bank), cfg.f_low)
    lead = tmax + 4 * cfg.pad_s                      # corrupted start of the matched filter
    d_start, d_end = start - lead, end + 4 * cfg.pad_s
    data, gates, ifos = {}, {}, []
    for ifo in cfg.ifos:
        try:
            ts = fetch(ifo, d_start, d_end)
        except Exception as e:                       # no data for this detector
            gates[ifo] = [f"no data: {type(e).__name__}"]
            continue
        if not np.all(np.isfinite(ts.numpy())):
            gates[ifo] = ["non-finite data"]
            continue
        ts, g = _condition(ts, cfg.f_low, cfg.gate_threshold)
        data[ifo], gates[ifo] = ts, g
        ifos.append(ifo)
    if len(ifos) < 2:
        return SubsegmentResult(start, end, ifos, 0.0, start, 0.0, {}, [], 1.0, gates,
                                "fewer than two detectors with data")

    fs = cfg.sample_rate
    n_pool = max(int(round(cfg.pool_s * fs)), 1)
    delays = _delays(ifos, cfg.ra_deg, cfg.dec_deg, 0.5 * (start + end))
    stilde, psds = {}, {}
    for ifo in ifos:
        ts = data[ifo]
        psd = pycbc.psd.welch(ts, seg_len=int(16 * fs), seg_stride=int(8 * fs),
                              avg_method="median")
        psd = pycbc.psd.interpolate(psd, 1.0 / ts.duration)
        psds[ifo] = pycbc.psd.inverse_spectrum_truncation(psd, int(4 * fs),
                                                          low_frequency_cutoff=cfg.f_low)
        stilde[ifo] = ts.to_frequencyseries()

    n_on = None
    pooled = {ifo: [] for ifo in ifos}               # per template, pooled |snr| aligned to geocentre
    for mc in bank:
        m = _mass_equal(mc)
        for ifo in ifos:
            st = stilde[ifo]
            hp, _ = get_fd_waveform(approximant="TaylorF2", mass1=m, mass2=m,
                                    delta_f=st.delta_f, f_lower=cfg.f_low)
            hp.resize(len(st))
            snr = matched_filter(hp, st, psd=psds[ifo], low_frequency_cutoff=cfg.f_low)
            t0 = float(snr.start_time)
            # geocentre time tg corresponds to detector time tg + delay
            i0 = int(round((start + delays[ifo] - t0) * fs))
            i1 = int(round((end + delays[ifo] - t0) * fs))
            x = np.abs(snr.numpy()[i0:i1]).astype(np.float32)
            p = _pool(x, n_pool)
            n_on = len(p) if n_on is None else min(n_on, len(p))
            pooled[ifo].append(p)
    P = {ifo: np.stack([p[:n_on] for p in pooled[ifo]]) for ifo in ifos}   # (n_templates, n_bins)

    def network(shifts):
        tot = None
        for k, ifo in enumerate(ifos):
            a = P[ifo] if shifts[k] == 0 else np.roll(P[ifo], shifts[k], axis=1)
            tot = a.astype(np.float64) ** 2 if tot is None else tot + a ** 2
        return np.sqrt(tot)

    net = network([0] * len(ifos))
    it, ib = np.unravel_index(int(np.argmax(net)), net.shape)
    loud = float(net[it, ib])
    loud_gps = start + (ib + 0.5) * n_pool / fs
    rho_ifo = {ifo: float(P[ifo][it, ib]) for ifo in ifos}

    step = max(int(round(cfg.slide_step_s / (n_pool / fs))), 1)
    n_slides = min(cfg.n_slides, max(n_on // step - 1, 0))
    bkg = []
    for k in range(1, n_slides + 1):
        # detector j slid by j*k*step bins: distinct relative shifts for all pairs
        bkg.append(float(network([j * k * step for j in range(len(ifos))]).max()))
    p = targeted_pvalue(loud, bkg) if bkg else 1.0
    return SubsegmentResult(start, end, ifos, loud, loud_gps, float(bank[it]), rho_ifo, bkg, p,
                            gates, f"{len(bank)} templates, {n_slides} slides, delays "
                            + ", ".join(f"{i}:{delays[i] * 1e3:+.1f} ms" for i in ifos))


def run(cfg: TargetedSearchConfig, fetch: StrainFetcher, out_dir: Optional[str | Path] = None,
        bank: Optional[Sequence[float]] = None, progress: bool = True) -> TargetedSearchResult:
    bank = list(bank) if bank is not None else build_bank(
        cfg.mchirp_min, cfg.mchirp_max, cfg.f_low, cfg.min_match, cfg.sample_rate, cfg.max_templates)
    edges = np.arange(cfg.gps_start, cfg.gps_end, cfg.subsegment_s).tolist() + [cfg.gps_end]
    subs = []
    for a, b in zip(edges[:-1], edges[1:]):
        if b - a < 16:
            continue
        r = search_subsegment(cfg, bank, fetch, a, b)
        subs.append(r)
        if progress:
            print(f"[{a:.0f}, {b:.0f}] rho_net={r.loudest_rho_net:.2f} at {r.loudest_gps:.2f} "
                  f"Mc={r.loudest_mchirp:.4f} p={r.p_value:.4f} {r.notes}", flush=True)
    valid = [s for s in subs if s.background]
    p_all = combine_subsegments([s.p_value for s in valid]) if valid else 1.0
    best = max(subs, key=lambda s: s.loudest_rho_net) if subs else None
    loudest = asdict(best) if best else {}
    loudest.pop("background", None)
    res = TargetedSearchResult(asdict(cfg), len(bank), subs, p_all, loudest)
    if out_dir is not None:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        res.to_json(out / f"{cfg.label}_targeted_result.json")
    return res
