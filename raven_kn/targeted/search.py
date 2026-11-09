"""Light targeted search for a compact-binary signal at a known sky position.

Per sub-segment of the on-source window (default 1024 s):

1. strain from GWOSC, local files, or any callable; resampled, high-passed;
   optional gating of loud glitches (``pycbc.strain.detect_loud_glitches``);
2. PSD by median Welch over the sub-segment's data;
3. matched filter of each template (non-spinning equal-mass TaylorF2 over a
   chirp-mass range) in each detector;
4. each detector's |SNR| shifted by its time delay from the geocentre for the
   known (RA, Dec) at that time, max-pooled in short bins (absorbs the
   rounding of the delays), and combined as rho_net = sqrt(sum_d rho_d^2)
   where at least two detectors exceed ``coinc_threshold`` (else 0);
5. the loudest rho_net over templates and the sub-segment's on-source time;
6. background: cyclic time slides of the other detectors (steps larger than
   the light travel time) give the loudest rho_net of noise-only trials of
   the same length, after removing the loudest on-source candidate
   (hierarchical removal); p-value = (n_louder + 1)/(n_slides + 1), plus an
   exponential-tail extrapolation below that floor.

Sub-segment p-values are combined with a trials factor
(:func:`raven_kn.treatment.far.combine_subsegments`).

Signal consistency: where a detector's SNR exceeds ``chisq_threshold``, the
power chi-square (``chisq_bins`` bins) is computed and the SNR replaced by
PyCBC's reweighted SNR (newSNR, q = 6, n = 2), on-source and in the slides.

Coherent statistic (``statistic="coherent"``, as PyGRB): each detector's
template is shifted by its exact (fractional-sample) delay, so the complex
SNRs z_d line up at the geocentre sample by sample. At samples where at least
two detectors exceed ``coinc_threshold``, z is projected onto the plane of
the network's response to the two polarisations, w_+ = sigma_d F+_d and
w_x = sigma_d Fx_d (sigma_d the template's sensitivity in detector d, F the
antenna patterns at the known sky position; the plane does not depend on
the polarisation angle): rho_coh^2 = z^H P z, rho_null^2 = |z|^2 -
rho_coh^2. With two detectors P = 1 (rho_coh = rho_net, no null stream);
with three, rho_null (2 degrees of freedom in noise) rejects candidates
whose phases and amplitudes no astrophysical signal from that direction can
produce. Candidates are cut, as in PyGRB, if rho_null > ``null_min`` (or
``null_min + null_grad (rho_coh - null_step)`` above ``null_step``). Per
detector newSNR reweighting scales z_d before projection. Time slides shift
whole samples; the loudest on-source candidate is removed from them by 1 s
blocks.

Limitations: zero-spin equal-mass templates. Use PyGRB for production; this
is for fast assessment and validation.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import asdict, dataclass, field, replace
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
    chisq_bins: Optional[int] = 16          # None: no chi-square reweighting
    chisq_threshold: float = 5.0             # compute chi-square where SNR exceeds this
    max_chisq_points: int = 500              # per template and detector (loudest kept)
    hierarchical_removal: bool = True        # mask the loudest on-source candidate in the slides
    coinc_threshold: float = 4.0             # a candidate needs >= 2 detectors above this (0: off)
    removal_window_s: float = 1.0
    removal_min_second_snr: float = 6.0      # remove only if the 2nd-loudest detector is this loud
    max_templates: int = 20000
    n_offsource: int = 0                     # off-source trials per sub-segment (0: slides only)
    offsource_gap_s: float = 64.0            # gap between the on-source window and off-source trials
    statistic: str = "network"               # "network" (incoherent) or "coherent" (known sky position)
    null_min: float = 5.25                   # coherent: null-SNR cut (PyGRB's --null-min/-grad/-step)
    null_grad: float = 0.2
    null_step: float = 20.0
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


def file_fetcher(files: dict[str, str | Path], sample_rate: int = 2048,
                 fallback: Optional[StrainFetcher] = None) -> StrainFetcher:
    """Local strain files per detector: GWOSC hdf5, or frame files given as
    ``path:channel`` (e.g. ``L-L1_LOSC_CLN_4_V1-1187007040-2048.gwf:L1:LOSC-STRAIN``).
    Detectors without a file use ``fallback`` (e.g. :func:`gwosc_fetcher`)."""
    def fetch(ifo, start, end):
        if ifo not in files:
            if fallback is None:
                raise KeyError(f"no file for {ifo}")
            return fallback(ifo, start, end)
        spec = str(files[ifo])
        if spec.endswith((".hdf5", ".h5")):
            from gwpy.timeseries import TimeSeries
            ts = TimeSeries.read(spec, format="hdf5.gwosc").crop(start, end).to_pycbc()
        else:
            from pycbc.frame import read_frame
            path, channel = spec.split(":", 1)
            ts = read_frame(path, channel, start_time=start, end_time=end)
        from pycbc.filter import resample_to_delta_t
        return resample_to_delta_t(ts, 1.0 / sample_rate)
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
    p_value_extrapolated: Optional[float] = None
    extrapolation_reach: Optional[float] = None
    offsource: list[float] = field(default_factory=list)
    p_value_offsource: Optional[float] = None
    statistic: str = "network"               # what loudest_rho_net ranks by
    loudest_rho_coh: Optional[float] = None  # coherent statistic only
    loudest_rho_null: Optional[float] = None
    loudest_rho_network: Optional[float] = None


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


def extrapolated_pvalue(loudest: float, background: Sequence[float], tail_frac: float = 0.1) -> float:
    """p-value of ``loudest`` from an exponential fit to the top ``tail_frac``
    of the background loudest values: P(X > x) = q exp(-(x - x_q)/lambda)
    for x > x_q, x_q the (1 - q) quantile, lambda the mean excess (maximum
    likelihood). Gives a p-value below the 1/(N+1) floor of the slides,
    as an extrapolation (Jamsin et al. extrapolate their background
    similarly). Returns the empirical value when ``loudest`` is inside the
    sampled range."""
    b = np.sort(np.asarray(background, float))
    n = len(b)
    k = max(int(round(tail_frac * n)), 5)
    xq = b[-k - 1] if n > k else b[0]
    lam = float(np.mean(b[-k:] - xq))
    if loudest <= b[-1] or lam <= 0:
        return targeted_pvalue(loudest, b)
    return float(k / n * np.exp(-(loudest - xq) / lam))


def extrapolation_reach(loudest: float, background: Sequence[float]) -> float:
    """How far beyond the sampled background the extrapolation goes, in units
    of the background's spread above its median: values >> 3 mean the
    extrapolated p-value is a formal number only."""
    b = np.asarray(background, float)
    if len(b) == 0:
        return float("inf")
    spread = max(float(np.max(b) - np.median(b)), 1e-6)
    return float((loudest - np.max(b)) / spread)


def slide_shifts(n_ifos: int, n_steps: int, n_slides: int) -> list[tuple[int, ...]]:
    """Distinct time-slide offsets (in steps) for each detector, first fixed.

    Two detectors: (0, k), k = 1 .. n_steps - 1. Three or more: the other
    detectors take independent offsets with all relative offsets non-zero,
    so many more slides fit in one sub-segment (n_steps^2 instead of n_steps).
    """
    out = []
    if n_ifos < 2 or n_steps < 2:
        return out
    if n_ifos == 2:
        return [(0, k) for k in range(1, min(n_steps, n_slides + 1))]
    rng = np.random.default_rng(12345)
    seen = set()
    tries = 0
    while len(out) < n_slides and tries < 50 * n_slides:
        tries += 1
        sh = (0,) + tuple(int(x) for x in rng.integers(1, n_steps, n_ifos - 1))
        rel = [(a - b) % n_steps for i, a in enumerate(sh) for b in sh[i + 1:]]
        if 0 in rel or sh in seen:
            continue
        seen.add(sh)
        out.append(sh)
    return out


def search_subsegment(cfg: TargetedSearchConfig, bank: Sequence[float], fetch: StrainFetcher,
                      start: float, end: float) -> SubsegmentResult:
    import pycbc.psd
    from pycbc.events.ranking import newsnr
    from pycbc.filter import matched_filter_core
    from pycbc.vetoes import power_chisq_at_points_from_precomputed, power_chisq_bins
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
        if len(ts) % 2:                              # even length: FFT and PSD sizes agree
            ts = ts[:len(ts) - 1]
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
        stilde[ifo] = ts.to_frequencyseries()
        psd = pycbc.psd.interpolate(psd, stilde[ifo].delta_f)
        if len(psd) != len(stilde[ifo]):
            psd.resize(len(stilde[ifo]))
        psds[ifo] = pycbc.psd.inverse_spectrum_truncation(psd, int(4 * fs),
                                                          low_frequency_cutoff=cfg.f_low)
    if cfg.statistic == "coherent":
        return _coherent_subsegment(cfg, bank, ifos, stilde, psds, delays, gates, start, end)
    if cfg.statistic != "network":
        raise ValueError(f"unknown statistic {cfg.statistic!r}")

    n_on = None
    n_chisq = {ifo: 0 for ifo in ifos}
    pooled = {ifo: [] for ifo in ifos}               # per template, pooled |snr| aligned to geocentre
    for mc in bank:
        m = _mass_equal(mc)
        for ifo in ifos:
            st = stilde[ifo]
            hp, _ = get_fd_waveform(approximant="TaylorF2", mass1=m, mass2=m,
                                    delta_f=st.delta_f, f_lower=cfg.f_low)
            hp.resize(len(st))
            snr_u, corr, norm = matched_filter_core(hp, st, psd=psds[ifo],
                                                    low_frequency_cutoff=cfg.f_low)
            t0 = float(snr_u.start_time)
            # geocentre time tg corresponds to detector time tg + delay
            i0 = int(round((start + delays[ifo] - t0) * fs))
            i1 = int(round((end + delays[ifo] - t0) * fs))
            su = snr_u.numpy()
            x = (np.abs(su[i0:i1]) * norm).astype(np.float32)
            if cfg.chisq_bins:
                idx = np.flatnonzero(x > cfg.chisq_threshold)
                if len(idx) > cfg.max_chisq_points:
                    idx = idx[np.argsort(x[idx])[-cfg.max_chisq_points:]]
                if len(idx):
                    bins = power_chisq_bins(hp, cfg.chisq_bins, psds[ifo], cfg.f_low)
                    gidx = (idx + i0).astype(np.uint32)
                    chisq = power_chisq_at_points_from_precomputed(corr, su[gidx], norm, bins, gidx)
                    red = np.asarray(chisq) / (2 * cfg.chisq_bins - 2)
                    x[idx] = newsnr(x[idx], red).astype(np.float32)
                    n_chisq[ifo] += len(idx)
            p = _pool(x, n_pool)
            n_on = len(p) if n_on is None else min(n_on, len(p))
            pooled[ifo].append(p)
    P = {ifo: np.stack([p[:n_on] for p in pooled[ifo]]) for ifo in ifos}   # (n_templates, n_bins)

    def network(shifts):
        tot, n_above = None, None
        for k, ifo in enumerate(ifos):
            a = P[ifo] if shifts[k] == 0 else np.roll(P[ifo], shifts[k], axis=1)
            tot = a.astype(np.float64) ** 2 if tot is None else tot + a ** 2
            if cfg.coinc_threshold > 0:
                above = (a >= cfg.coinc_threshold).astype(np.int8)
                n_above = above if n_above is None else n_above + above
        net = np.sqrt(tot)
        if n_above is not None:
            net[n_above < 2] = 0.0           # coincidence: at least two detectors
        return net

    net = network([0] * len(ifos))
    it, ib = np.unravel_index(int(np.argmax(net)), net.shape)
    loud = float(net[it, ib])
    loud_gps = start + (ib + 0.5) * n_pool / fs
    rho_ifo = {ifo: float(P[ifo][it, ib]) for ifo in ifos}

    second = sorted(rho_ifo.values())[-2] if len(rho_ifo) > 1 else 0.0
    removed = cfg.hierarchical_removal and second >= cfg.removal_min_second_snr
    if removed:
        # a loud coincident candidate would otherwise enter every slide through
        # its single-detector SNRs (PyCBC "hierarchical removal"). Only for
        # candidates loud in >= 2 detectors: a single-detector glitch must stay
        # in the background, where its slid copies show it is not significant.
        w = max(int(round(cfg.removal_window_s * fs / n_pool)), 1)
        for ifo in ifos:
            P[ifo][:, max(ib - w, 0):ib + w + 1] = 0.0
    step = max(int(round(cfg.slide_step_s / (n_pool / fs))), 1)
    shifts = slide_shifts(len(ifos), n_on // step, cfg.n_slides)
    n_slides = len(shifts)
    bkg = [float(network([s * step for s in sh]).max()) for sh in shifts]
    p = targeted_pvalue(loud, bkg) if bkg else 1.0
    p_ext = extrapolated_pvalue(loud, bkg) if len(bkg) >= 50 else None
    reach = extrapolation_reach(loud, bkg) if bkg else None
    res = SubsegmentResult(start, end, ifos, loud, loud_gps, float(bank[it]), rho_ifo, bkg, p,
                            gates, f"{len(bank)} templates, {n_slides} slides, "
                            + ("no chi-square" if not cfg.chisq_bins else
                               f"newSNR ({cfg.chisq_bins}-bin chi-square at "
                               + ", ".join(f"{i}:{n_chisq[i]}" for i in ifos) + " points)")
                            + ", delays " + ", ".join(f"{i}:{delays[i] * 1e3:+.1f} ms" for i in ifos)
                            + (", loudest removed from background" if removed else ""))
    res.p_value_extrapolated = p_ext
    res.extrapolation_reach = reach
    if reach is not None and reach > 3:
        res.notes += (f"; loudest is {reach:.0f} background spreads beyond the loudest slide: "
                      "extrapolated p-value only formal, quote p < 1/(n_slides+1)")
    return res


def _antenna(ifos, ra_deg, dec_deg, gps):
    import lal
    gmst = lal.GreenwichMeanSiderealTime(lal.LIGOTimeGPS(float(gps)))
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    return {i: lal.ComputeDetAMResponse(lal.cached_detector_by_prefix[i].response, ra, dec, 0.0, gmst)
            for i in ifos}


def coherent_projector(sigmas: Sequence[float], fplus: Sequence[float],
                       fcross: Sequence[float]) -> np.ndarray:
    """Real projector onto span{sigma F+, sigma Fx} (detector space)."""
    w = np.stack([np.asarray(sigmas) * np.asarray(fplus), np.asarray(sigmas) * np.asarray(fcross)], 1)
    return w @ np.linalg.pinv(w.T @ w) @ w.T


def null_cut(rho_coh: np.ndarray, rho_null: np.ndarray, null_min: float = 5.25,
             null_grad: float = 0.2, null_step: float = 20.0) -> np.ndarray:
    """True where the candidate fails PyGRB's null-SNR cut."""
    thr = np.where(rho_coh < null_step, null_min, null_min + null_grad * (rho_coh - null_step))
    return rho_null > thr


def _coherent_subsegment(cfg, bank, ifos, stilde, psds, delays, gates, start, end):
    from pycbc.events.ranking import newsnr
    from pycbc.filter import matched_filter_core
    from pycbc.vetoes import power_chisq_at_points_from_precomputed, power_chisq_bins
    from pycbc.waveform import get_fd_waveform

    fs = cfg.sample_rate
    L = int(round((end - start) * fs))
    ant = _antenna(ifos, cfg.ra_deg, cfg.dec_deg, 0.5 * (start + end))
    fp = np.array([ant[i][0] for i in ifos])
    fx = np.array([ant[i][1] for i in ifos])
    step = max(int(round(cfg.slide_step_s * fs)), 1)
    shifts = [(0,) * len(ifos)] + slide_shifts(len(ifos), L // step, cfg.n_slides)
    sh = np.array(shifts, dtype=np.int64) * step                       # (n_slides+1, n_ifos) samples
    blk = max(int(round(cfg.removal_window_s * fs)), 1)
    n_blk = -(-L // blk)
    blockmax = np.zeros((len(shifts), n_blk))
    best = {"stat": 0.0}
    min_coinc = 2 if cfg.coinc_threshold > 0 else 1
    thr = cfg.coinc_threshold if cfg.coinc_threshold > 0 else cfg.chisq_threshold
    n_chisq = {i: 0 for i in ifos}
    n_null_cut = 0
    t0 = None
    for it, mc in enumerate(bank):
        m = _mass_equal(mc)
        z, sig = [], []
        for ifo in ifos:
            st = stilde[ifo]
            hp, _ = get_fd_waveform(approximant="TaylorF2", mass1=m, mass2=m,
                                    delta_f=st.delta_f, f_lower=cfg.f_low)
            hp.resize(len(st))
            t0 = float(st.start_time)
            x = (start + delays[ifo] - t0) * fs
            i0 = int(round(x))
            tau = (x - i0) / fs                     # fractional delay: shift the template
            hp = hp * np.exp(-2j * np.pi * hp.sample_frequencies.numpy() * tau)
            snr_u, corr, norm = matched_filter_core(hp, st, psd=psds[ifo],
                                                    low_frequency_cutoff=cfg.f_low)
            su = snr_u.numpy()
            zi = (su[i0:i0 + L] * norm).astype(np.complex64)
            if cfg.chisq_bins:
                a = np.abs(zi)
                idx = np.flatnonzero(a > cfg.chisq_threshold)
                if len(idx) > cfg.max_chisq_points:
                    idx = idx[np.argsort(a[idx])[-cfg.max_chisq_points:]]
                if len(idx):
                    bins = power_chisq_bins(hp, cfg.chisq_bins, psds[ifo], cfg.f_low)
                    gidx = (idx + i0).astype(np.uint32)
                    chisq = power_chisq_at_points_from_precomputed(corr, su[gidx], norm, bins, gidx)
                    red = np.asarray(chisq) / (2 * cfg.chisq_bins - 2)
                    zi[idx] *= (newsnr(a[idx], red) / a[idx]).astype(np.float32)
                    n_chisq[ifo] += len(idx)
            z.append(zi)
            sig.append(4 * st.delta_f / norm)       # sqrt(sigmasq): matched_filter_core's norm
        P = coherent_projector(sig, fp, fx)
        above = [np.flatnonzero(np.abs(zi) >= thr) for zi in z]
        for k in range(len(shifts)):
            pos = np.concatenate([(a + sh[k, d]) % L for d, a in enumerate(above)])
            if not len(pos):
                continue
            u, c = np.unique(pos, return_counts=True)
            cand = u[c >= min_coinc]
            if not len(cand):
                continue
            Z = np.stack([z[d][(cand - sh[k, d]) % L] for d in range(len(ifos))])   # (n_ifo, n_cand)
            net2 = np.sum(np.abs(Z) ** 2, axis=0)
            coh2 = np.real(np.einsum("ic,ij,jc->c", Z.conj(), P, Z))
            coh = np.sqrt(np.clip(coh2, 0, None))
            null = np.sqrt(np.clip(net2 - coh2, 0, None))
            stat = coh.copy()
            if len(ifos) >= 3:
                cut = null_cut(coh, null, cfg.null_min, cfg.null_grad, cfg.null_step)
                stat[cut] = 0.0
                if k == 0:
                    n_null_cut += int(cut.sum())
            np.maximum.at(blockmax[k], cand // blk, stat)
            if k == 0:
                j = int(np.argmax(stat))
                if stat[j] > best["stat"]:
                    best = {"stat": float(stat[j]), "it": it, "ib": int(cand[j]),
                            "coh": float(coh[j]), "null": float(null[j]),
                            "net": float(np.sqrt(net2[j])),
                            "rho_ifo": {i: float(abs(Z[d, j])) for d, i in enumerate(ifos)}}
    loud = best["stat"]
    if loud > 0:
        ib, rho_ifo = best["ib"], best["rho_ifo"]
        loud_gps, loud_mc = start + ib / fs, float(bank[best["it"]])
    else:
        ib, rho_ifo, loud_gps, loud_mc = 0, {}, start, 0.0
    second = sorted(rho_ifo.values())[-2] if len(rho_ifo) > 1 else 0.0
    removed = loud > 0 and cfg.hierarchical_removal and second >= cfg.removal_min_second_snr
    bkg = []
    for k in range(1, len(shifts)):
        b = blockmax[k].copy()
        if removed:
            for d in range(len(ifos)):
                c = ((ib + sh[k, d]) % L) // blk
                b[max(c - 1, 0):c + 2] = 0.0
        bkg.append(float(b.max()))
    p = targeted_pvalue(loud, bkg) if bkg else 1.0
    p_ext = extrapolated_pvalue(loud, bkg) if len(bkg) >= 50 else None
    reach = extrapolation_reach(loud, bkg) if bkg else None
    res = SubsegmentResult(start, end, ifos, loud, loud_gps, loud_mc, rho_ifo, bkg, p, gates,
                           f"coherent, {len(bank)} templates, {len(bkg)} slides, "
                           + ("no chi-square" if not cfg.chisq_bins else
                              f"newSNR ({cfg.chisq_bins}-bin chi-square at "
                              + ", ".join(f"{i}:{n_chisq[i]}" for i in ifos) + " points)")
                           + (f", {n_null_cut} zero-lag samples failed the null cut"
                              if len(ifos) >= 3 else ", 2 detectors: no null stream")
                           + ", delays " + ", ".join(f"{i}:{delays[i] * 1e3:+.1f} ms" for i in ifos)
                           + (", loudest removed from background" if removed else ""),
                           statistic="coherent")
    if loud > 0:
        res.loudest_rho_coh, res.loudest_rho_null = best["coh"], best["null"]
        res.loudest_rho_network = best["net"]
    res.p_value_extrapolated = p_ext
    res.extrapolation_reach = reach
    if reach is not None and reach > 3:
        res.notes += (f"; loudest is {reach:.0f} background spreads beyond the loudest slide: "
                      "extrapolated p-value only formal, quote p < 1/(n_slides+1)")
    return res


def offsource_windows(cfg: TargetedSearchConfig, length: float, n: int) -> list[tuple[float, float]]:
    """``n`` windows of ``length`` s outside the whole on-source window,
    alternately before and after it, ``offsource_gap_s`` away from it and
    from each other's neighbours (PyGRB-style off-source trials)."""
    out, k_before, k_after = [], 0, 0
    for i in range(n):
        if i % 2 == 0:
            b = cfg.gps_start - cfg.offsource_gap_s - k_before * length
            out.append((b - length, b))
            k_before += 1
        else:
            a = cfg.gps_end + cfg.offsource_gap_s + k_after * length
            out.append((a, a + length))
            k_after += 1
    return out


def offsource_background(cfg: TargetedSearchConfig, bank: Sequence[float], fetch: StrainFetcher,
                         length: float, n: int) -> list[float]:
    """Loudest zero-lag network statistic in each off-source trial of the
    same length and sky position as an on-source sub-segment. Unlike time
    slides, this samples the real non-stationarity of the data around the
    window (at the cost of a full matched filter per trial). Trials with
    fewer than two detectors are skipped."""
    c = replace(cfg, n_slides=0)
    out = []
    for a, b in offsource_windows(cfg, length, n):
        r = search_subsegment(c, bank, fetch, a, b)
        if len(r.ifos) >= 2:
            out.append(r.loudest_rho_net)
    return out


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
        if cfg.n_offsource > 0 and len(r.ifos) >= 2:
            r.offsource = offsource_background(cfg, bank, fetch, b - a, cfg.n_offsource)
            if r.offsource:
                r.p_value_offsource = targeted_pvalue(r.loudest_rho_net, r.offsource)
                r.notes += f"; {len(r.offsource)} off-source trials"
        subs.append(r)
        if progress:
            print(f"[{a:.0f}, {b:.0f}] rho_net={r.loudest_rho_net:.2f} at {r.loudest_gps:.2f} "
                  f"Mc={r.loudest_mchirp:.4f} p={r.p_value:.4f}"
                  + (f" p_off={r.p_value_offsource:.3f}" if r.p_value_offsource is not None else "")
                  + f" {r.notes}", flush=True)
    valid = [s for s in subs if s.background]
    # with off-source trials, the more conservative of the two p-values per sub-segment
    p_all = combine_subsegments([max(s.p_value, s.p_value_offsource or 0.0) for s in valid]) \
        if valid else 1.0
    best = max(subs, key=lambda s: s.loudest_rho_net) if subs else None
    loudest = asdict(best) if best else {}
    loudest.pop("background", None)
    loudest.pop("offsource", None)
    res = TargetedSearchResult(asdict(cfg), len(bank), subs, p_all, loudest)
    if out_dir is not None:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        res.to_json(out / f"{cfg.label}_targeted_result.json")
    return res
