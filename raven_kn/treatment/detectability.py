"""Is a GW search worth it? Expected network SNR of the kilonova's merger.

For each (D_L, iota) sample of the EM posterior, a merger time is drawn in
the t0 window and a polarisation angle uniformly; the SNR in each detector is

    rho_det = 8 * (D_h,det / D_L) * sqrt(F+^2 ((1 + cos^2 iota)/2)^2 + Fx^2 cos^2 iota)

with D_h the optimal horizon (overhead, face-on, SNR 8) and D_h = 2.264 x
the BNS range (sky/orientation-averaged, 1.4+1.4 Msun) quoted in detector
summary pages. F+, Fx come from LAL at the candidate's position and the
drawn time (Earth rotation over a multi-day window matters). Without LAL a
sky-averaged F^2 = 2/5 is used.

The probability that the network SNR exceeds the targeted-search threshold
is the figure of merit. A hard distance cut (default 500 Mpc) can be added:
beyond it the O4/O5 BNS sensitivity makes the search pointless.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

from .distance_inclination import DistanceInclinationPosterior
from .t0 import MergerTimeWindow

HORIZON_OVER_RANGE = 2.264

# Approximate O4 BNS ranges (Mpc). Override with the live values of the
# detector summary pages for the epoch of the window.
DEFAULT_RANGES_MPC = {"H1": 150.0, "L1": 160.0, "V1": 55.0}

# Chirp mass of the 1.4 + 1.4 Msun binary that quoted BNS ranges refer to.
BNS_CHIRP_MASS = 1.2188


def scale_ranges(ranges_mpc: dict[str, float], chirp_mass_msun: Optional[float]) -> dict[str, float]:
    """Ranges for another chirp mass: inspiral amplitude scales as Mc^(5/6).

    Valid for inspiral-dominated signals (sub-solar to a few Msun chirp
    mass), e.g. the chirp-mass bin of an LVK alert (HasSSM candidates:
    0.1-0.87 Msun). Not for heavy binary black holes, whose merger lies in
    band. None: ranges unchanged (BNS).
    """
    if chirp_mass_msun is None:
        return dict(ranges_mpc)
    if chirp_mass_msun <= 0:
        raise ValueError("chirp mass must be positive")
    f = (chirp_mass_msun / BNS_CHIRP_MASS) ** (5 / 6)
    return {k: v * f for k, v in ranges_mpc.items()}


def _antenna_grid(ifos: Sequence[str], ra_deg: float, dec_deg: float,
                  times: np.ndarray) -> Optional[dict[str, np.ndarray]]:
    """F+ and Fx at psi = 0 for each detector and time: shape (n_t, 2)."""
    try:
        import lal
    except ImportError:
        return None
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    out = {}
    for ifo in ifos:
        resp = lal.cached_detector_by_prefix[ifo].response
        f = np.empty((len(times), 2))
        for i, t in enumerate(times):
            gmst = lal.GreenwichMeanSiderealTime(lal.LIGOTimeGPS(float(t)))
            f[i] = lal.ComputeDetAMResponse(resp, ra, dec, 0.0, gmst)
        out[ifo] = f
    return out


def _in_segments(t: np.ndarray, segs: Sequence[Sequence[float]]) -> np.ndarray:
    on = np.zeros(t.shape, bool)
    for a, b in segs:
        on |= (t >= a) & (t < b)
    return on


@dataclass
class Detectability:
    p_detect: float
    snr_median: float
    snr_90: tuple[float, float]
    snr_threshold: float
    ranges_mpc: dict[str, float]
    data_fraction: dict[str, float]
    antenna: str
    worth_search: bool
    reasons: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {k: (list(v) if isinstance(v, tuple) else v)
                for k, v in self.__dict__.items()}


def network_snr_samples(post: DistanceInclinationPosterior, window: MergerTimeWindow,
                        ra_deg: float, dec_deg: float,
                        ranges_mpc: dict[str, float] = DEFAULT_RANGES_MPC,
                        segments: Optional[dict[str, Sequence[Sequence[float]]]] = None,
                        duty_cycle: Optional[dict[str, float]] = None,
                        n_times: int = 96, seed: Optional[int] = 1):
    """Network SNR for each posterior sample (and per-detector data fraction)."""
    rng = np.random.default_rng(seed)
    n = len(post.d_mpc)
    ifos = list(ranges_mpc)
    tgrid = np.linspace(window.gps_start, window.gps_end, n_times)
    grid = _antenna_grid(ifos, ra_deg, dec_deg, tgrid)
    it = rng.integers(0, n_times, n)
    t = tgrid[it]
    psi = rng.uniform(0, np.pi, n)
    c2p, s2p = np.cos(2 * psi), np.sin(2 * psi)
    ci = post.cos_iota
    ap, ax = (1 + ci ** 2) / 2, ci
    snr2 = np.zeros(n)
    frac = {}
    for ifo in ifos:
        dh = HORIZON_OVER_RANGE * ranges_mpc[ifo]
        if grid is not None:
            fp0, fx0 = grid[ifo][it, 0], grid[ifo][it, 1]
            fp, fx = fp0 * c2p + fx0 * s2p, -fp0 * s2p + fx0 * c2p
            amp2 = (fp * ap) ** 2 + (fx * ax) ** 2
        else:
            amp2 = 0.2 * (ap ** 2 + ax ** 2)          # <F+^2> = <Fx^2> = 1/5
        on = np.ones(n, bool)
        if segments and ifo in segments:
            on = _in_segments(t, segments[ifo])
        elif duty_cycle and ifo in duty_cycle:
            on = rng.random(n) < duty_cycle[ifo]
        frac[ifo] = float(on.mean())
        snr2 += on * (8 * dh / post.d_mpc) ** 2 * amp2
    return np.sqrt(snr2), frac, ("LAL" if grid is not None else "sky-averaged")


def assess_detectability(post: DistanceInclinationPosterior, window: MergerTimeWindow,
                         ra_deg: float, dec_deg: float,
                         ranges_mpc: dict[str, float] = DEFAULT_RANGES_MPC,
                         snr_threshold: float = 6.0, p_min: float = 0.1,
                         d_max_useful_mpc: Optional[float] = 500.0,
                         segments=None, duty_cycle=None, seed: Optional[int] = 1,
                         chirp_mass_msun: Optional[float] = None) -> Detectability:
    """``ranges_mpc`` are BNS (1.4 + 1.4 Msun) ranges; ``chirp_mass_msun``
    rescales them to the source class searched for (see :func:`scale_ranges`)."""
    ranges_mpc = scale_ranges(ranges_mpc, chirp_mass_msun)
    snr, frac, antenna = network_snr_samples(post, window, ra_deg, dec_deg, ranges_mpc,
                                             segments, duty_cycle, seed=seed)
    p = float(np.mean(snr >= snr_threshold))
    d_med = float(np.median(post.d_mpc))
    reasons, worth = [], True
    if d_max_useful_mpc is not None and d_med > d_max_useful_mpc:
        worth = False
        reasons.append(f"median D_L {d_med:.0f} Mpc > {d_max_useful_mpc:.0f} Mpc cut")
    if p < p_min:
        worth = False
        reasons.append(f"P(network SNR >= {snr_threshold}) = {p:.3f} < {p_min}")
    if all(f == 0 for f in frac.values()):
        worth = False
        reasons.append("no detector data in the window")
    if post.constraints == []:
        reasons.append("WARNING: no distance information; prior only")
    if worth:
        reasons.append(f"P(network SNR >= {snr_threshold}) = {p:.3f}")
    if chirp_mass_msun is not None:
        reasons.append(f"ranges scaled to chirp mass {chirp_mass_msun:g} Msun")
    return Detectability(p, float(np.median(snr)),
                         tuple(float(x) for x in np.quantile(snr, [0.05, 0.95])),
                         snr_threshold, dict(ranges_mpc), frac, antenna, worth, reasons)
