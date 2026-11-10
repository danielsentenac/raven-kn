"""Is the candidate a kilonova? Evidence required before a GW search is triggered.

A GRB or a single optical detection is not enough to spend a targeted GW
search: most optical transients in a GRB error box are unrelated supernovae
or AGN flares. The search is triggered only once at least one of these holds:

* **spectrum**: a spectroscopic classification containing "kilonova"
  (TNS, GCN circular);
* **light curve**, from >= 2 detection epochs at least ``min_epoch_sep_d``
  apart, with a distance for the candidate:

  - peak absolute magnitude inside ``abs_mag_range`` (kilonovae: about
    -17 to -13, widened by 0.5 mag each side for distance and real-time
    photometry errors: the first GCN magnitude of AT2017gfo, i = 16.0, was
    1 mag brighter than the published one; supernovae peak at -16 to -20
    days later, fast blue optical transients at -19 to -21), and
  - fast evolution in the first ``max_age_d`` days after discovery: fading
    of at least ``min_decline`` mag/day in one band, or reddening (blue minus
    red colour increasing) of at least ``min_reddening`` mag/day, each at
    >= ``min_sigma`` (errors floored at 0.03 mag for cross-calibration).
    Normal supernovae rise or fade at < 0.1 mag/day in the first days.

* **classifier score**: an upstream kilonova classifier score ``kn_score``
  >= ``min_score`` (e.g. a Fink or ALeRCE filter), if given.

Magnitudes are not corrected for Galactic extinction (a small shift of the
absolute magnitude, no effect on rates of change).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Optional

import numpy as np

# effective wavelength (micron) by band name, for colours
_WAVELENGTH = {"u": 0.36, "B": 0.44, "g": 0.48, "c": 0.53, "V": 0.55, "r": 0.62, "R": 0.64,
               "o": 0.68, "w": 0.61, "i": 0.75, "I": 0.80, "z": 0.87, "y": 0.96, "Y": 1.03,
               "J": 1.25, "H": 1.65, "K": 2.2, "Ks": 2.15}


def band_wavelength(band: str) -> Optional[float]:
    """Effective wavelength of a band name such as 'ps1::g', 'ztfr', 'sdssu',
    'lsstg', 'r', 'Ks' (case matters only for J, H, K, B, V, R, I, Y)."""
    b = band.split("::")[-1]
    for pre in ("2mass", "ztf", "lsst", "sdss", "ps1", "atlas", "decam", "swope"):
        if b.lower().startswith(pre) and len(b) > len(pre):
            b = b[len(pre):]
            break
    if b in _WAVELENGTH:
        return _WAVELENGTH[b]
    if b.lower() in _WAVELENGTH:
        return _WAVELENGTH[b.lower()]
    return None


@dataclass
class KilonovaEvidence:
    confirmed: bool
    by: list[str] = field(default_factory=list)       # "spectrum", "light curve", "classifier score"
    reasons: list[str] = field(default_factory=list)
    n_epochs: int = 0
    peak_abs_mag: Optional[float] = None
    max_decline_mag_per_day: Optional[float] = None
    max_reddening_mag_per_day: Optional[float] = None

    def to_dict(self) -> dict:
        return asdict(self)


def _slope(t: np.ndarray, y: np.ndarray, err: np.ndarray) -> tuple[float, float]:
    """Weighted least-squares slope and its 1-sigma error."""
    w = 1 / err ** 2
    tm = np.sum(w * t) / np.sum(w)
    sxx = np.sum(w * (t - tm) ** 2)
    slope = float(np.sum(w * (t - tm) * y) / sxx)
    return slope, float(1 / np.sqrt(sxx))


def _err(p, floor: float) -> float:
    e = p.mag_err if p.mag_err is not None and np.isfinite(p.mag_err) else 0.1
    return max(float(e), floor)


def _epochs(mjd: np.ndarray, sep_d: float) -> int:
    n, last = 0, None
    for t in np.sort(mjd):
        if last is None or t - last >= sep_d:
            n, last = n + 1, t
    return n


def _decline_rates(dets, t_first: float, max_age_d: float, min_dt_d: float,
                   err_floor: float = 0.03) -> dict[str, tuple[float, float]]:
    """Per band: fading rate (mag/day, > 0 = fading) and its error, from the
    band's brightest point to its later points, by weighted least squares."""
    out = {}
    bands = {p.band for p in dets}
    for b in bands:
        pts = sorted((p for p in dets if p.band == b and p.mjd - t_first <= max_age_d),
                     key=lambda p: p.mjd)
        if len(pts) < 2:
            continue
        k = int(np.argmin([p.mag for p in pts]))
        after = pts[k:]
        t = np.array([p.mjd for p in after])
        if len(after) < 2 or t[-1] - t[0] < min_dt_d:
            continue
        out[b] = _slope(t, np.array([p.mag for p in after]),
                        np.array([_err(p, err_floor) for p in after]))
    return out


def _reddening_rate(dets, t_first: float, max_age_d: float, pair_dt_d: float,
                    min_dt_d: float, err_floor: float = 0.03) -> Optional[tuple[float, float]]:
    """Rate of change (mag/day) of the colour bluer - redder band, and its
    error, from near-simultaneous pairs (within ``pair_dt_d``) of any two
    bands at least 0.1 micron apart; > 0 = reddening. The pair with the most
    significant reddening."""
    pts = [p for p in dets if p.mjd - t_first <= max_age_d and band_wavelength(p.band)]
    by_band: dict[str, list] = {}
    for p in pts:
        by_band.setdefault(p.band, []).append(p)
    best = None
    names = sorted(by_band, key=band_wavelength)
    for i, blue in enumerate(names):
        for red in names[i + 1:]:
            if band_wavelength(red) - band_wavelength(blue) < 0.1:
                continue
            col = []
            for pb in by_band[blue]:
                near = [pr for pr in by_band[red] if abs(pr.mjd - pb.mjd) <= pair_dt_d]
                if near:
                    pr = min(near, key=lambda q: abs(q.mjd - pb.mjd))
                    col.append((0.5 * (pb.mjd + pr.mjd), pb.mag - pr.mag,
                                np.hypot(_err(pb, err_floor), _err(pr, err_floor))))
            if len(col) < 2:
                continue
            t, c, e = np.array(col).T
            if t.max() - t.min() < min_dt_d:
                continue
            rate = _slope(t, c, e)
            if best is None or rate[0] / rate[1] > best[0] / best[1]:
                best = rate
    return best


def kilonova_evidence(kn, d_mpc: Optional[float] = None,
                      abs_mag_range: tuple[float, float] = (-17.5, -12.5),
                      min_decline: float = 0.3, min_reddening: float = 0.3,
                      max_age_d: float = 4.0, min_epoch_sep_d: float = 0.25,
                      pair_dt_d: float = 0.2, min_sigma: float = 3.0,
                      min_score: Optional[float] = 0.8) -> KilonovaEvidence:
    """Evidence that ``kn`` is a kilonova (see the module docstring).
    ``d_mpc``: luminosity distance (median of the EM posterior)."""
    ev = KilonovaEvidence(confirmed=False)
    cls = (getattr(kn, "classification", None) or "").lower()
    if "kilonova" in cls:
        ev.by.append("spectrum")
        ev.reasons.append(f"spectroscopic classification '{kn.classification}'")
    score = getattr(kn, "kn_score", None)
    if min_score is not None and score is not None:
        if score >= min_score:
            ev.by.append("classifier score")
        ev.reasons.append(f"classifier score {score:.2f} ({'>=' if score >= min_score else '<'} "
                          f"{min_score})")

    dets = kn.detections
    if not dets:
        ev.reasons.append("no detection")
        ev.confirmed = bool(ev.by)
        return ev
    t_first = dets[0].mjd
    ev.n_epochs = _epochs(np.array([p.mjd for p in dets]), min_epoch_sep_d)
    rates = _decline_rates(dets, t_first, max_age_d, min_epoch_sep_d)
    # a rate counts only if significant: rate >= threshold and rate / error >= min_sigma
    fade = [(r, e, b) for b, (r, e) in rates.items() if r / e >= min_sigma]
    if fade:
        r, e, b = max(fade)
        ev.max_decline_mag_per_day, decline_txt = r, f"fading {r:.2f} +- {e:.2f} mag/d in {b}"
    else:
        decline_txt = None
    red = _reddening_rate(dets, t_first, max_age_d, pair_dt_d, min_epoch_sep_d)
    if red is not None and red[0] / red[1] >= min_sigma:
        ev.max_reddening_mag_per_day = red[0]
        red_txt = f"reddening {red[0]:.2f} +- {red[1]:.2f} mag/d"
    else:
        red_txt = None
    if d_mpc is not None and d_mpc > 0:
        dm = 5 * np.log10(d_mpc * 1e6 / 10)
        ev.peak_abs_mag = float(min(p.mag for p in dets) - dm)

    lc_ok, why = True, []
    if ev.n_epochs < 2:
        lc_ok = False
        why.append(f"{ev.n_epochs} detection epoch(s), need 2 at least {min_epoch_sep_d} d apart")
    if ev.peak_abs_mag is None:
        lc_ok = False
        why.append("no distance: absolute magnitude unknown")
    elif not abs_mag_range[0] <= ev.peak_abs_mag <= abs_mag_range[1]:
        lc_ok = False
        why.append(f"peak M = {ev.peak_abs_mag:.1f} outside the kilonova range "
                   f"[{abs_mag_range[0]}, {abs_mag_range[1]}]")
    else:
        why.append(f"peak M = {ev.peak_abs_mag:.1f} in the kilonova range")
    fast = []
    if ev.max_decline_mag_per_day is not None and ev.max_decline_mag_per_day >= min_decline:
        fast.append(decline_txt)
    if ev.max_reddening_mag_per_day is not None and ev.max_reddening_mag_per_day >= min_reddening:
        fast.append(red_txt)
    if fast:
        why.append(" and ".join(fast) + f" within {max_age_d:g} d")
    elif ev.n_epochs >= 2:
        lc_ok = False
        why.append(f"no significant fast evolution ({decline_txt or 'no fading rate'}, "
                   f"{red_txt or 'no reddening rate'} at >= {min_sigma:g} sigma; need >= "
                   f"{min_decline} / {min_reddening} mag/d)")
    if lc_ok:
        ev.by.append("light curve")
    ev.reasons.append("light curve: " + "; ".join(why))
    ev.confirmed = bool(ev.by)
    return ev
