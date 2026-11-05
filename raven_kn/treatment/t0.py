"""Merger-time (t0) window: the coincidence window [tl, th] of RAVEN.

For GRBs RAVEN uses a fixed window of seconds. For a kilonova the window
comes from the data, in decreasing order of precision:

1. a GRB (or other prompt signal) associated with the kilonova: seconds;
2. a t0 posterior from a light-curve fit (NMMA, Redback, ...): Jamsin et al.
   (arXiv:2610.28775) find 0.3--0.7 d (90 %) at best with LSST, ~5.5 d
   typically;
3. a bracket between the last constraining non-detection and the first
   detection (no model fit; the fallback, and a sanity check on 2).

Windows from several methods are intersected.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

from ..connection import KNCandidate

DAY = 86400.0


@dataclass
class MergerTimeWindow:
    """GPS interval in which the merger is believed to lie."""
    gps_start: float
    gps_end: float
    method: str
    credible_level: Optional[float] = None
    notes: str = ""
    # data that anchor the merger time: last constraining non-detection
    # (early side) and first detection (late side), GPS
    anchor_early_gps: Optional[float] = None
    anchor_late_gps: Optional[float] = None

    @property
    def duration(self) -> float:
        return self.gps_end - self.gps_start

    def relative_to(self, gps_ref: float) -> tuple[float, float]:
        """(tl, th) in seconds relative to ``gps_ref``, RAVEN convention."""
        return self.gps_start - gps_ref, self.gps_end - gps_ref

    def intersect(self, other: "MergerTimeWindow") -> "MergerTimeWindow":
        lo, hi = max(self.gps_start, other.gps_start), min(self.gps_end, other.gps_end)
        if lo >= hi:
            raise ValueError(f"Inconsistent windows: {self} and {other}")
        return MergerTimeWindow(lo, hi, f"{self.method}&{other.method}",
                                notes="intersection")


# GRB: RAVEN searches GW candidates from 5 s before to 1 s after a GRB
# (merger precedes the burst; GRB 170817A came 1.74 s after merger).
GRB_WINDOW = (-5.0, 1.0)


def window_from_grb(grb_gps: float, tl: float = GRB_WINDOW[0],
                    th: float = GRB_WINDOW[1]) -> MergerTimeWindow:
    return MergerTimeWindow(grb_gps + tl, grb_gps + th, "GRB", 1.0)


def window_from_samples(t0_gps_samples: Sequence[float],
                        credible_level: float = 0.9,
                        pad_s: float = 0.0,
                        model_sys_frac: float = 0.0,
                        model_sys_min_h: float = 1.0,
                        anchor_early_gps: Optional[float] = None,
                        anchor_late_gps: Optional[float] = None) -> MergerTimeWindow:
    """Central credible interval of a t0 posterior (light-curve fit), widened
    by a model-systematic margin.

    The statistical interval assumes the kilonova model is right. Models
    differ most at early times, so the error on t0 grows with how far the
    fit extrapolates from the data that anchor it:

    * late side:  ``model_sys_frac * (first detection - t0_median)``
    * early side: ``model_sys_frac * (t0_median - last constraining
      non-detection)``, or the late-side gap if there is no such limit,

    each at least ``model_sys_min_h``. Calibration: with Bu2019lm on the real
    AT2017gfo light curve (first point +17 h, no pre-merger limit) the
    median t0 is 40 h early and the 90 % interval misses the merger; on
    GRB 211211A (first point +6 h) it is 28 h early. 0.75 is the smallest
    round fraction recovering both (0.5 misses GRB 211211A by 0.5 h). Two
    events, one model: to be revisited with cross-model injections (Jamsin
    et al., arXiv:2610.28775).
    """
    s = np.asarray(t0_gps_samples, float)
    a = (1 - credible_level) / 2
    lo, med, hi = np.quantile(s, [a, 0.5, 1 - a])
    early = late = 0.0
    notes = ""
    if model_sys_frac > 0 and anchor_late_gps is not None:
        gap_late = max(anchor_late_gps - med, 0.0)
        gap_early = (max(med - anchor_early_gps, 0.0) if anchor_early_gps is not None
                     else gap_late)
        floor = model_sys_min_h * 3600
        late = max(model_sys_frac * gap_late, floor)
        early = max(model_sys_frac * gap_early, floor)
        notes = (f"model-systematic margin -{early / 3600:.1f} h / +{late / 3600:.1f} h "
                 f"(frac {model_sys_frac}, stat. width {(hi - lo) / 3600:.1f} h)")
    return MergerTimeWindow(lo - pad_s - early, hi + pad_s + late, "lc_fit", credible_level,
                            notes, anchor_early_gps, anchor_late_gps)


def distance_modulus(d_mpc: float | np.ndarray) -> float | np.ndarray:
    return 5 * np.log10(np.asarray(d_mpc) * 1e6 / 10)


def window_from_bracket(kn: KNCandidate,
                        d_mpc_upper: Optional[float] = None,
                        early_abs_mag: float = -15.0,
                        t_visible_d: float = 0.25,
                        max_age_d: float = 7.0) -> MergerTimeWindow:
    """Bracket the merger between the last *constraining* non-detection and
    the first detection.

    A non-detection at t_UL with limiting magnitude m_lim constrains the
    merger only if a kilonova at the candidate's distance would have been
    brighter than m_lim at t_UL, had the merger already happened. Model
    (deliberately simple, to be calibrated on Jamsin et al. injections):
    from ``t_visible_d`` after merger until the first detection, the
    kilonova is at least as bright as ``early_abs_mag`` (AT2017gfo was at
    M ~ -16 in the optical ~11 h after merger). A limit is constraining if
    ``m_lim > early_abs_mag + DM(d_mpc_upper)``; then t0 > t_UL - t_visible.

    Without a constraining limit the early side is set to ``max_age_d``
    before the first detection (worst case in Jamsin et al.: about a week).

    ``d_mpc_upper`` should be a conservative (e.g. 95th percentile)
    distance; if unknown, any limit deeper than the first detection
    magnitude is taken as constraining.
    """
    fd = kn.first_detection
    if fd is None:
        raise ValueError(f"{kn.name}: no detection to bracket")
    t_det = fd.gps
    lo = t_det - max_age_d * DAY
    used = None
    for ul in reversed(kn.upper_limits):          # latest first
        if ul.gps >= t_det:
            continue
        if d_mpc_upper is not None:
            needed = early_abs_mag + distance_modulus(d_mpc_upper)
        else:
            needed = fd.mag
        if ul.mag > needed:
            cand = ul.gps - t_visible_d * DAY
            if cand > lo:
                lo, used = cand, ul
            break
    notes = (f"last constraining non-detection {used.survey} {used.band}<{used.mag:.1f} "
             f"at MJD {used.mjd:.4f}" if used else
             f"no constraining non-detection: early edge = first detection - {max_age_d} d")
    return MergerTimeWindow(lo, t_det, "bracket", None, notes,
                            used.gps if used else None, t_det)


def estimate_window(kn: KNCandidate, d_mpc_upper: Optional[float] = None,
                    credible_level: float = 0.9, model_sys_frac: float = 0.75,
                    model_sys_min_h: float = 1.0, **bracket_kw) -> MergerTimeWindow:
    """Best available window: GRB > (light-curve fit + model margin) & bracket > bracket."""
    if kn.grb is not None:
        return window_from_grb(kn.grb.gps)
    win = None
    if kn.first_detection is not None:
        win = window_from_bracket(kn, d_mpc_upper=d_mpc_upper, **bracket_kw)
    if kn.t0_samples_gps:
        fit = window_from_samples(
            kn.t0_samples_gps, credible_level, model_sys_frac=model_sys_frac,
            model_sys_min_h=model_sys_min_h,
            anchor_early_gps=win.anchor_early_gps if win else None,
            anchor_late_gps=win.anchor_late_gps if win else None)
        win = fit if win is None else _safe_intersect(fit, win)
    if win is None:
        raise ValueError(f"{kn.name}: no GRB, photometry or t0 samples")
    return win


def _safe_intersect(fit: MergerTimeWindow, bracket: MergerTimeWindow) -> MergerTimeWindow:
    try:
        w = fit.intersect(bracket)
        w.notes = "; ".join(x for x in (fit.notes, bracket.notes) if x)
        w.anchor_early_gps, w.anchor_late_gps = bracket.anchor_early_gps, bracket.anchor_late_gps
        return w
    except ValueError:
        # The fit disagrees with the data bracket: keep the bracket, flag it
        # (a contaminant often gives an implausible t0, cf. AT2025ulz).
        bracket.notes += "; WARNING light-curve t0 outside the photometric bracket"
        return bracket
