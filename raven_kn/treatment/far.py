"""Joint false alarm rates, with the same conventions as ``ligo.raven.search.coinc_far``.

Two situations:

* **Coincidence with an existing GW candidate** (a superevent from the
  online pipelines, possibly low-significance, inside the t0 window):
  RAVEN's temporal FAR ``(th - tl) * R_ext * FAR_gw``, divided by the sky
  overlap integral and, new here, by a distance overlap factor (the
  kilonova distance is known, the GW sky map gives a distance per
  direction). For kilonovae ``th - tl`` is the t0 window (hours to days)
  instead of seconds, which is why a narrow window matters.
* **Targeted search result** (PyCBC/PyGRB-like on-source/off-source): the
  p-value of the loudest on-source trigger from the off-source background,
  then converted into a joint FAR with RAVEN's "SubGRBTargeted" product
  formula, the kilonova classifier FAR playing the GRB FAR's role.

``kn_rate`` is the rate of *candidates* passed by the kilonova filter (true
kilonovae plus contaminants), not the astrophysical rate. Default: ~30
candidates per month (order of magnitude quoted for current broker
filters; set it from the filter actually used).
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

SECONDS_PER_YEAR = 365.25 * 86400
DEFAULT_KN_CANDIDATE_RATE = 30 * 12 / SECONDS_PER_YEAR      # Hz
# Transients spectroscopically classified as kilonovae (TNS "Kilonova" or a
# circular): a handful per year at most (AT2017gfo, AT2025ulz, a few GRB
# kilonovae). 5 / yr is a conservative upper value.
SPEC_KN_RATE = 5 / SECONDS_PER_YEAR                          # Hz
SPEC_KN_CLASSES = ("kilonova", "kn")


def candidate_rate(kn, base_rate: float = DEFAULT_KN_CANDIDATE_RATE,
                   spec_kn_rate: float = SPEC_KN_RATE,
                   score_survival: Optional[Sequence[tuple[float, float]]] = None) -> tuple[float, str]:
    """Rate R_KN of external triggers at least as kilonova-like as ``kn``,
    the analogue of RAVEN's GRB rate in the temporal FAR.

    * spectroscopically classified kilonova (``kn.classification``):
      ``spec_kn_rate``;
    * else, with a classifier score and the filter's score distribution
      ``score_survival`` = [(score, fraction of filter candidates scoring
      >= score), ...]: ``base_rate`` x that fraction (interpolated; this is
      the ``far_kn`` of :func:`targeted_far`);
    * else ``base_rate``, the rate of everything the filter passes.

    Returns (rate in Hz, explanation).
    """
    cls = (getattr(kn, "classification", None) or "").strip().lower()
    if cls in SPEC_KN_CLASSES:
        return spec_kn_rate, f"spectroscopically classified kilonova: {spec_kn_rate * SECONDS_PER_YEAR:g}/yr"
    score = getattr(kn, "kn_score", None)
    if score is not None and score_survival:
        sc, fr = np.array(sorted(score_survival)).T
        frac = float(np.clip(np.interp(score, sc, fr), 0.0, 1.0))
        frac = max(frac, 1.0 / (base_rate * SECONDS_PER_YEAR))     # at least one candidate per year
        return base_rate * frac, (f"filter rate x fraction scoring >= {score:g} ({frac:.3g}): "
                                  f"{base_rate * frac * SECONDS_PER_YEAR:.3g}/yr")
    return base_rate, f"filter candidate rate {base_rate * SECONDS_PER_YEAR:.3g}/yr"


def coinc_far_kn(window_s: float, far_gw: float, kn_rate: float = DEFAULT_KN_CANDIDATE_RATE,
                 sky_overlap: Optional[float] = None,
                 distance_overlap: Optional[float] = None) -> dict:
    """RAVEN-style joint FAR for a GW candidate inside a kilonova t0 window.

    Returns the same keys as ``ligo.raven.search.coinc_far`` plus the
    distance-informed one.
    """
    if window_s <= 0:
        raise ValueError("window must be positive")
    temporal = window_s * kn_rate * far_gw
    spatiotemporal = temporal / sky_overlap if sky_overlap else None
    full = None
    if spatiotemporal is not None and distance_overlap:
        full = spatiotemporal / distance_overlap
    return {"temporal_coinc_far": temporal,
            "spatiotemporal_coinc_far": spatiotemporal,
            "spatiotemporal_distance_coinc_far": full,
            "skymap_overlap": sky_overlap,
            "distance_overlap": distance_overlap,
            "window_s": window_s, "kn_rate": kn_rate, "far_gw": far_gw}


def targeted_pvalue(loudest_onsource: float, offsource_loudest: Sequence[float]) -> float:
    """p-value of the loudest on-source statistic against the loudest
    statistic of each off-source trial (time slides / off-source segments
    of the same length). Conservative (n_louder + 1) / (n + 1)."""
    off = np.asarray(offsource_loudest)
    return float((np.sum(off >= loudest_onsource) + 1) / (len(off) + 1))


def combine_subsegments(pvalues: Sequence[float]) -> float:
    """Long windows are split in sub-segments (Earth rotation, Jamsin et al.);
    trials factor for the minimum p-value over independent sub-segments."""
    p = np.asarray(pvalues)
    return float(1 - (1 - p.min()) ** len(p))


def targeted_far(window_s: float, far_kn: float, far_gw: float,
                 far_kn_thresh: float = DEFAULT_KN_CANDIDATE_RATE,
                 far_gw_thresh: float = 2 / 86400.) -> float:
    """RAVEN 'SubGRBTargeted' formula: the product of two uniform FARs mapped
    to its CDF, z (1 - ln(z / z_max)), with z = T far_kn far_gw.

    ``far_kn``: rate of filter candidates at least as kilonova-like as this
    one (from the classifier score distribution). ``far_gw``: FAR of the
    targeted-search trigger (e.g. p-value / window).
    """
    z = window_s * far_kn * far_gw
    z_max = window_s * far_kn_thresh * far_gw_thresh
    z = min(z, z_max)
    return float(z * (1 - np.log(z / z_max)))
