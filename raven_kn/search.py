"""RAVEN-style orchestration for kilonova triggers.

Same entry points as ``ligo.raven.search`` (``query``, ``search``,
``coinc_far``, ``calc_signif_gracedb``) with the kilonova t0 window in place
of RAVEN's fixed [tl, th], plus :func:`assess`, which decides whether a
targeted GW search is worth launching and writes its inputs.

GraceDB messages use the prefix ``RAVEN-KN:`` and the ``ext_coinc`` tag, as
RAVEN does with ``RAVEN:``.
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Optional

import numpy as np

from .connection import GraceDbConnector, KNCandidate
from .skymap import (distance_overlap, distance_overlap_valid, kn_skymap_3d, read_gw_skymap,
                     skymap_overlap_integral, write_skymap)
from .treatment.detectability import (DEFAULT_RANGES_MPC, HORIZON_OVER_RANGE,
                                      assess_detectability, scale_ranges)
from .treatment.distance_inclination import DistanceInclinationPosterior, posterior_for
from .treatment.far import DEFAULT_KN_CANDIDATE_RATE, coinc_far_kn
from .treatment.t0 import MergerTimeWindow, estimate_window

TAGS = ["ext_coinc"]


def query(gracedb: GraceDbConnector, window: MergerTimeWindow,
          far_max: Optional[float] = None) -> list[dict]:
    """Superevents whose t_0 lies in the merger-time window."""
    return gracedb.superevents(window.gps_start, window.gps_end, far_max=far_max)


def search(kn: KNCandidate, window: MergerTimeWindow, gracedb: GraceDbConnector,
           far_max: Optional[float] = None, upload: bool = True) -> list[dict]:
    """Look for GW superevents in the kilonova's window; log the result to
    GraceDB (on the kilonova's External event if it has one, and on each
    superevent found)."""
    found = query(gracedb, window, far_max)
    span = f"[{window.gps_start:.1f}, {window.gps_end:.1f}] (GPS, {window.duration / 3600:.2f} h, {window.method})"
    if upload:
        if not found and kn.graceid:
            gracedb.write_log(kn.graceid, f"RAVEN-KN: no superevent in t0 window {span}", tags=TAGS)
        for se in found:
            gracedb.write_log(se["superevent_id"],
                              f"RAVEN-KN: kilonova candidate {kn.name} (RA {kn.ra:.5f}, Dec {kn.dec:.5f}) "
                              f"has this superevent inside its t0 window {span}", tags=TAGS)
            if kn.graceid:
                gracedb.write_log(kn.graceid,
                                  f"RAVEN-KN: superevent {se['superevent_id']} in t0 window {span}",
                                  tags=TAGS)
    return found


def coinc_far(kn: KNCandidate, superevent: dict, window: MergerTimeWindow,
              gw_skymap: Optional[dict | str | Path] = None,
              post: Optional[DistanceInclinationPosterior] = None,
              kn_rate: float = DEFAULT_KN_CANDIDATE_RATE,
              ranges_mpc: dict[str, float] = DEFAULT_RANGES_MPC,
              chirp_mass_msun: Optional[float] = None) -> dict:
    """Joint FAR of a superevent and a kilonova candidate (RAVEN keys + distance).

    The distance factor is referenced to the network horizon (2.264 x the
    largest range), for the superevent's chirp mass if given (else BNS).
    If the kilonova lies beyond that horizon the factor is not meaningful
    (and the association is implausible on that ground alone).
    """
    far_gw = superevent.get("far")
    if not far_gw:
        raise ValueError("RAVEN-KN: superevent has no FAR")
    sky = dist = None
    dist_valid = None
    if gw_skymap is not None:
        gw = read_gw_skymap(gw_skymap) if not isinstance(gw_skymap, dict) else gw_skymap
        sky = skymap_overlap_integral(gw["PROB"], kn.ra, kn.dec, gw.get("nested", True))
        if post is not None and "DISTMU" in gw:
            d_ref = HORIZON_OVER_RANGE * max(scale_ranges(ranges_mpc, chirp_mass_msun).values())
            dist = distance_overlap(gw, kn.ra, kn.dec, post, d_ref_mpc=d_ref)
            dist_valid = distance_overlap_valid(post, d_ref)
    out = coinc_far_kn(window.duration, far_gw, kn_rate, sky, dist)
    out["distance_overlap_valid"] = dist_valid
    if dist_valid is False:
        out["distance_overlap_note"] = ("kilonova beyond the reference horizon: distance factor "
                                        "approximate; association disfavoured by detectability")
    out.update(superevent_id=superevent.get("superevent_id"), external_event=kn.graceid or kn.name)
    return out


def calc_signif_gracedb(kn: KNCandidate, superevent: dict, window: MergerTimeWindow,
                        gracedb: GraceDbConnector, **kw) -> dict:
    """Compute :func:`coinc_far` and upload ``coincidence_far_kn.json``."""
    res = coinc_far(kn, superevent, window, **kw)
    payload = json.dumps(res)
    msg = f"RAVEN-KN: joint FAR(s) in Hz with kilonova candidate {kn.name}"
    gracedb.write_log(superevent["superevent_id"], msg, "coincidence_far_kn.json", payload, TAGS)
    if kn.graceid:
        gracedb.write_log(kn.graceid, f"RAVEN-KN: joint FAR(s) in Hz with {superevent['superevent_id']}",
                          "coincidence_far_kn.json", payload, TAGS)
    return res


def assess(kn: KNCandidate, out_dir: Optional[str | Path] = None,
           ranges_mpc: dict[str, float] = DEFAULT_RANGES_MPC,
           snr_threshold: float = 6.0, p_min: float = 0.1,
           d_max_useful_mpc: Optional[float] = 500.0,
           segments=None, duty_cycle=None, d_max_prior: float = 1500.0,
           nside: int = 512, model_sys_frac: float = 0.75,
           chirp_mass_msun: Optional[float] = None) -> dict:
    """Targeted-search trigger decision for one kilonova candidate.

    Steps: (1) distance/inclination posterior from EM data; (2) t0 window,
    with the distance used to decide which non-detections constrain it;
    (3) expected GW network SNR over the window; (4) decision.
    If ``out_dir`` is given, write the decision, the 3D sky map, the
    posterior samples and the targeted-search request there.
    """
    post = posterior_for(kn, d_max=d_max_prior)
    d95 = float(np.quantile(post.d_mpc, 0.95)) if post.constraints else None
    window = estimate_window(kn, d_mpc_upper=d95, model_sys_frac=model_sys_frac)
    det = assess_detectability(post, window, kn.ra, kn.dec, ranges_mpc, snr_threshold,
                               p_min, d_max_useful_mpc, segments, duty_cycle,
                               chirp_mass_msun=chirp_mass_msun)
    decision = {
        "candidate": kn.name, "ra": kn.ra, "dec": kn.dec,
        "t0_window": {**asdict(window), "duration_h": window.duration / 3600},
        "distance_inclination": post.summary(),
        "search_priors": post.search_priors(),
        "detectability": det.to_dict(),
        "trigger_targeted_search": det.worth_search,
    }
    if out_dir is not None:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        stem = kn.name.replace(" ", "_")
        cols = kn_skymap_3d(kn.ra, kn.dec, post if post.constraints else None,
                            kn.pos_err_arcsec, nside)
        fits = write_skymap(out / f"{stem}_skymap.fits", cols, kn.name, kn.gpstime,
                            post if post.constraints else None)
        np.savez(out / f"{stem}_dist_incl_samples.npz", luminosity_distance=post.d_mpc,
                 cos_iota=post.cos_iota)
        request = {
            "candidate": kn.name, "ra_deg": kn.ra, "dec_deg": kn.dec,
            "gps_start": window.gps_start, "gps_end": window.gps_end,
            "priors": decision["search_priors"], "skymap": fits.name,
            "mchirp_range": ([0.9 * chirp_mass_msun, 1.1 * chirp_mass_msun] if chirp_mass_msun
                             else [1.0, 1.5]),
            "samples": f"{stem}_dist_incl_samples.npz",
            "note": "split windows longer than a few hours into sub-segments (Earth rotation)",
        }
        decision["files"] = {"skymap": fits.name, "samples": request["samples"]}
        if det.worth_search:
            (out / f"{stem}_targeted_search_request.json").write_text(json.dumps(request, indent=2))
            decision["files"]["targeted_search_request"] = f"{stem}_targeted_search_request.json"
        (out / f"{stem}_decision.json").write_text(json.dumps(decision, indent=2))
    return decision
