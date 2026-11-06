"""RAVEN-KN applied to S250818k / AT2025ulz (ZTF25abjmnps).

Public inputs:
- S250818k (GraceDB): t0 = GPS 1439515224.030 (2025-08-18 01:20:06 UTC),
  FAR = 6.8e-8 Hz (~2.1/yr, low significance), H1 L1 V1; GstLAL-SSM and
  PyCBC Live; HasSSM > 99 %, chirp mass bin 0.1-0.87 Msun (GCN 41437).
  Sky maps: bayestar.fits.gz and Bilby.offline0.fits.gz (GCN 41440).
- AT2025ulz (GCN 41414): RA 237.9757129, Dec +30.9023146; first ZTF epoch
  2025-08-18 04:02 UTC (r = 21.29), pre-discovery limits all shallower than
  the transient (no constraining non-detection). Host redshift 0.084840
  (DESI, arXiv:2510.23723; Keck 0.0848, GCN 41436).

Usage: python examples/s250818k_case.py --maps DIR   (DIR holds the two FITS files;
       download from https://gracedb.ligo.org/api/superevents/S250818k/files/<name>)
"""
import argparse
import json
from pathlib import Path

import numpy as np
from astropy.time import Time

from raven_kn.connection import KNCandidate, PhotometryPoint, Redshift
from raven_kn.skymap import distance_overlap, read_gw_skymap, skymap_overlap_integral
from raven_kn.treatment.detectability import (BNS_CHIRP_MASS, HORIZON_OVER_RANGE,
                                              assess_detectability, scale_ranges)
from raven_kn.treatment.distance_inclination import posterior_for
from raven_kn.treatment.far import coinc_far_kn
from raven_kn.treatment.t0 import window_from_bracket

GW_T0 = 1439515224.030
GW_FAR = 6.81022356082428e-08
RA, DEC = 237.9757129, 30.9023146


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--maps", required=True)
    p.add_argument("--ranges", default="H1=150,L1=160,V1=55",
                   help="approximate O4c BNS ranges (Mpc); replace with the day's values")
    a = p.parse_args()
    ranges_bns = {k: float(v) for k, v in (kv.split("=") for kv in a.ranges.split(","))}

    first = Time("2025-08-18T04:02:00", scale="utc").mjd
    kn = KNCandidate("AT2025ulz", RA, DEC, pos_err_arcsec=0.5, source="ZTF/GCN 41414",
                     photometry=[PhotometryPoint(first, "r", 21.29, 0.13, survey="ZTF")],
                     redshift=Redshift(0.084840, 1e-5, source="DESI host"))
    post = posterior_for(kn)
    win = window_from_bracket(kn)           # no constraining limit -> 7 d before first detection
    out = {"distance_from_redshift_mpc": post.summary()["luminosity_distance_mpc"],
           "t0_window": {"start_minus_gw_h": (win.gps_start - GW_T0) / 3600,
                         "end_minus_gw_h": (win.gps_end - GW_T0) / 3600,
                         "duration_h": win.duration / 3600, "notes": win.notes,
                         "contains_gw_time": bool(win.gps_start <= GW_T0 <= win.gps_end)}}

    for name in ("bayestar.fits.gz", "Bilby.offline0.fits.gz"):
        gw = read_gw_skymap(Path(a.maps) / name)
        import healpy as hp
        nside = hp.npix2nside(len(gw["PROB"]))
        ipix = hp.ang2pix(nside, RA, DEC, lonlat=True, nest=True)
        order = np.argsort(gw["PROB"])[::-1]
        cred = np.empty_like(gw["PROB"])
        cred[order] = np.cumsum(gw["PROB"][order])
        area90 = float((np.cumsum(gw["PROB"][order]) < 0.9).sum() * hp.nside2pixarea(nside, degrees=True))
        sky = skymap_overlap_integral(gw["PROB"], RA, DEC)
        res = {"area90_deg2": area90, "kn_credible_level": float(cred[ipix]),
               "gw_distance_at_kn_pixel_mpc": [float(gw["DISTMU"][ipix]), float(gw["DISTSIGMA"][ipix])],
               "sky_overlap": sky}
        for label, mc in (("BNS horizon", BNS_CHIRP_MASS), ("sub-solar Mc=0.5", 0.5)):
            dref = HORIZON_OVER_RANGE * max(scale_ranges(ranges_bns, mc).values())
            dist = distance_overlap(gw, RA, DEC, post, d_ref_mpc=dref)
            far = coinc_far_kn(win.duration, GW_FAR, sky_overlap=sky, distance_overlap=dist)
            res[f"joint_far ({label}, d_ref={dref:.0f} Mpc)"] = {
                "distance_overlap": dist,
                "temporal_far_per_yr": far["temporal_coinc_far"] * 3.156e7,
                "spatiotemporal_far_per_yr": far["spatiotemporal_coinc_far"] * 3.156e7,
                "with_distance_far_per_yr": (far["spatiotemporal_distance_coinc_far"] or np.nan) * 3.156e7}
        out[name] = res

    for label, mc in (("BNS 1.4+1.4", None), ("sub-solar Mc=0.5", 0.5), ("sub-solar Mc=0.87", 0.87)):
        det = assess_detectability(post, win, RA, DEC, ranges_bns, d_max_useful_mpc=None,
                                   chirp_mass_msun=mc)
        out[f"detectability ({label})"] = {"ranges_mpc": {k: round(v) for k, v in
                                                          scale_ranges(ranges_bns, mc).items()},
                                          "snr_median": det.snr_median, "snr_90": det.snr_90,
                                          "p_snr_ge_6": det.p_detect}
    print(json.dumps(out, indent=2, default=float))


if __name__ == "__main__":
    main()
