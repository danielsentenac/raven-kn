"""End-to-end check of the NMMA wrapper on AT2017gfo / GW170817.

Case A: real AT2017gfo PS1 g,r,i,z photometry (NMMA example file
        example_files/lightcurves/AT2017gfo_reduced.dat), first point +0.7 d.
Case B: same + synthetic g, r non-detections 3 h before the merger
        (Jamsin et al. arXiv:2610.28775 "what if" scenario).
Both use the NGC 4993 redshift and the VLBI viewing angle as priors.

Usage:
  python examples/validate_at2017gfo.py run A|B --data AT2017gfo_reduced.dat --workdir runs/
  python examples/validate_at2017gfo.py report --workdir runs/ --data AT2017gfo_reduced.dat
NMMA executable: $RAVEN_KN_NMMA_EXE (default ~/miniforge3/envs/nmma/bin/lightcurve-analysis).
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from astropy.time import Time

from raven_kn.connection import (KNCandidate, PhotometryPoint, Redshift,
                                 ViewingAngleMeasurement, gps_to_mjd)
from raven_kn.treatment import lightcurve_fit as lcf
from raven_kn.treatment.t0 import estimate_window, window_from_bracket

T0_GPS = 1187008882.43            # GW170817 merger
BANDS = {"ps1::g": "g", "ps1::r": "r", "ps1::i": "i", "ps1::z": "z"}
EXE = os.path.expanduser(os.environ.get("RAVEN_KN_NMMA_EXE",
                                        "~/miniforge3/envs/nmma/bin/lightcurve-analysis"))


def candidate(data: str, case: str) -> KNCandidate:
    phot = []
    for line in Path(data).read_text().split("\n"):
        if not line.strip():
            continue
        t, f, m, e = line.split()
        if f not in BANDS:
            continue
        e = float(e)
        phot.append(PhotometryPoint(Time(t, format="isot", scale="utc").mjd, BANDS[f], float(m),
                                    None if np.isinf(e) else e, upper_limit=bool(np.isinf(e)),
                                    survey="PS1"))
    if case == "B":
        t0 = gps_to_mjd(T0_GPS)
        phot += [PhotometryPoint(t0 - 3 / 24, "g", 24.5, upper_limit=True, survey="PS1"),
                 PhotometryPoint(t0 - 3 / 24 + 0.02, "r", 24.3, upper_limit=True, survey="PS1")]
    return KNCandidate(f"AT2017gfo_{case}", 197.450374, -23.381495, photometry=phot,
                       redshift=Redshift(0.009783, 1e-5, source="NGC 4993"),
                       viewing_angle=ViewingAngleMeasurement(20.0, 5.0, "VLBI"))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("action", choices=["run", "report"])
    p.add_argument("case", nargs="?", choices=["A", "B"])
    p.add_argument("--data", required=True)
    p.add_argument("--workdir", default="runs")
    p.add_argument("--nlive", default="512")
    p.add_argument("--model", default="Bu2019lm")
    p.add_argument("--fit-days", type=float, default=3.0, help="<= 0: all data")
    a = p.parse_args()
    if a.action == "run":
        kn = candidate(a.data, a.case)
        kn, s = lcf.fit(kn, Path(a.workdir) / a.case, executable=EXE,
                        interpolation_type="sklearn_gp", model=a.model,
                        fit_days=a.fit_days if a.fit_days > 0 else None, svd_path=os.path.expanduser("~/nmma_models"),
                        extra_args=["--nlive", a.nlive, "--local-only"])
        print("done", a.case)
        return
    rows = {}
    for case in ("A", "B"):
        kn = candidate(a.data, case)
        br = window_from_bracket(kn, d_mpc_upper=55)
        try:
            lcf.import_result(kn, Path(a.workdir) / case)
        except FileNotFoundError as e:
            rows[case] = {"error": str(e)}
            continue
        t0 = np.asarray(kn.t0_samples_gps)
        lo, med, hi = np.quantile(t0, [0.05, 0.5, 0.95])
        w = estimate_window(kn, d_mpc_upper=55)
        rows[case] = {
            "t0_median_minus_true_h": (med - T0_GPS) / 3600,
            "t0_90pc_width_h": (hi - lo) / 3600,
            "true_t0_in_90pc": bool(lo <= T0_GPS <= hi),
            "bracket_width_h": br.duration / 3600,
            "final_window_h": w.duration / 3600,
            "true_t0_in_final_window": bool(w.gps_start <= T0_GPS <= w.gps_end),
            "window_notes": w.notes,
            "D_L_median_mpc": float(np.median(kn.lc_fit_samples["luminosity_distance"])),
            "theta_obs_median_deg": float(np.degrees(np.median(kn.lc_fit_samples["inclination"]))),
        }
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
