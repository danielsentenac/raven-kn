"""Check the NMMA wrapper on the GRB 211211A kilonova (z = 0.0763).

Truth: the merger time is taken as the Swift-BAT / Fermi-GBM trigger,
2021-12-11T13:09:59 UTC (GCN 31202). The GRB is HIDDEN from the fit, as if
the kilonova had been found by an optical survey alone.

Data: NMMA example file example_files/lightcurves/GRB211211A.txt (tag
v1.0.2); only optical/NIR points (PS1 g,r,i,z, Bessell r,i, 2MASS J,H,Ks)
are used; UVOT/X-ray/radio are dropped. The early optical points include
afterglow light, which a kilonova-only model does not describe.

Case C: host redshift as distance prior (spectroscopy available).
Case D: no distance information (kilonova-only distance).

Usage:
  python examples/validate_grb211211a.py run C|D --data GRB211211A.txt --workdir runs/
  python examples/validate_grb211211a.py report --data GRB211211A.txt --workdir runs/
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from astropy.time import Time

from raven_kn.connection import KNCandidate, PhotometryPoint, Redshift
from raven_kn.treatment import lightcurve_fit as lcf
from raven_kn.treatment.distance_inclination import posterior_for
from raven_kn.treatment.t0 import estimate_window, window_from_bracket

T0_GPS = float(Time("2021-12-11T13:09:59", format="isot", scale="utc").gps)
RA, DEC = 212.2920, 27.8886            # optical afterglow, approximate
Z_HOST = 0.0763
KEEP = {"ps1::g": ("PS1", "g"), "ps1::r": ("PS1", "r"), "ps1::i": ("PS1", "i"),
        "ps1::z": ("PS1", "z"), "bessellr": ("", "bessellr"), "besselli": ("", "besselli"),
        "2massj": ("", "2massj"), "2massh": ("", "2massh"), "2massks": ("", "2massks")}
EXE = os.path.expanduser(os.environ.get("RAVEN_KN_NMMA_EXE",
                                        "~/miniforge3/envs/nmma/bin/lightcurve-analysis"))


def candidate(data: str, case: str) -> KNCandidate:
    phot, seen = [], set()
    for line in Path(data).read_text().split("\n"):
        if not line.strip():
            continue
        t, f, m, e = line.split()
        if f not in KEEP or (t, f) in seen:          # file has duplicated rows
            continue
        seen.add((t, f))
        survey, band = KEEP[f]
        e = float(e)
        phot.append(PhotometryPoint(Time(t, format="isot", scale="utc").mjd, band, float(m),
                                    None if np.isinf(e) else e, upper_limit=bool(np.isinf(e)),
                                    survey=survey))
    kn = KNCandidate(f"GRB211211A_{case}", RA, DEC, photometry=phot)
    if case == "C":
        kn.redshift = Redshift(Z_HOST, 1e-4, source="host galaxy")
    return kn


def main():
    p = argparse.ArgumentParser()
    p.add_argument("action", choices=["run", "report"])
    p.add_argument("case", nargs="?", choices=["C", "D"])
    p.add_argument("--data", required=True)
    p.add_argument("--workdir", default="runs")
    p.add_argument("--nlive", default="512")
    p.add_argument("--model", default="Bu2019lm")
    p.add_argument("--fit-days", type=float, default=3.0, help="<= 0: all data")
    a = p.parse_args()
    if a.action == "run":
        kn = candidate(a.data, a.case)
        lcf.fit(kn, Path(a.workdir) / a.case, executable=EXE, interpolation_type="sklearn_gp", model=a.model,
                        fit_days=a.fit_days if a.fit_days > 0 else None,
                svd_path=os.path.expanduser("~/nmma_models"),
                extra_args=["--nlive", a.nlive, "--local-only"])
        print("done", a.case)
        return
    d_true = float(np.median(posterior_for(candidate(a.data, "C")).d_mpc))
    rows = {}
    for case in ("C", "D"):
        kn = candidate(a.data, case)
        br = window_from_bracket(kn, d_mpc_upper=500)
        try:
            lcf.import_result(kn, Path(a.workdir) / case)
        except FileNotFoundError as e:
            rows[case] = {"error": str(e)}
            continue
        t0 = np.asarray(kn.t0_samples_gps)
        lo, med, hi = np.quantile(t0, [0.05, 0.5, 0.95])
        r = {"first_detection_after_true_t0_h": (kn.first_detection.gps - T0_GPS) / 3600,
             "t0_median_minus_true_h": (med - T0_GPS) / 3600,
             "t0_90pc_width_h": (hi - lo) / 3600,
             "true_t0_in_90pc": bool(lo <= T0_GPS <= hi),
             "bracket_width_h": br.duration / 3600}
        for f in (0.0, 0.5, 0.75):
            w = estimate_window(kn, d_mpc_upper=500, model_sys_frac=f)
            r[f"window_frac{f}_h"] = w.duration / 3600
            r[f"true_t0_in_window_frac{f}"] = bool(w.gps_start <= T0_GPS <= w.gps_end)
        dl = np.asarray(kn.lc_fit_samples["luminosity_distance"])
        r["D_L_median_mpc"] = float(np.median(dl))
        r["D_L_90pc_mpc"] = [float(x) for x in np.quantile(dl, [0.05, 0.95])]
        r["D_L_from_redshift_mpc"] = d_true
        r["theta_obs_median_deg"] = float(np.degrees(np.median(kn.lc_fit_samples["inclination"])))
        rows[case] = r
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
