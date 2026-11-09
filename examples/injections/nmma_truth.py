"""Apparent magnitudes of NMMA SVD kilonova models (run with the nmma env's python).

  python nmma_truth.py SVD_PATH in.json out.json

in.json: list of {"model": "Bu2019lm"|"Ka2017", "params": {...model parameters,
luminosity_distance (Mpc), [inclination_EM (rad)]}, "filters": ["ps1::g", ...],
"t_days": [...]} with times in days after the merger (observer frame).
out.json: list of {filter: [mag, ...]} (inf outside the model's time range).
"""
import json
import os
import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")
svd_path, fin, fout = sys.argv[1:4]
from astropy.cosmology import Planck18, z_at_value  # noqa: E402
import astropy.units as u  # noqa: E402
from nmma.em.model import SVDLightCurveModel  # noqa: E402

jobs = json.load(open(fin))
models, out = {}, []
for job in jobs:
    key = (job["model"], tuple(job["filters"]))
    if key not in models:
        models[key] = SVDLightCurveModel(job["model"], svd_path=os.path.expanduser(svd_path),
                                         interpolation_type="sklearn_gp", filters=job["filters"])
    m = models[key]
    p = m.parameter_conversion(dict(job["params"], timeshift=0.0))
    p["redshift"] = float(z_at_value(Planck18.luminosity_distance,
                                     p["luminosity_distance"] * u.Mpc))
    m.redshift_func = lambda q: q["redshift"]
    m.em_parameter_setup(p, combine_params=False)
    obs_t, lc = m.gen_detector_lc(p)
    t = np.asarray(job["t_days"], float)
    out.append({f: np.interp(t, obs_t, np.asarray(lc[f]), left=np.inf, right=np.inf).tolist()
                for f in job["filters"]})
json.dump(out, open(fout, "w"))
