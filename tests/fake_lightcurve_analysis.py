#!/usr/bin/env python
"""Stand-in for NMMA's ``lightcurve-analysis``: same arguments, writes a
bilby-style ``<label>_posterior_samples.dat`` whose t0 is the GW170817
merger +- 0.1 d (FAKE_T0_MJD env var), within the prior bounds."""
import argparse
import os
from pathlib import Path

import numpy as np

p = argparse.ArgumentParser()
for a in ("--em-model", "--light-curve-data", "--trigger-time", "--prior-file", "--outdir",
          "--label", "--em-tmin", "--em-tmax", "--interpolation-type", "--sampler",
          "--em-error-budget", "--svd-path"):
    p.add_argument(a)
p.add_argument("--plot", action="store_true")
a = p.parse_args()

# the inputs must be readable the way NMMA reads them
from astropy.time import Time  # noqa: E402
for line in Path(a.light_curve_data).read_text().split("\n"):
    if line:
        t, f, m, e = line.split()
        Time(t), float(m), float(e)

trigger = float(a.trigger_time)
t0 = float(os.environ["FAKE_T0_MJD"])
rng = np.random.default_rng(3)
ts = rng.normal(t0 - trigger, 0.1, 4000)
d = rng.normal(43, 3, 4000)
inc = np.abs(rng.normal(0.35, 0.08, 4000))
out = Path(a.outdir)
out.mkdir(parents=True, exist_ok=True)
np.savetxt(out / f"{a.label}_posterior_samples.dat", np.c_[ts, d, inc],
           header="timeshift luminosity_distance inclination_EM", comments="")
