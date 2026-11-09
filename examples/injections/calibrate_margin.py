"""Calibrate the model-systematic margin of the light-curve t0 window with injections.

Synthetic kilonovae are drawn from one NMMA model (Bu2019lm or Ka2017),
observed with a ZTF-like survey, and fitted by RAVEN-KN's NMMA wrapper with
both models (self- and cross-model). For each fit, the statistical 90 %
interval [lo, hi] and the anchors (first detection, last constraining
non-detection) give the smallest margin fraction that covers the true
merger time:

    f_needed = (truth - hi) / (first detection - t0_median)          if truth > hi
             = (lo - truth) / (t0_median - early anchor, or late gap) if truth < lo
             = 0                                                      otherwise

as in ``raven_kn.treatment.t0.window_from_samples``.

Survey (ZTF-like cadence and depth, PS1 g/r/i model bands): epochs every night at a random phase (merger at t = 0), g/r/i 0.03 d
apart, 5-sigma depth 21.5, nights lost with probability 0.3, epochs from
-3 d to +6 d; in half of the injections the survey has no pre-merger
coverage (late discovery, as AT2017gfo). Distances uniform in volume in
40-150 Mpc; host redshift known (informed distance prior, as in practice).

  python calibrate_margin.py make   --out runs_inj [--n 12 --seed 1 --survey ztf|rubin]
  python calibrate_margin.py run    --out runs_inj [--workers 6]
  python calibrate_margin.py report --out runs_inj
"""
import argparse
import json
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from raven_kn.connection import (KNCandidate, PhotometryPoint, Redshift, gps_to_mjd,
                                 load_candidate, save_candidate)
from raven_kn.treatment import lightcurve_fit as lcf
from raven_kn.treatment.t0 import estimate_window, window_from_bracket

NMMA_PY = os.path.expanduser("~/miniforge3/envs/nmma/bin/python")
EXE = os.path.expanduser(os.environ.get("RAVEN_KN_NMMA_EXE",
                                        "~/miniforge3/envs/nmma/bin/lightcurve-analysis"))
SVD = os.path.expanduser("~/nmma_models")
FILTERS = {"ps1::g": "g", "ps1::r": "r", "ps1::i": "i"}
T_BASE = 1400000000.0
HERE = Path(__file__).parent


def draw_params(model, rng):
    if model == "Bu2019lm":
        return {"log10_mej_dyn": rng.uniform(-2.8, -1.8), "log10_mej_wind": rng.uniform(-1.9, -1.0),
                "KNphi": rng.uniform(15, 75),
                "inclination_EM": float(np.arccos(rng.uniform(np.cos(np.radians(70)), 1)))}
    return {"log10_mej": rng.uniform(-2.5, -1.3), "log10_vej": rng.uniform(-1.3, -0.7),
            "log10_Xlan": rng.uniform(-6, -2)}


SURVEY = {"depth": 21.5, "cadence_d": 1.0, "lost": 0.3, "bands_per_epoch": 3,
          "d_min": 40.0, "d_max": 150.0, "fit_days": 3.0}
PRESETS = {
    "ztf": {},
    # Rubin WFD-like: deeper, a field revisited every ~3 nights, 2 filters per night
    "rubin": {"depth": 24.5, "cadence_d": 3.0, "lost": 0.2, "bands_per_epoch": 2, "d_max": 300.0,
              "fit_days": 6.5},
}


def survey_times(rng, pre_merger):
    cad = SURVEY["cadence_d"]
    phase = rng.uniform(0, cad)
    first = -3 if pre_merger else int(rng.integers(0, 2))
    epochs = [k * cad + phase for k in range(first, int(6 / cad) + 2) if k * cad + phase <= 6.0]
    epochs = [t for t in epochs if rng.uniform() > SURVEY["lost"]]
    times, bands = [], []
    for t in epochs:
        fs = list(FILTERS)
        if SURVEY["bands_per_epoch"] < len(fs):
            fs = sorted(rng.choice(fs, SURVEY["bands_per_epoch"], replace=False),
                        key=list(FILTERS).index)
        for j, f in enumerate(fs):
            times.append(t + 0.03 * j)
            bands.append(f)
    return np.array(times), bands


def make(out: Path, n: int, seed: int):
    rng = np.random.default_rng(seed)
    out.mkdir(parents=True, exist_ok=True)
    meta, jobs = [], []
    while len(meta) < 2 * n:
        model = "Bu2019lm" if len(meta) % 2 == 0 else "Ka2017"
        lo, hi = SURVEY["d_min"], SURVEY["d_max"]
        d = (lo ** 3 + rng.uniform() * (hi ** 3 - lo ** 3)) ** (1 / 3)
        par = dict(draw_params(model, rng), luminosity_distance=d)
        pre = bool(rng.uniform() < 0.5)
        t, b = survey_times(rng, pre)
        meta.append({"truth_model": model, "params": par, "pre_merger_coverage": pre,
                     "t": t.tolist(), "bands": b})
        jobs.append({"model": model, "params": par, "filters": list(FILTERS), "t_days": t.tolist()})
    fin, fout = out / "truth_in.json", out / "truth_out.json"
    fin.write_text(json.dumps(jobs))
    subprocess.run([NMMA_PY, str(HERE / "nmma_truth.py"), SVD, str(fin), str(fout)], check=True,
                   stdout=subprocess.DEVNULL)
    mags = json.loads(fout.read_text())
    kept = []
    for i, (m, lc) in enumerate(zip(meta, mags)):
        t0_gps = T_BASE + i * 1e6
        phot, ndet = [], 0
        for k, (t, f) in enumerate(zip(m["t"], m["bands"])):
            mag = lc[f][k] if t > 0 else np.inf
            mjd = gps_to_mjd(t0_gps + t * 86400)
            if np.isfinite(mag) and mag < SURVEY["depth"]:
                sig = float(np.hypot(0.03, 0.2 * 10 ** (0.4 * (mag - SURVEY["depth"]))))
                phot.append(PhotometryPoint(mjd, FILTERS[f], float(mag + rng.normal(0, sig)), sig,
                                            survey="PS1"))
                ndet += 1
            else:
                phot.append(PhotometryPoint(mjd, FILTERS[f], SURVEY["depth"], upper_limit=True,
                                            survey="PS1"))
        first = min((p.mjd for p in phot if not p.upper_limit), default=None)
        n_early = sum(1 for p in phot if not p.upper_limit and first is not None and p.mjd < first + SURVEY["fit_days"])
        if n_early < 3:
            continue
        from astropy.cosmology import Planck18, z_at_value
        import astropy.units as u
        z = float(z_at_value(Planck18.luminosity_distance, m["params"]["luminosity_distance"] * u.Mpc))
        name = f"inj{i:03d}_{m['truth_model']}"
        kn = KNCandidate(name, 150.0, 10.0, photometry=phot, redshift=Redshift(z, 1e-4, source="host"))
        save_candidate(kn, out / f"{name}.json")
        kept.append({"name": name, "t0_gps": t0_gps, **{k: m[k] for k in
                     ("truth_model", "params", "pre_merger_coverage")}})
    (out / "injections.json").write_text(json.dumps(kept, indent=1))
    (out / "survey.json").write_text(json.dumps(SURVEY, indent=1))
    print(f"{len(kept)} injections kept of {len(meta)}")


def _fit(args):
    out, name, model, fit_days = args
    kn = load_candidate(out / f"{name}.json")
    wd = out / "fits" / f"{name}__{model}"
    if (wd / "nmma_out").exists() and any((wd / "nmma_out").glob("*posterior*")):
        return name, model, "exists"
    try:
        lcf.fit(kn, wd, executable=EXE, interpolation_type="sklearn_gp", model=model,
                svd_path=SVD, fit_days=fit_days, extra_args=["--nlive", "512", "--local-only"])
        return name, model, "ok"
    except Exception as e:      # keep going; reported
        return name, model, f"failed: {e}"


def run(out: Path, workers: int):
    inj = json.loads((out / "injections.json").read_text())
    sv = json.loads((out / "survey.json").read_text()) if (out / "survey.json").exists() else SURVEY
    tasks = [(out, i["name"], m, sv.get("fit_days", 3.0)) for i in inj for m in ("Ka2017", "Bu2019lm")]
    with ProcessPoolExecutor(workers) as ex:
        for name, model, status in ex.map(_fit, tasks):
            print(name, model, status, flush=True)


def report(out: Path):
    inj = json.loads((out / "injections.json").read_text())
    rows = []
    for i in inj:
        for model in ("Bu2019lm", "Ka2017"):
            kn = load_candidate(out / f"{i['name']}.json")
            try:
                lcf.import_result(kn, out / "fits" / f"{i['name']}__{model}")
            except Exception:
                continue
            w0 = window_from_bracket(kn)        # anchors (bracket alone: max age 7 d)
            s = np.asarray(kn.t0_samples_gps)
            lo, med, hi = np.quantile(s, [0.05, 0.5, 0.95])
            truth = i["t0_gps"]
            late_gap = max(w0.anchor_late_gps - med, 0)
            early_gap = max(med - w0.anchor_early_gps, 0) if w0.anchor_early_gps else late_gap
            if truth > hi:
                need = (truth - hi) / max(late_gap, 1)
            elif truth < lo:
                need = (lo - truth) / max(early_gap, 1)
            else:
                need = 0.0
            # final window as used by RAVEN-KN (fit + margin, intersected with the bracket)
            win = {f: estimate_window(kn, model_sys_frac=f) for f in (0.0, 0.5, 0.75, 1.0)}
            rows.append({"injection": i["name"], "truth_model": i["truth_model"], "fit_model": model,
                         "pre_merger_coverage": i["pre_merger_coverage"],
                         "d_mpc": i["params"]["luminosity_distance"],
                         "first_det_h": (w0.anchor_late_gps - truth) / 3600,
                         "t0_med_minus_true_h": (med - truth) / 3600,
                         "stat_width_h": (hi - lo) / 3600,
                         "late_gap_h": late_gap / 3600, "f_needed": need,
                         "bracket_h": w0.duration / 3600,
                         **{f"width_h_f{f}": w.duration / 3600 for f, w in win.items()},
                         **{f"covered_f{f}": bool(w.gps_start <= truth <= w.gps_end)
                            for f, w in win.items()}})
    (out / "report.json").write_text(json.dumps(rows, indent=1))
    print(f"{'injection':22s} {'fit':9s} pre  1stdet_h  bias_h  stat_h  gap_h  f_need")
    for r in rows:
        print(f"{r['injection']:22s} {r['fit_model']:9s} {'y' if r['pre_merger_coverage'] else 'n'}  "
              f"{r['first_det_h']:8.1f} {r['t0_med_minus_true_h']:7.1f} {r['stat_width_h']:7.1f} "
              f"{r['late_gap_h']:6.1f} {r['f_needed']:6.2f}")
    for label, sel in (("self-model", lambda r: r["truth_model"] == r["fit_model"]),
                       ("cross-model", lambda r: r["truth_model"] != r["fit_model"]),
                       ("all", lambda r: True)):
        f = np.array([r["f_needed"] for r in rows if sel(r)])
        if len(f) == 0:
            continue
        cov = {x: float(np.mean(f <= x)) for x in (0, 0.25, 0.5, 0.75, 1.0, 1.5)}
        print(f"{label:12s} n={len(f):3d}  coverage vs frac: " +
              "  ".join(f"{k}:{v:.2f}" for k, v in cov.items()) +
              f"   90% quantile of f_needed: {np.quantile(f, 0.9):.2f}")
        sub = [r for r in rows if sel(r)]
        for x in (0.0, 0.5, 0.75, 1.0):
            print(f"{'':12s} final window (with bracket), frac {x}: coverage "
                  f"{np.mean([r[f'covered_f{x}'] for r in sub]):.2f}, median width "
                  f"{np.median([r[f'width_h_f{x}'] for r in sub]):.1f} h "
                  f"(bracket alone {np.median([r['bracket_h'] for r in sub]):.1f} h)")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("action", choices=["make", "run", "report"])
    p.add_argument("--out", default="runs_inj")
    p.add_argument("--n", type=int, default=12)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--survey", choices=sorted(PRESETS), default="ztf",
                   help="survey preset used by 'make' (stored with the injections)")
    a = p.parse_args()
    out = Path(a.out)
    SURVEY.update(PRESETS[a.survey])
    {"make": lambda: make(out, a.n, a.seed), "run": lambda: run(out, a.workers),
     "report": lambda: report(out)}[a.action]()


if __name__ == "__main__":
    main()
