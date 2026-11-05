# RAVEN-KN

Kilonova-triggered coincidence and targeted-GW-search trigger pipeline, built
with the layout of LVK RAVEN (`ligo.raven` 3.3) so that it can later be merged
into it or called from gwcelery.

Direction: **kilonova candidate → is a GW search worth it? → merger-time window,
sky position, distance and inclination priors → coincidence with existing GW
candidates (joint FAR) and/or a targeted search of archived strain.**

## Modules

| RAVEN-KN | Role | RAVEN analogue |
|---|---|---|
| `raven_kn.connection` | `KNCandidate` (photometry incl. non-detections, redshift, distance, viewing angle, GRB, light-curve-fit samples); connectors: Fink (ZTF API), GraceDB, GCN JSON, local JSON | `gracedb_events.ExtTrig` / `SE` |
| `raven_kn.treatment.t0` | merger-time window [tl, th]: GRB (seconds) > light-curve-fit t0 posterior ∩ photometric bracket (last *constraining* non-detection → first detection) | fixed `tl, th` |
| `raven_kn.treatment.lightcurve_fit` | NMMA wrapper: candidate → NMMA data/prior files and `lightcurve-analysis` command → posterior read back as `t0_samples_gps` (+ D_L, ι) | — (new) |
| `raven_kn.treatment.distance_inclination` | joint (D_L, cos ι) posterior from spectroscopic redshift, host distance, VLBI/afterglow viewing angle, prompt GRB, light-curve-fit samples; search priors | — (new) |
| `raven_kn.treatment.detectability` | expected network SNR over the window (LAL antenna patterns at the candidate position, Earth rotation, duty cycle/segments); decision: trigger or not | — (new) |
| `raven_kn.treatment.far` | joint FARs: `(th−tl) R_KN FAR_GW` / sky overlap / **distance overlap**; targeted-search p-value, sub-segment trials, RAVEN "SubGRBTargeted" product formula | `search.coinc_far` |
| `raven_kn.skymap` | kilonova 2D/3D HEALPix map (PROB, DISTMU, DISTSIGMA, DISTNORM), GW map reader, sky overlap (RAVEN's own function if installed), distance overlap | `search.skymap_overlap_integral` |
| `raven_kn.search` | `query`, `search`, `coinc_far`, `calc_signif_gracedb` (GraceDB logs `RAVEN-KN:`, tag `ext_coinc`, file `coincidence_far_kn.json`), `assess` | `search.*` |

`KNCandidate.to_raven_event_dict()` gives a GraceDB-External-event-shaped dict
accepted by `ligo.raven.gracedb_events.ExtTrig(..., event_dict=..., use_radec=True)`.

## Use

```bash
python -m venv --system-site-packages .venv && .venv/bin/pip install -e .[test]
.venv/bin/python examples/make_examples.py
.venv/bin/raven-kn assess examples/gw170817_like.json --out examples/out --ranges H1=150,L1=160,V1=55
.venv/bin/raven-kn coinc  examples/gw170817_like.json --superevent-far 1e-7 --gw-skymap bayestar.fits.gz
.venv/bin/python -m pytest -q -p no:anyio tests
```

`assess` writes `<name>_decision.json`, `<name>_skymap.fits` (3D),
`<name>_dist_incl_samples.npz` and, if the search is worth it,
`<name>_targeted_search_request.json` (position, GPS window, D_L and
inclination priors).

Example results (examples/, photometry synthetic, O4-like ranges):

| candidate | window | D_L (90 %) | θ_obs | network SNR | trigger |
|---|---|---|---|---|---|
| AT2017gfo-like, UL 3 h before merger | 19.4 h (bracket) | 44.6 [37, 52] Mpc | 21° [13, 29] | ~51 | yes |
| same + GRB | 6 s | same | same | ~68 | yes |
| z = 0.15 host | 7 d (no constraining UL) | 737 Mpc | unconstrained | ~2 | no (> 500 Mpc, P(SNR≥6)=0) |

## NMMA light-curve fit

```bash
raven-kn nmma examples/gw170817_like.json --workdir runs/gfo --dry-run   # write inputs, print command
raven-kn nmma examples/gw170817_like.json --workdir runs/gfo --save cand_fit.json          # run locally
raven-kn nmma examples/gw170817_like.json --workdir runs/gfo --import-only --save cand_fit.json  # after a cluster run
raven-kn assess cand_fit.json --out out/      # window now = NMMA t0 interval ∩ photometric bracket
```

Interface checked against the NMMA 1.0.1 source (not yet run with NMMA itself):
- data lines `isot filter mag err`, non-detections with `err = inf`; survey bands mapped to NMMA filters (`lsstg`, `ztfr`, `ps1::i`, ...);
- NMMA drops data before `--trigger-time`, so the trigger is placed before the earliest non-detection kept; merger = trigger + `timeshift` (days);
- `timeshift` prior: from first detection − `max_age` (7 d) to first detection − (1+z_max)·`em_tmin`, as required by NMMA's model/data time check;
- priors: model parameters from NMMA's prior files (Bu2019lm, Bu2023Ye built in; any other model, e.g. Bu2026, via `--model-prior`), distance as a truncated Gaussian from the EM distance posterior, `inclination_EM` from the viewing angle (else `Sine`);
- with these informed priors the fitted D_L/ι are **not** reused as an extra constraint (`lc_fit_informed`), to avoid counting the redshift twice; use `--uninformed` to get an independent kilonova distance.

## Validation with real NMMA (v1.0.2, Bu2019lm, sklearn_gp, pymultinest 512 live points)

`examples/validate_at2017gfo.py` on NMMA's AT2017gfo PS1 g,r,i,z light curve (first point +0.7 d),
NGC 4993 redshift and VLBI angle as priors; truth = GW170817 merger time. About 4 min per run.

| case | t0 median − true | 90 % stat. window | truth inside | window with model margin (frac 0.5) | truth inside |
|---|---|---|---|---|---|
| A: real data only | −41.5 h | 38.4 h | **no** | 96.8 h (±29.2 h) | yes |
| B: + synthetic g, r non-detections 3 h before merger | −1.9 h | 8.6 h | yes | 19.0 h (−1.0 / +9.4 h) | yes |

A is biased early, as Jamsin et al. found for Bu2019lm on AT2017gfo (too high early opacity → too long a rise);
the pre-merger non-detection removes the bias and narrows the window 4.5x. B's non-detection is placed
relative to the known merger time (best case by construction). Results in `examples/validation/report.json`.

Environment notes: NMMA runs from the conda env `nmma` (`--executable` or `$RAVEN_KN_NMMA_EXE`);
models in `~/nmma_models` (`python -c "from nmma.core.gitlab import main; main()" --model Bu2019lm --filters ...`,
the `svdmodel-download` entry point is broken in v1.0.2); some sncosmo bandpasses redirect to SVO's dead
`http://` URL and had to be fetched once from `https://svo2.cab.inta-csic.es/theory/fps/` into sncosmo's cache.

## Assumptions to revisit (all are parameters)

- **Model-systematic margin** on a light-curve t0 (`window_from_samples`, `--model-sys-frac`, default 0.5):
  the 90 % interval is widened by 0.5 x the extrapolation gap from the median t0 to the data that anchor it
  (first detection on the late side; last constraining non-detection, else the same gap, on the early side),
  at least 1 h, then intersected with the photometric bracket. 0.5 is the smallest round value that recovers
  GW170817 in case A: calibrated on one event and one model; to be redone with cross-model injections.
- **Bracket model** (`window_from_bracket`): kilonova at least as bright as
  M = −15 from 0.25 d after merger; a non-detection constrains only if deeper
  than that at the 95th-percentile distance; 7 d default early edge. To be
  calibrated on the Jamsin et al. (arXiv:2610.28775) injections; the
  light-curve fit (NMMA/Redback t0 samples, `t0_samples_gps`) should replace
  it when available.
- **BNS ranges**: defaults H1 150, L1 160, V1 55 Mpc are approximate O4 values;
  pass the live ranges and science segments of the window. Horizon = 2.264 × range.
- **Thresholds**: targeted network SNR 6, P_min 0.1, hard cut 500 Mpc.
- **Prompt GRB → θ_obs ≲ 30°** soft cut; jet-structure model would be better.
- **Redshift → distance** at fixed Planck18 cosmology with 300 km/s peculiar
  velocity (do not use it for an H0 measurement).
- **Kilonova candidate rate** R_KN = 30/month (rate of filter candidates, incl.
  contaminants): set from the actual broker filter.
- **Distance overlap** referenced to the network horizon; meaningful only for
  GW maps with distance layers.
- Fink LSST endpoint/schema not checked; ZTF API column names used.

## Not yet done

- Targeted-search runner (PyCBC/PyGRB config from the request JSON, sub-segmentation).
- Archive queries for the last non-detection in other surveys (ATLAS, ZTF forced photometry, GOTO, BlackGEM).
- Model-systematics: a t0 bias like case A would make a targeted search miss the merger; try Bu2026 (Jamsin et al.) / other models, and widen windows by a model-systematic term.
- GCN Kafka listener and live GraceDB/segment (DQSegDB) queries.

## License and credits

GPL-3.0-or-later (see `LICENSE`), as LVK RAVEN (`ligo-raven`) and NMMA, which this code follows or borrows from:
the built-in Bu2019lm / Bu2023Ye model priors are copied from NMMA's prior files, and
`examples/validation/AT2017gfo_reduced.dat` is NMMA's `example_files/lightcurves/AT2017gfo_reduced.dat`
(https://github.com/nuclear-multimessenger-astronomy/nmma, tag v1.0.2).
Merger-time methodology: Jamsin et al., arXiv:2610.28775.
