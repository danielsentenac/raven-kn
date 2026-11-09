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
| `raven_kn.connection` | `KNCandidate` (photometry incl. non-detections, redshift, distance, viewing angle, GRB, light-curve-fit samples, spectroscopic class); connectors: Fink, ALeRCE, ANTARES, ATLAS, GraceDB, GCN JSON, NED/NED-D, GLADE+, DES Y6, TNS, local JSON | `gracedb_events.ExtTrig` / `SE` |
| `raven_kn.treatment.t0` | merger-time window [tl, th]: GRB (seconds) > light-curve-fit t0 posterior ∩ photometric bracket (last *constraining* non-detection → first detection) | fixed `tl, th` |
| `raven_kn.treatment.lightcurve_fit` | NMMA wrapper: candidate → NMMA data/prior files and `lightcurve-analysis` command → posterior read back as `t0_samples_gps` (+ D_L, ι) | — (new) |
| `raven_kn.treatment.distance_inclination` | joint (D_L, cos ι) posterior from spectroscopic redshift, host distance, VLBI/afterglow viewing angle, prompt GRB, light-curve-fit samples; search priors | — (new) |
| `raven_kn.treatment.detectability` | expected network SNR over the window (LAL antenna patterns at the candidate position, Earth rotation, duty cycle/segments); decision: trigger or not | — (new) |
| `raven_kn.treatment.far` | joint FARs: `(th−tl) R_KN FAR_GW` / sky overlap / **distance overlap**; targeted-search p-value, sub-segment trials, RAVEN "SubGRBTargeted" product formula | `search.coinc_far` |
| `raven_kn.skymap` | kilonova 2D/3D HEALPix map (PROB, DISTMU, DISTSIGMA, DISTNORM), GW map reader, sky overlap (RAVEN's own function if installed), distance overlap | `search.skymap_overlap_integral` |
| `raven_kn.listener` | GCN Kafka listener: `igwn.gwalert` (multi-order sky map decoded in the alert) and unified-schema GRB notices matched against a store of kilonova candidates, both ways in time | `listen` / `handle_*` |
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

## Validation with real NMMA (v1.0.2, sklearn_gp, pymultinest 512 live points)

Two real kilonovae with a known merger time, the prompt signal hidden from the fit:
AT2017gfo (`examples/validate_at2017gfo.py`, PS1 g,r,i,z, first point +17 h; truth GW170817) and
GRB 211211A (`examples/validate_grb211211a.py`, PS1/Bessell/2MASS, first point +6.2 h; truth the
Swift-BAT/Fermi-GBM trigger 2021-12-11T13:09:59, GCN 31202). Distance from the host redshift
(and VLBI angle for AT2017gfo) as priors, priors clipped to the model training range. ~4-10 min per run.

| event | model, data fitted | t0 median − true | 90 % interval | truth in 90 % | window with 0.75 margin | truth inside |
|---|---|---|---|---|---|---|
| AT2017gfo | Bu2019lm, all data | −40.1 h | 32 h | no | 89 h | yes |
| AT2017gfo | Bu2019lm, first 3 d | −18.4 h | 46 h | no | 90 h | yes |
| AT2017gfo | Ka2017, first 3 d | −27.4 h | 87 h | no | 138 h | yes |
| AT2017gfo + g,r non-detections 3 h before merger (synthetic) | Bu2019lm, all data | −1.3 h | 9 h | yes | 19 h | yes |
| GRB 211211A | Bu2019lm, all data | −28.2 h | 25 h | no | 74 h | yes |
| GRB 211211A | Bu2019lm, first 3 d | −31.8 h | 26 h | no | 78 h | yes |
| GRB 211211A | Ka2017, first 3 d | −26.2 h | 38 h | no | 83 h | yes |
| GRB 211211A, no distance prior | Bu2019lm, all data | −25.9 h | 24 h | no | — | — |

Findings:
- Without a recent pre-merger non-detection, both models put the merger **~1 day too early** and the
  statistical interval misses it; the model-systematic margin (default 0.75) is what recovers it.
  A non-detection a few hours before the merger removes the bias (−1.3 h) and gives a 19 h window.
- Bu2019lm cannot follow AT2017gfo's slow evolution: at the true t0 it is up to ~1 mag too bright in i,z at
  1.5–2.5 d and ~1 mag too faint after 6 d (chi2 profile over fixed t0); the fit makes the source older,
  with parameters at the grid edges. Fitting only the first 3 days (`fit_days`, default) halves the bias
  on AT2017gfo but not on GRB 211211A, whose early optical light includes afterglow.
- The kilonova-only distance is unreliable: 165 Mpc (90 %: 115–274) for GRB 211211A at 357 Mpc.
  Use a spectroscopic redshift for the "is a search worth it" decision.
- Fixed on the way: priors copied from NMMA's `Bu2019lm.prior` exceeded the model's training range
  (50–80 % of samples extrapolated); NMMA crashes on non-detection-only filters (`nmma_launcher.py`
  works around it); the timeshift bound now follows the distance prior's redshift.

Reports: `examples/validation/report*.json`, `examples/validation_grb211211a/report*.json`.

### Injection calibration of the t0 window (`examples/injections/calibrate_margin.py`)

29 synthetic kilonovae (Bu2019lm and Ka2017, parameters inside the training ranges, 40–150 Mpc uniform in
volume, known host redshift), observed with a ZTF-like survey (nightly g/r/i epochs at a random phase,
5σ depth 21.5, 30 % of nights lost, pre-merger coverage in half of the cases), each fitted with both models
(58 NMMA fits, `--nlive 512`, first 3 days):

| fits | statistical 90 % interval covers t0 | with margin 0.5 / 0.75 / 1.0 | final window (fit ∩ bracket), 0.75: coverage, median width |
|---|---|---|---|
| same model as the truth | 29/29 | 29/29 | 29/29, 59 h |
| other model | 23/29 (79 %) | 28/29 / 28/29 / 29/29 | 28/29, 59 h |
| *Rubin-like (`--survey rubin`: depth 24.5, revisit 3 d, 2 of g/r/i, 40–300 Mpc, fit 6.5 d), 27 injections:* | | | |
| same model as the truth | 26/27 | 27/27 | 27/27, 83 h |
| other model | 19/27 (70 %) | 27/27 / 27/27 / 27/27 | 27/27, 83 h |
| *ZTF + Rubin on the same field (`--survey mixed`: both cadences, independent phases and lost nights, each limit at its own depth, 40–300 Mpc, fit 6.5 d), 30 injections:* | | | |
| same model as the truth | 26/30 | 30/30 / 30/30 / 30/30 | 30/30, 66 h |
| other model | 20/30 (67 %) | 28/30 / 30/30 / 30/30 | 30/30, 101 h |

- Self-model posteriors are well calibrated; the margin is a model-systematic correction, as intended. The
  one cross-model case needing more than 0.75 (0.86) is a Ka2017 kilonova fitted with Bu2019lm, first seen
  at +31 h with no pre-merger coverage.
- The injections exposed a **bracket bug**: a g non-detection minutes before the first r detection of a
  red kilonova was taken as "before the merger", and two faint kilonovae were not yet at M = −15 after 6 h;
  9/29 brackets excluded the true merger, by up to 65 h. Fixed (same-night rule, 0.5 d), 29/29 now.
- At ZTF depth, the fit with its margin is about as wide as the bracket (median 59 h); without margin it
  would be 36 h but miss 10 % of mergers.
- At Rubin depth the margin 0.75 covers everything (largest needed fraction 0.37), but the window is
  *wider* than at ZTF depth (83 h): with a 3-day revisit the last non-detection is typically 3 days before
  discovery. Depth helps only with recency: a narrow window needs a deep limit hours before the merger
  (Jamsin et al.), i.e. Rubin combined with nightly surveys (ZTF, ATLAS) or deep-drilling fields.
- Rubin + ZTF confirms it. Median window (margin 0.75) when the field was observed before the merger:
  ZTF alone 35 h (but only to 150 Mpc), Rubin alone 83 h, **both 50 h** (first detection at +20 h instead
  of +41 h); for late discoveries all three stay at 130–140 h. Margin 0.75 still covers 60/60 fits (largest
  needed fraction 0.62, cross-model); shallow ZTF limits taken after the merger, before Rubin's first
  detection, are correctly ignored by the bracket (not deeper than the first detection).

Reports: `examples/injections/runs_inj/report.json` (ZTF-like), `runs_inj_rubin/report.json` (Rubin-like),
`runs_inj_mixed/report.json` (ZTF + Rubin).

Environment notes: NMMA runs from the conda env `nmma` (`--executable` or `$RAVEN_KN_NMMA_EXE`);
models in `~/nmma_models` (`python -c "from nmma.core.gitlab import main; main()" --model Bu2019lm --filters ...`,
the `svdmodel-download` entry point is broken in v1.0.2); some sncosmo bandpasses redirect to SVO's dead
`http://` URL and had to be fetched once from `https://svo2.cab.inta-csic.es/theory/fps/` into sncosmo's cache.

## Targeted GW search (`raven_kn.targeted`, `raven-kn targeted`)

Light PyCBC-based search at the kilonova position over its t0 window, following Jamsin et al.:
- template bank: equal-mass, non-spinning TaylorF2 over a chirp-mass range, spaced with the analytic
  chirp-mass metric calibrated on direct matches (minimal match 0.97: 181 templates for 1.18–1.22 Msun,
  ~2150 for 1.0–1.5, ~10000 for sub-solar 0.4–0.9);
- per sub-segment (default 1024 s): GWOSC, local hdf5 or frame files (`IFO=path:channel`), gating of loud
  glitches, median-Welch PSD, matched filter;
- signal consistency: 16-bin power chi-square where a detector's SNR > 5, SNR replaced by newSNR;
- each detector aligned by its time delay for the known sky position; network SNR sqrt(sum rho_d^2)
  counted only where >= 2 detectors exceed 4 (coincidence), maximised over templates;
- background: time slides (independent shifts per detector, many per sub-segment); the loudest on-source
  candidate is removed from the background only if loud in two detectors (2nd-loudest SNR >= 6), so a
  single-detector glitch stays in its own background; p-value, exponential-tail extrapolation flagged as
  formal when far beyond the sampled background; sub-segments combined with a trials factor.
- off-source trials (`--n-offsource N`): the same search, zero-lag, in N windows of the sub-segment length
  outside the whole on-source window (alternately before/after, 64 s gap), as PyGRB does; it samples the
  real non-stationarity around the window that slides of on-source data cannot. With trials, each
  sub-segment uses the larger of the slide and off-source p-values.

```bash
raven-kn targeted --request out/<name>_targeted_search_request.json --out targeted_out \
    --ca-bundle /etc/pki/tls/certs/ca-bundle.crt      # if Python's SSL fails against gwosc.org
```

GW170817, H1 and V1 from GWOSC, L1 from the glitch-subtracted release (LIGO-P1700349,
`L-L1_LOSC_CLN_4_V1-1187007040-2048.gwf`); 181 templates, 1000 slides:

| window | network SNR | H1 / L1 / V1 | time, chirp mass | p-value | background max / median |
|---|---|---|---|---|---|
| containing the merger | 30.9 (published 32.4) | 18.6 / 24.5 / 3.2 | 1187008882.426, 1.1976 | < 1/1001 (22 background spreads beyond the loudest slide) | 8.5 / 7.5 |
| 15 min earlier (noise) | 7.4 | 4.0 / 5.4 / 3.1 | — | 0.62 | 8.5 / 7.5 |

With the GWOSC GWTC-1 L1 data instead (`GW170817_signal/noise_targeted_result.json`), the glitch 1.1 s
before the merger is gated out with part of the signal and the chi-square penalises the damaged L1 signal:
network SNR 22.6 (H1 / L1 / V1 = 18.5 / 12.6 / 3.2), still p < 1/1001 against a background of max 8.3;
the noise window gives p = 0.48. Before hierarchical removal was added, the H1 signal leaked into every
slide (background median 18.7). Tests use synthetic Gaussian noise
with injections and a sine-Gaussian glitch. Results: `examples/targeted_gw170817/`.

**Coherent statistic** (`--statistic coherent`, as PyGRB): templates shifted by each detector's exact
delay, complex SNRs projected on the network's response to the two polarisations at the known sky position
(coherent SNR), the remainder is the null stream, cut as in PyGRB (`null_min` 5.25, `null_grad` 0.2,
`null_step` 20). Same GWOSC data, 181 templates, 1000 slides:

| window | coherent SNR (null) | network SNR | p-value | background max / median |
|---|---|---|---|---|
| containing the merger | 22.4 (3.1) at 1187008882.429, Mc 1.1976 | 22.6 | < 1/1001 | 7.5 / 6.6 (network: 8.3 / 7.5) |
| 15 min earlier (noise) | 6.6 (2.1) | 7.5 | 0.34 (network 0.48) | 7.5 / 6.5 (network: 8.8 / 7.5) |

The signal keeps its SNR while the noise background drops by about one unit (null cut, exact-sample
coincidence instead of 8 ms pooling). On synthetic data a 3-detector injection gives coherent SNR 29.5 with
null 1.1, and the same injection with V1's sign flipped (loud everywhere, but no source at that position
gives those phases) falls from 57 to 8 through the null cut. Writing these tests showed that the synthetic
injections had used degrees where `project_wave` expects radians (wrong injected sky position, hidden by the
incoherent statistic's pooling); fixed.
Remaining limitations: incoherent statistic, zero-spin equal-mass templates. Production: PyGRB.

## Redshift, distance and inclination: spectroscopy and host galaxies (`raven-kn host`, `raven-kn tns`)

The "is a search worth it" cut (500 Mpc) and the GW distance/inclination priors rest on the EM distance.
A kilonova spectrum rarely gives a precise redshift by itself (broad features); the redshift comes from the
host galaxy's spectrum. The connection module now fetches it:

| source | gives | access |
|---|---|---|
| **NED** (TAP `NEDTAP.objdir`) | galaxies with redshift around the candidate: SDSS, DESI, 6dF, literature; spectroscopic/photometric flag | public |
| **NED-D** | redshift-independent distances of the host (SBF, TRGB, Cepheids, SNe Ia, FP, Tully-Fisher); GW-based entries are dropped | public |
| **GLADE+** (VizieR VII/291, mirrors with fallback) | LVK's galaxy catalogue: CMB-frame z (spec. or photometric), B/K magnitudes, chance-coincidence probability | public |
| **DES Y6 Gold** (NOIRLab Data Lab TAP `des_dr2.y6_gold`) | southern 5,000 deg², ~450 M galaxies (20× GLADE+ there, complete well beyond 150 Mpc): DNF photo-z, r magnitude, half-light radius | public |
| **TNS** | spectroscopic classification (e.g. "Kilonova"), transient/host redshift, list of spectra (ASCII/FITS) | bot ID + API key |

`find_host` keeps a spectroscopic redshift already in the candidate (from TNS, a GCN circular or a spectrum);
otherwise takes NED's nearest galaxy with a spectroscopic redshift within 50 kpc projected, else GLADE+ (lowest
chance-coincidence probability, may be photometric), else DES Y6 Gold in its footprint (photometric); for a
host at z < 0.03 it adds the best NED-D distance.
Redshifts are converted to the CMB frame; the 300 km/s peculiar-velocity error is added downstream.

| event | result | truth |
|---|---|---|
| AT2017gfo | NGC 4993 (NED), z_cmb 0.0109 spec., 2.2 kpc; SBF 40.7 ± 1.5 Mpc (NED-D) | 40.7 Mpc (SBF) |
| GRB 211211A | WISEA J140910.48+275320.6 (NED), z 0.0763 spec. | z 0.0763 (not in GLADE+) |
| AT2025ulz | GLADE+ 15877920, 0.9″, **photometric** z 0.102 ± 0.036 (NED has no host redshift yet) | DESI/Keck z 0.0848 |

DES Y6 DNF photo-z against 109 NED spectroscopic galaxies at 0.01 < z < 0.15 (7 DES fields): NMAD
Δz/(1+z) = 0.013, bias < 0.007, 1% outliers; the catalogue's `dnf_zsigma` is 1.2–1.5× too small, so the error
used is 1.5 × `dnf_zsigma` (floor 0.01). ESO 243-49 (z 0.0224), transient placed 8″ off: DES host at 7.5″,
photo-z 0.022 ± 0.013, P_cc 8e-4 (P_cc in r with the half-light radius, Bloom et al. 2002).

So: NED first (spectroscopic completeness beyond 100 Mpc), GLADE+, then DES Y6 as fallbacks; for a new event the
spectroscopic host redshift typically arrives hours to days later through TNS or GCN circulars, and replaces
a photometric one.

**Inclination.** Spectra constrain the viewing angle only through radiative-transfer models (blue/red
component ratio, line velocities, polarisation), and model-dependently. The constraints that matter come from
the afterglow (VLBI centroid motion, light-curve modelling: GW170817 θ_obs ≈ 15–25°) or from kilonova
light-curve/spectral fits. RAVEN-KN takes them as `ViewingAngleMeasurement(mean, std, method)`; nothing is
extracted automatically from spectra.

## Last non-detection from survey archives

The t0 window depends most on the latest deep image of the position before the first detection.
`connection.py` merges photometry from several sources into one candidate:
- `ALeRCEConnector` — ZTF detections and non-detections (`diffmaglim`) from the public ALeRCE API
  (tested live: ZTF21abdpqpq, last non-detection r > 20.5 1.98 d before the first detection);
- `ANTARESConnector` — the NOIRLab ANTARES broker (public, no token), found by cone search: the same ZTF
  detections and limits as ALeRCE (ZTF21abdpqpq: identical 3 + 91 points, merged without duplicates), so
  either can stand in when the other is down; LSST records are read from ANTARES' generic fields;
- `ATLASForcedPhotometryConnector` — ATLAS forced photometry (cyan/orange, 5σ limits) via the
  fallingstar-data.com queue API; needs a personal token (written from the API guide, tested on mocked
  responses only);
- `FinkConnector` (ZTF API);
- `FinkLSSTConnector` — Rubin light curves from the Fink/LSST API (`api.lsst.fink-portal.org`, checked live
  2026-11): alert sources and forced photometry (DiaForcedSource) in nJy, converted to AB magnitudes
  (S/N >= 5) or 5σ limits, TAI -> UTC; `conesearch` for Rubin objects at a known position. Fink's forced
  photometry currently covers only the visits of the object's own detections, so pre-discovery Rubin limits
  are not yet available there. NMMA fits of Rubin data need the `lsst*` filter models (not downloaded here);
- `last_nondetection(kn, depth_mag)` — latest (optionally deep enough) limit before the first detection,
  over all merged surveys.

## PyGRB configuration (`raven-kn pygrb`)

The production targeted search is PyGRB (coherent SNR, chi-square and null-stream vetoes). It is built for
on-source windows of seconds to minutes, so the kilonova window is split into sub-segments; per
sub-segment RAVEN-KN writes an event file and a `pycbc_make_offline_grb_workflow` command
(`run_pygrb.sh`), with one shared HDF template bank. PyCBC's config parser rejects an option defined in two
files, so:
- the event file *replaces* the per-GRB file of the LVK base configuration: it is built from a template
  (e.g. LVK's `GRB170817A.ini`) whose `[workflow]` event options (trigger name/time, ra/dec as
  `"<value> rad"`, sky-error, start/end time) are replaced and whose other sections (injection settings)
  are kept, self-references pointed to the new file;
- on-/off-source durations and the bank are passed as `--config-overrides`
  (`workflow-exttrig_segments:on-before/on-after/min-duration/max-duration`,
  `workflow-tmpltbank:tmpltbank-pregenerated-bank`).

```bash
raven-kn pygrb --request out/<name>_targeted_search_request.json --out pygrb_configs \
    --base-config analysis.ini executables_common.ini executables_condorpool.ini injections.ini \
                  postprocessing.ini cats.ini o2-data.ini \
    --event-template GRB170817A.ini --bank-path /path/on/cluster/<name>_bank.hdf
```

**Tested on CIT** (2026-11-06) with LVK's O2 test configuration (pycbc-config `O4/pygrb/dev/o2_tests`,
IGWN conda PyCBC 2.8.1; the CVMFS RHEL-8 PyCBC virtualenv does not run on CIT's AlmaLinux 10): a 1024 s
on-source window around GW170817 with the RAVEN-KN bank generates a complete workflow (29 jobs: coherent
H1/L1/V1 matched filtering on O2 cleaned C02 frames over 4 bank splits, injections, clustering, results
page). Not yet submitted. Lessons: PyGRB aborts if the on-source contains vetoed time (a CAT2 H1 second
444 s after the merger stopped the first attempt), so for long windows pass the vetoed intervals
(`--vetoes FILE`, all detectors merged): sub-segments are then placed in the veto-free stretches
(`veto_free_subsegments`, 1 s guard, pieces < 16 s dropped), or let RAVEN-KN query them:
`--veto-definer VDF.xml --veto-categories 1 2 --science-flag H1=H1:DMT-ANALYSIS_READY:1 ...` (gwpy +
DQSegDB, `raven_kn.targeted.segments`; padding as in the veto definer; tested offline, and on the O2
veto definer file at CIT — DQSegDB was unreachable on 2026-11-09); the `--output` directory must exist before generation (`run_pygrb.sh` creates it).

## Real-time listener (`raven_kn.listener`, `raven-kn listen`)

```
export GCN_CLIENT_ID=... GCN_CLIENT_SECRET=...          # from gcn.nasa.gov; pip install gcn-kafka
raven-kn listen --store kn_store/                        # igwn.gwalert + GRB topics
raven-kn add candidate.json --store kn_store/            # a new or updated kilonova candidate
```

- **GW alert** (`igwn.gwalert`): the base64 multi-order FITS map is decoded (UNIQ → order, pixel;
  flattened to NESTED nside 512 for the overlaps, exact row looked up at the candidate). Every stored
  candidate whose t0 window contains the merger time gets the joint FAR (temporal / sky / distance), written
  to `kn_store/coinc/<superevent>__<candidate>.json` with the GW conditional-distance parameters at the
  candidate. Updates replace the stored alert (an alert without a map keeps the previous map); retractions
  delete it.
- **GRB notice**, GCN unified JSON (`trigger_time`, `ra`, `dec`, `ra_dec_error`: Swift-BAT GUANO,
  Einstein Probe WXT, SVOM ECLAIRs) or classic VOEvent XML (Fermi-GBM flight/ground/final positions, with
  the 3.7° GBM systematic added in quadrature): attached to candidates within 2× the error radius whose
  first detection follows within 3 days; their window becomes RAVEN's [−5, +1] s. With 1× the radius the
  GW170817 GBM ground position (24° from NGC 4993) would be missed. Test notices are ignored.
- **Late candidate**: `raven-kn add` matches it against GW alerts of the last 14 days.
- **LVK RAVEN coincidence**: when an alert carries `external_coinc` (LVK's own GW + GRB coincidence:
  observatory, time difference, temporal and sky-position joint FARs, combined GW × GRB map), it is kept
  with the superevent, also across later updates that omit it, and copied into every RAVEN-KN result as
  `raven_external_coinc`, with the combined map at the candidate (`combined_map_at_kn`). The sky factor
  uses the GW-only map unless `--use-combined-skymap` is given.

Nothing is uploaded to GraceDB from the listener (use `calc_signif_gracedb` once results are vetted).
Tested on synthetic alerts only (`tests/test_listener.py`); not yet run against the live Kafka stream.

## Real case: S250818k / AT2025ulz

`examples/s250818k_case.py` runs RAVEN-KN on the public S250818k sky maps (GraceDB; FAR 6.8e-8 Hz,
HasSSM > 99 %) and AT2025ulz (GCN 41414; host z = 0.0848, DESI). Pre-discovery ZTF limits were shallower
than the transient, so the t0 window is the 7-day fallback.

- Host distance: 400 Mpc (Planck18; 370–386 Mpc for H0 = 73–70), ±5 Mpc from peculiar velocity.
- GW distance along the host direction (conditional posterior from the 3D map, `conditional_distance`):
  272 ± 74 Mpc, 90 % 153–395 Mpc (BAYESTAR); 266 ± 70 Mpc, 90 % 153–383 Mpc (Bilby); P(D > 400 Mpc) = 3–4 %,
  i.e. a ~1.8–2σ tension. (DISTMU/DISTSIGMA, 218 / ~80 Mpc, are ansatz parameters, not mean and sd.)
- Joint FAR: 14.8/yr (time), 0.22/yr (time + sky), 1.8–2.5/yr with the distance factor. The host lies beyond
  the reference horizon (362 Mpc BNS, 172 Mpc for Mc = 0.5), so the distance factor is only approximate
  (`distance_overlap_valid` = False); the detectability argument below is the cleaner statement.
- Expected network SNR at 400 Mpc: 3.3 for a BNS, 1.6–2.5 for chirp masses 0.5–0.87 Msun (approximate O4c
  ranges): no targeted search would be triggered.
- Not done: GW parameter estimation with sky position and distance fixed to the host (needs O4c strain data
  or LVK posterior samples), the proper test of whether the GW data accommodate 400 Mpc.

SN 2025ulz was classified as a type IIb supernova; the physical link is still debated.

## Assumptions to revisit (all are parameters)

- **Model-systematic margin** on a light-curve t0 (`window_from_samples`, `--model-sys-frac`, default 0.75):
  the 90 % interval is widened by 0.75 x the extrapolation gap from the median t0 to the data that anchor it
  (first detection on the late side; last constraining non-detection, else the same gap, on the early side),
  at least 1 h, then intersected with the photometric bracket. 0.75 is the smallest round value recovering
  the merger in all the real-data fits above (two events, two models), and covers 97 % of cross-model
  injections (below); 1.0 covers all of them.
- **Bracket model** (`window_from_bracket`): kilonova at least as bright as
  M = −15 from 0.5 d after merger; a non-detection constrains only if deeper
  than that at the 95th-percentile distance; a non-detection in another band
  less than 0.25 d before the first detection is ignored (colour, not age);
  7 d default early edge. Calibrated on the injections below (it was 0.25 d
  and had no same-night rule: 9/29 brackets excluded the merger).
- **BNS ranges**: defaults H1 150, L1 160, V1 55 Mpc are approximate O4 values;
  pass the live ranges and science segments of the window. Horizon = 2.264 × range.
- **Thresholds**: targeted network SNR 6, P_min 0.1, hard cut 500 Mpc.
- **Source class**: ranges are BNS (1.4+1.4 Msun) ranges; `--chirp-mass` / `chirp_mass_msun` rescales them
  as Mc^(5/6) (inspiral-dominated sources only), e.g. to the chirp-mass bin of a sub-solar (HasSSM) alert.
  The distance factor of the joint FAR uses the same scaled horizon.
- **Prompt GRB → θ_obs ≲ 30°** soft cut; jet-structure model would be better.
- **Redshift → distance** at fixed Planck18 cosmology with 300 km/s peculiar
  velocity (do not use it for an H0 measurement).
- **Kilonova candidate rate** R_KN, the analogue of RAVEN's GRB rate (`candidate_rate`): RAVEN-KN applies no
  kilonova classifier itself; it takes what upstream filters pass. R_KN = 5/yr for a spectroscopically
  classified kilonova (`classification`, e.g. from TNS); else, with a classifier score `kn_score` and the
  filter's score distribution (`score_survival`), 30/month x the fraction of filter candidates scoring at
  least as high (floor 1/yr); else 30/month, the rate of everything the filter passes (incl. contaminants).
  Both numbers are placeholders to set from the actual filter and from TNS statistics. A spectroscopic
  classification lowers the joint FAR by 72x.
- **Distance overlap** referenced to the network horizon; meaningful only for
  GW maps with distance layers.

## Not yet done

- Compare the PyGRB result with the light search; live DQSegDB test of the veto query.
- Coherent statistic and null stream are left to PyGRB: the light search is a fast first pass on coincident network SNR (a coherent statistic through the time slides would keep the complex SNR of every template and detector, ~9 GB for the GW170817 bank). PyGRB run on GW170817 at CC-IN2P3 (Slurm) pending.
- Injection calibration with mixed surveys (Rubin + nightly ZTF/ATLAS) and a third model (e.g. Bu2026).
- ZTF forced photometry (IRSA ZFPS), GOTO, BlackGEM archive connectors; live test of the ATLAS connector.
- Listener: run on the live GCN stream (needs GCN credentials); live GraceDB/segment (DQSegDB) queries.

## License and credits

GPL-3.0-or-later (see `LICENSE`), as LVK RAVEN (`ligo-raven`) and NMMA, which this code follows or borrows from:
the built-in Bu2019lm / Ka2017 / Bu2023Ye model priors are adapted from NMMA's prior files, and
`examples/validation/AT2017gfo_reduced.dat` and `examples/validation_grb211211a/GRB211211A.txt` are NMMA's
`example_files/lightcurves/AT2017gfo_reduced.dat` and `GRB211211A.txt`
(https://github.com/nuclear-multimessenger-astronomy/nmma, tag v1.0.2).
Merger-time methodology: Jamsin et al., arXiv:2610.28775.
