# Light search around O4c sub-solar-mass candidates

Scripts that re-examine three public O4c superevents with a sub-solar-mass component
(S250720j, S250818k, S251112cm) with the light targeted search, at a fixed sky position, on
LVK data at CIT. They need LVK credentials and data access; no strain, pipeline trigger
parameter or result is distributed here.

1. `cp events.example.json events.json` and fill, for each event, the template of the pipeline
   event to re-examine (detector-frame `mchirp`, `m1`, `m2`, `spin1z`, `spin2z`, from GraceDB).
   The time, window and sky position given are public.
2. `python prepare.py` (CIT login node, IGWN environment, SciToken): science segments and
   about 17 minutes of strain per event, written to small frames, for two kinds of data:
   - `raw`: LIGO RAW frames (`CAL-DELTAL_EXTERNAL_DQ` converted to strain) + Virgo `HoftOnline`;
   - `hoft`: LIGO calibrated `GDS-CALIB_STRAIN_CLEAN` + Virgo `HoftOnline`.
3. `condor_submit light_search.sub` (a copy or link of the `raven_kn` package must be in the
   directory): 12 jobs, `light_search.sh <event> <raw|hoft> <tmpl|equal> <coherent|network>`.
   - `tmpl`: templates with the mass ratio and aligned spins of the pipeline event
     (`--mass-ratio`, `--spin1z`, `--spin2z`), chirp mass within 1 %;
   - `equal`: equal-mass non-spinning templates (the default bank), chirp mass within 1 %.
   Each job runs 20 off-source trials.
   `condor_submit light_search_polarisation.sub` adds, on the calibrated strain with the pipeline
   template, the coherent statistic restricted to a circular polarisation (`circular`, `face_on`,
   `face_away`), to compare with the unrestricted one.
4. `python summary.py`: one line per result (statistic, per-detector SNR, chirp mass, time
   offset from `t0`, off-source maximum and p-value).
