# Ongoing work (state on 2026-11-18)

Goal: use what the electromagnetic alerts give (inclination, distance, redshift) as priors that restrict
the targeted GW search and the parameter estimation of its candidate.

## Done and pushed

| Prior | What it does | Validated on GW170817 |
|---|---|---|
| Inclination | coherent statistic restricted to a circular polarisation (`--polarisation`, selected by the request when the viewing angle is below 60°) | yes, `docs/examples/gw170817-search.md` |
| Distance | candidates ranked by the likelihood marginalised over the EM distance and inclination (`--distance-prior`, opt-in) | yes, same page: statistic 29.79 against 30.07, no off-source candidate in 100 trials |
| Redshift | chirp-mass range of the request in the detector frame, source range times (1 + z) | unit tests only |
| Parameter estimation | Bilby on the loudest candidate, sky fixed, with or without the EM priors (`--pe em\|gw\|both`) | **no, see below** |

## Running now

Parameter estimation of GW170817 on CIT, directory `~/raven_kn_em_priors`, two HTCondor jobs
(`pe.sub`: `em` and `gw`), started 2026-11-18. The first run took 2 h and 3.7 h.

- The first run (kept in `~/raven_kn_em_priors/first_run`) is wrong: both distance posteriors sit on the
  lower prior bound and the log Bayes factors (664, 13942) exceed what an SNR of 32 allows (about 525).
  Bilby's optimiser of the relative-binning reference point diverged.
- Fix in `raven_kn/targeted/pe.py`: `maximise()` finds the reference point with the exact likelihood
  (analytic over distance and phase, FFT over time) and Bilby's optimiser is off. Tested on a synthetic
  injection only (`tests/test_pe.py`, 4 tests pass; full suite not rerun after the fix).
- The copy of `raven_kn` in `~/raven_kn_em_priors` has the fix.

## Next steps

1. Read the two summaries when the jobs end: `SSS17a_pe_em_pe_summary.json`, `SSS17a_pe_gw_pe_summary.json`.
   Check before trusting them:
   - `maximum_likelihood.snr` close to the search value (about 30) and a distance near 40 Mpc;
   - `distance_at_prior_edge` is null;
   - log Bayes factor below about 525.
2. If they are right: add the with / without EM priors comparison (distance, viewing angle) to
   `docs/examples/gw170817-search.md`, copy the three `search_*_targeted_result.json` and the two PE
   summaries to `examples/em_priors_gw170817/`, run the full test suite, commit.
3. If they are wrong again: compare the exact and the relative-binning likelihood at the reference point
   on the real data (`maximise()` returns both), then look at the PSD (`_interferometers`).
4. Open points of the distance prior:
   - the face-on distance estimate of GW170817 is 54 Mpc, at the upper edge of the 30–54 Mpc prior; check
     whether the estimate is biased high;
   - the sensitivity gain is not measured (empty background, p-value at its floor 1/101): needs
     injections of weak signals at the prior distance;
   - the request of the example has no joint distance-inclination samples file, so the runs used the
     intervals; the samples path (`prior_samples`) is untested on real data;
   - decide whether the distance prior becomes automatic when the distance comes from a spectroscopic
     redshift.
5. Sub-solar-mass O4c events: the polarisation runs are finished on CIT (`~/raven_kn_ssm_o4c`). Their
   results and the pipeline template values are not public and stay out of this repository.

## How to run on CIT

- Search and PE jobs run on CIT with HTCondor, not on the workstation.
- `ssh cit` goes through `ssh.igwn.org`; the connection drops after some hours, reopen it with
  `ssh -fN cit` (`ssh-add ~/.ssh/id_rsa` first if the agent refuses the key).
- Example to follow: `examples/em_priors_gw170817/` (`prepare.sh`, then `condor_submit jobs.sub`, then
  `condor_submit pe.sub`). Job logs appear only when a job ends.
- To update the code there: `rsync -a --exclude __pycache__ raven_kn cit:raven_kn_em_priors/`.
