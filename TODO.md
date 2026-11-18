# Ongoing work (state on 2026-11-18)

Goal: use what the electromagnetic alerts give (inclination, distance, redshift) as priors that restrict
the targeted GW search and the parameter estimation of its candidate.

## Done and pushed

| Prior | What it does | Validated on GW170817 |
|---|---|---|
| Inclination | coherent statistic restricted to a circular polarisation (`--polarisation`, selected by the request when the viewing angle is below 60°) | yes, `docs/examples/gw170817-search.md` |
| Distance | candidates ranked by the likelihood marginalised over the EM distance and inclination (`--distance-prior`, opt-in) | yes, same page: statistic 29.79 against 30.07, no off-source candidate in 100 trials |
| Redshift | chirp-mass range of the request in the detector frame, source range times (1 + z) | unit tests only |
| Parameter estimation | Bilby on the loudest candidate, sky fixed, with or without the EM priors (`--pe em\|gw\|both`) | yes, same page: 40.6 (32.1–46.9) Mpc with the EM priors, 39.6 (26.7–46.6) Mpc without |

## Parameter estimation: what was wrong and how it was fixed

- First run (kept on CIT in `~/raven_kn_em_priors/first_run`): both distance posteriors on the lower prior
  bound and log Bayes factors (664, 13942) above what the SNR allows. Bilby's optimiser of the
  relative-binning reference point diverged.
- Fix in `raven_kn/targeted/pe.py`: `maximise()` finds the reference point with the exact likelihood
  (analytic over distance and phase, FFT over time) and Bilby's optimiser is off.
- Second run (2026-11-18, 3.6 h and 4.7 h): maximum-likelihood SNR 33.5, log Bayes factors 537.7 and
  532.5 below SNR²/2 = 562.6, `distance_at_prior_edge` null in both. Full test suite passes (132 tests).

## Next steps

1. Open points of the parameter estimation:
   - the reference point of `maximise()` is at 25 Mpc, viewing angle 59° and mass ratio 1.0, the inclined
     end of the degeneracy and the mass-ratio bound: not studied further;
   - the gain from the EM priors is modest with intervals; redo with joint distance-inclination samples.
2. Open points of the distance prior:
   - the face-on distance estimate of GW170817 is 54 Mpc, at the upper edge of the 30–54 Mpc prior; check
     whether the estimate is biased high;
   - the sensitivity gain is not measured (empty background, p-value at its floor 1/101): needs
     injections of weak signals at the prior distance;
   - the request of the example has no joint distance-inclination samples file, so the runs used the
     intervals; the samples path (`prior_samples`) is untested on real data;
   - decide whether the distance prior becomes automatic when the distance comes from a spectroscopic
     redshift.
3. Sub-solar-mass O4c events: the polarisation runs are finished on CIT (`~/raven_kn_ssm_o4c`). Their
   results and the pipeline template values are not public and stay out of this repository.

## How to run on CIT

- Search and PE jobs run on CIT with HTCondor, not on the workstation.
- `ssh cit` goes through `ssh.igwn.org`; the connection drops after some hours, reopen it with
  `ssh -fN cit` (`ssh-add ~/.ssh/id_rsa` first if the agent refuses the key).
- Example to follow: `examples/em_priors_gw170817/` (`prepare.sh`, then `condor_submit jobs.sub`, then
  `condor_submit pe.sub`). Job logs appear only when a job ends.
- To update the code there: `rsync -a --exclude __pycache__ raven_kn cit:raven_kn_em_priors/`.
