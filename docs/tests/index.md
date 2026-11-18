# Strategy tests

A strategy is the set of rules that takes an optical transient to a GW search: which evidence is asked
before searching, how the merger-time window is set, which statistic and which electromagnetic priors the
search uses. A strategy test replays a known event with one strategy and records what it would have done:
when the search is triggered, whether the window contains the merger, how the signal stands against the
background.

These replays are kept apart from the live stream. A test feeds archived alerts and photometry, in the
order they arrived, into a store of its own; nothing is read from or written to the store of the
[listener](../listener.md), and the GW alert of the event is withheld.

| Test | Inputs | Reference result |
|---|---|---|
| [GW170817 from the GRB and the kilonova](gw170817-reverse.md) | Fermi-GBM notice, real AT2017gfo light curve, GWOSC strain | search triggered 29.5 h after the GRB, GW170817 found at reweighted coherent SNR 30.3, 2 ms from the merger time, p < 0.01 |

## What a test varies

| Stage | Rule | Option of `examples/gw170817_reverse.py` | In the pipeline |
|---|---|---|---|
| Trigger | kilonova evidence required before searching | `--no-kn-evidence` | `raven-kn assess --no-kn-evidence` |
| Trigger | separation of the two epochs of the light-curve test (0.25 d) | `--min-epoch-sep` | `kilonova_evidence(min_epoch_sep_d=)` |
| Trigger | minimum fading rate (0.3 mag/day) and its significance (3σ) | `--min-decline`, `--min-sigma` | `kilonova_evidence(min_decline=, min_sigma=)` |
| Search | detection statistic | `--statistic coherent\|network` | `raven-kn targeted --statistic` |
| Search | chirp-mass range | `--mchirp` | `raven-kn targeted --mchirp-min --mchirp-max` |
| Search | polarisation and distance prior from the EM data | `examples/em_priors_gw170817/` | `raven-kn targeted --polarisation --distance-prior` |
| Follow-up | parameter estimation with or without the EM priors | `examples/em_priors_gw170817/` | `raven-kn targeted --pe em\|gw\|both` |

The trigger rules are tried in seconds (`--skip-search`). The searches need the strain and, for a
background, HTCondor jobs.

## Limits

One event is one test: GW170817 is close, loud and attached to a GRB, so it says little about a faint
kilonova without a GRB, where the window is hours long. The synthetic kilonovae of the
[injection calibration](../examples/injections.md) cover the window; no replay of a second real event
exists yet.
