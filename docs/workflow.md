# Workflow overview

![RAVEN-KN workflow](raven-kn-workflow.svg)

RAVEN-KN works in the opposite direction to most multi-messenger pipelines: it starts from the
**electromagnetic counterpart** and asks what it implies for the GW data.

```
kilonova candidate → is a GW search worth it? → merger-time window, sky position, distance, inclination
                   → coincidence with existing GW alerts (joint FAR) and/or targeted search of archived strain
```

## 1. Inputs

Four kinds of alerts arrive as live streams; the rest is looked up when a candidate exists.

| Alert stream | Gives | Source and topics | Command |
|---|---|---|---|
| Kilonova candidates, ZTF | new candidate or new photometry, Fink kilonova score | Fink Kafka livestream: `fink_kn_candidates_ztf`, `fink_early_kn_candidates_ztf`, `fink_rate_based_kn_candidates_ztf` | `raven-kn fink-listen --survey ztf` |
| Transient candidates, Rubin | new candidate or new photometry (no kilonova score) | Fink Kafka livestream; no kilonova topic yet, closest are `fink_fast_transient_gvom_lsst` and `fink_extragalactic_new_candidate_lsst` | `raven-kn fink-listen --survey lsst --topics ...` |
| GRB and X-ray notices | prompt-emission time and position | GCN Kafka: Fermi-GBM (flight, ground, final positions), Swift-BAT (GUANO, on-board positions), SVOM-ECLAIRs, Einstein Probe WXT; Swift-XRT and SVOM-MXT afterglow positions, which replace the position of their burst | `raven-kn listen` |
| Neutrino alerts | arrival time and position | GCN Kafka: IceCube gold / bronze track alerts | `raven-kn listen` |
| GW alerts | merger time, sky map, distance, FAR | GCN Kafka: LVK `igwn.gwalert` (preliminary, initial, update, retraction) | `raven-kn listen` |

All streams fill the same store, and every new message is matched against what is already there,
whichever arrives first. Both Fink streams and the GCN topics were checked on live alerts in November 2026.

| Looked up per candidate | Gives | Sources |
|---|---|---|
| Optical light curves | detections, and non-detections with their depth | ZTF, Rubin and ATLAS through ALeRCE, ANTARES, Fink, and ZTF and ATLAS forced photometry |
| Spectra and host galaxy | classification, redshift, distance | TNS, NED / NED-D, GLADE+, DES Y6 |
| Viewing angle | inclination prior | external: afterglow modelling, VLBI, light-curve fits |

A candidate can also be added by hand (`raven-kn add`).

GRB notices usually arrive hours before the optical counterpart is found (GW170817: GRB 170817A at
12:41 UTC, AT2017gfo at 23:33 UTC). The listener therefore **stores every GRB notice** and attaches it to a
kilonova candidate found up to 3 days later within twice the notice's error radius. IceCube neutrino
alerts are handled the same way, within 1.4 times their 90 % error radius; an attached neutrino sets the
merger-time window to ±500 s around it (the standard GW–neutrino coincidence window), unless a GRB gives a
tighter one.

## 2. Connection: the kilonova candidate

All inputs are merged into one `KNCandidate` (`raven_kn.connection`): position, photometry including
non-detections, redshift and distance, viewing angle, GRB and neutrino associations, classification and classifier score,
and the samples of a light-curve fit. Photometry from several surveys is merged without duplicates, so
the latest deep non-detection can come from any of them.

The candidate can be exported as a GraceDB external-event dictionary that LVK RAVEN accepts directly.

## 3. Treatment

Three quantities are derived from the candidate. They are explained in detail on their own pages.

- **[Merger-time window t0](t0-window.md)**: the stretch of time in which the merger happened.
  Everything downstream depends on it. A GRB gives it to seconds; otherwise it comes from a kilonova
  light-curve fit, widened by a model-systematic margin, and cut by the last non-detection before the
  first detection.
- **[Distance and inclination](coincidence.md#distance-and-inclination)**: from the host redshift (best),
  a redshift-independent host distance, or the light-curve fit (least reliable); the viewing angle from
  external measurements.
- **[Kilonova rate](coincidence.md#kilonova-candidate-rate)** R_KN: how often candidates like this one
  appear by chance, the analogue of RAVEN's GRB rate.

## 4. Decision: is a GW search worth it?

`raven-kn assess` triggers a targeted search only when two conditions hold.

**① The transient is confirmed as a kilonova.** A GRB, or a single optical detection, is not enough: most
optical transients in a GRB error box are unrelated. Any one of these confirms it
(`raven_kn.treatment.kilonova_evidence`):

| Evidence | Condition |
|---|---|
| light curve | ≥ 2 detection epochs at least 6 h apart; peak absolute magnitude between −17.5 and −12.5 at the host distance; fading or reddening of at least 0.3 mag/day, significant at 3σ, within 4 days of discovery |
| spectrum | a spectroscopic classification "kilonova" (TNS, GCN circular) |
| classifier score | an upstream kilonova classifier score ≥ 0.8 |

Supernovae evolve at less than 0.1 mag/day in their first days, and fast blue optical transients peak at
−19 to −21, so the light-curve test rejects both.

**② A signal would be detectable.** The expected network SNR of a binary neutron star at the candidate's
distance and sky position, over the window (detector antenna patterns, Earth rotation, duty cycle), must
reach the threshold (default 6, probability ≥ 0.1), and the candidate must lie within 500 Mpc. A candidate
at z = 0.15 (737 Mpc) is rejected; an AT2017gfo-like one at 45 Mpc gives an expected SNR of about 50.

Until both hold, the candidate is **held** and re-assessed when new photometry or a classification arrives.
In the [GW170817 replay](examples/gw170817-reverse.md), the GRB alone and the first two nights of
photometry do not trigger; the fading on the second night does, 29.5 h after the GRB.
`--no-kn-evidence` turns the first condition off.

## 5a. With a GW alert: coincidence and joint FAR

When an LVK alert arrives whose merger time falls inside a stored candidate's window, RAVEN-KN computes a
joint FAR from the window length, the kilonova rate, the GW FAR, and the overlap of the two sky maps and
distance distributions. See [Distance, rate and joint FAR](coincidence.md).

## 5b. Without a GW alert: targeted search

The window, sky position and chirp-mass range (from the distance) become a search request. RAVEN-KN runs a
light PyCBC search itself, or writes a configuration for PyGRB, LVK's production targeted search. See
[Targeted GW search](targeted-search.md).

## 6. Real time

Two processes share one store. The [GCN listener](listener.md) (`raven-kn listen`) runs branch 5a for
every GW alert and every stored candidate, and attaches GRB notices and neutrino alerts to candidates,
whichever arrives first. The [Fink listener](listener.md#kilonova-candidates-from-fink)
(`raven-kn fink-listen`) creates and updates the candidates from the ZTF and Rubin alert streams and,
with `--assess`, re-runs the decision of step 4 on every update, so that a search is triggered as soon
as the light curve confirms a kilonova.
