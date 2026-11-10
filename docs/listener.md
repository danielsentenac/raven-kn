# Real-time listener

`raven-kn listen` subscribes to GCN Kafka and matches every incoming alert against a store of kilonova
candidates, in both time orders: an alert after a candidate, or a candidate after an alert.

```bash
export GCN_CLIENT_ID=... GCN_CLIENT_SECRET=...   # from gcn.nasa.gov; pip install gcn-kafka
raven-kn listen --store kn_store/                 # GW alerts and GRB notices
raven-kn add candidate.json --store kn_store/     # add or update a kilonova candidate
```

## What it does with each message

| Message | Topic | Action |
|---|---|---|
| GW alert | `igwn.gwalert` | Decodes the multi-order sky map in the alert. For every stored candidate whose window contains the merger time, computes the joint FAR and writes `kn_store/coinc/<superevent>__<candidate>.json`. Updates replace the stored alert; retractions delete it. |
| GRB notice, unified JSON | Swift-BAT GUANO, Einstein Probe WXT, SVOM-ECLAIRs | Attached to candidates within 2× the error radius whose first detection follows within 3 days; their window becomes [−5, +1] s around the GRB. |
| GRB notice, VOEvent | Fermi-GBM flight / ground / final positions, Swift-BAT positions | Same, with the 3.7° GBM systematic error added. |
| New candidate | `raven-kn add` | Matched against GW alerts of the last 14 days. |

If an LVK alert carries its own GW + GRB coincidence (`external_coinc`), it is kept and copied into every
RAVEN-KN result, with the combined sky map at the candidate.

Mock superevents (`MS…`, sent hourly outside observing runs) and test notices are skipped unless
`--include-mock` is given. A fixed Kafka consumer group (`--group-id`, default `raven-kn`) lets a
restarted listener resume where it stopped.

Nothing is uploaded to GraceDB from the listener.

## Live test

On 2026-11-10, 514 retained messages were replayed with no parse failure: 336 GW alerts (preliminary,
initial and retraction, one with `external_coinc`), 162 Fermi-GBM and 16 Einstein Probe notices. Real
SVOM-ECLAIRs and Swift-BAT notices parse too.

Outbound port 9092 is blocked on some campus networks, so the listener may have to run elsewhere, e.g. on a
cluster as an HTCondor local-universe job, or on a server inside a `screen` session with a restart loop.
