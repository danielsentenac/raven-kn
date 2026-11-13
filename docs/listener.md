# Real-time listener

`raven-kn listen` subscribes to GCN Kafka and matches every incoming alert against a store of kilonova
candidates, in both time orders: an alert after a candidate, or a candidate after an alert.

```bash
export GCN_CLIENT_ID=... GCN_CLIENT_SECRET=...   # from gcn.nasa.gov; pip install gcn-kafka
raven-kn listen --store kn_store/                 # GW alerts, GRB notices, neutrino alerts
raven-kn add candidate.json --store kn_store/     # add or update a kilonova candidate
```

## What it does with each message

| Message | Topic | Action |
|---|---|---|
| GW alert | `igwn.gwalert` | Decodes the multi-order sky map in the alert. For every stored candidate whose window contains the merger time, computes the joint FAR and writes `kn_store/coinc/<superevent>__<candidate>.json`. Updates replace the stored alert; retractions delete it. |
| GRB notice, unified JSON | Swift-BAT GUANO, Einstein Probe WXT, SVOM-ECLAIRs | Stored, and attached to candidates within 2× the error radius whose first detection follows within 3 days; their window becomes [−5, +1] s around the GRB. |
| GRB notice, VOEvent | Fermi-GBM flight / ground / final positions, Swift-BAT positions | Same, with the 3.7° GBM systematic error added. |
| X-ray afterglow position, VOEvent | Swift-XRT (`gcn.classic.voevent.SWIFT_XRT_POSITION`), SVOM-MXT (`gcn.notices.svom.voevent.mxt`) | Replaces the position of the stored Swift-BAT or SVOM-ECLAIRs burst with the same trigger number; the burst time is kept. Candidates are matched again within 2× the X-ray radius; one that the listener had attached to this burst and that lies outside is detached. |
| New candidate | `raven-kn add` | Matched against stored GRB notices and neutrino alerts (the latest one up to 3 days before its first detection) and against GW alerts of the last 14 days. |
| Neutrino alert | `gcn.notices.icecube.gold_bronze_track_alerts` | Stored, and attached to candidates within 1.4× its 90 % error radius whose first detection follows within 3 days; their window becomes ±500 s around the neutrino, unless a GRB gives a tighter one. Updates refine the time and position; retractions remove the alert and detach it. |

## Kilonova candidates from Fink

`raven-kn fink-listen` fills the same store from the Fink broker's Kafka livestream (ZTF or Rubin):

```bash
pip install .[fink]                      # fink-client
finkctl auth register -survey ztf -username ... -groupid ... -servers ...   # values sent by Fink
raven-kn fink-listen --survey ztf --store kn_store/ --find-host --assess
```

Each alert becomes or updates a candidate and goes through the same matching as `raven-kn add`
(stored GRB notices, neutrino alerts, GW alerts). `--find-host` looks the host galaxy up once per
candidate, which the light-curve kilonova test needs; `--assess` re-runs the search decision on every
update and logs when a search becomes worth it. For ZTF the default topics are Fink's three kilonova
streams; Fink has no Rubin kilonova stream yet, so Rubin topics are chosen with `--topics`
(`finkctl topic list -survey lsst`). The closest Rubin topics are
`fink_fast_transient_gvom_lsst` (new, brighter than 21, extragalactic, fast-evolving) and
`fink_extragalactic_new_candidate_lsst` (new, brighter than 24, rising or fading). Both surveys were
checked on live alerts in November 2026.

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

Swift-XRT and SVOM-MXT were added on 2026-11-13 from real notices (three XRT, one MXT). Their time is
that of the X-ray image, 81 s (XRT) and 131 s (MXT) after the burst in these notices, so it is never used
as a merger time: an X-ray position whose burst notice was not received is stored and attached to
nothing. Position radii are 3.6–4.7″ for XRT and 76″ for the MXT notice, against 3′ for BAT. The on-board
XRT notice of trigger 1519536 has the `Cosmic_Ray` flag set although the ground analysis confirmed its
position 72 min later, 11″ away; the flag is therefore not used.

Outbound port 9092 is blocked on some campus networks, so the listener may have to run elsewhere, e.g. on a
cluster as an HTCondor local-universe job, or on a server inside a `screen` session with a restart loop.
