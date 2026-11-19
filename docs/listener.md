# Real-time listener

`raven-kn listen` subscribes to GCN Kafka and matches every incoming alert against a store of kilonova
optical transients, in both time orders: an alert after an optical transient, or an optical transient after an alert.

```bash
export GCN_CLIENT_ID=... GCN_CLIENT_SECRET=...   # from gcn.nasa.gov; pip install gcn-kafka
raven-kn listen --store kn_store/                 # GW alerts, GRB notices, neutrino alerts, circulars
raven-kn add optical.json --store kn_store/     # add or update an optical transient
```

## What it does with each message

| Message | Topic | Action |
|---|---|---|
| GW alert | `igwn.gwalert` | Decodes the multi-order sky map in the alert. For every stored optical transient whose window contains the merger time, computes the joint FAR and writes `kn_store/coinc/<superevent>__<optical>.json`. Updates replace the stored alert; retractions delete it. |
| GRB notice, unified JSON | Swift-BAT GUANO, Einstein Probe WXT, SVOM-ECLAIRs | Stored, and attached to optical transients within 2× the error radius whose first detection follows within 3 days; their window becomes [−5, +1] s around the GRB. |
| GRB notice, VOEvent | Fermi-GBM flight / ground / final positions, Swift-BAT positions | Same, with the 3.7° GBM systematic error added. |
| X-ray afterglow position, VOEvent | Swift-XRT (`gcn.classic.voevent.SWIFT_XRT_POSITION`), SVOM-MXT (`gcn.notices.svom.voevent.mxt`) | Replaces the position of the stored Swift-BAT or SVOM-ECLAIRs burst with the same trigger number; the burst time is kept. Optical transients are matched again within 2× the X-ray radius; one that the listener had attached to this burst and that lies outside is detached. |
| New optical transient | `raven-kn add` | Matched against stored GRB notices and neutrino alerts (the latest one up to 3 days before its first detection) and against GW alerts of the last 14 days. |
| Neutrino alert | `gcn.notices.icecube.gold_bronze_track_alerts` | Stored, and attached to optical transients within 1.4× its 90 % error radius whose first detection follows within 3 days; their window becomes ±500 s around the neutrino, unless a GRB gives a tighter one. Updates refine the time and position; retractions remove the alert and detach it. |
| GCN Circular | `gcn.circulars` | Linked to the alerts of the store it is about, labelled (bands, optical detection or upper limits, redshift) and stored: see [GCN Circulars](#gcn-circulars). `--no-circulars` leaves them out. |

## The store

Everything the listener knows is in the `--store` directory, as JSON files that can be read and edited by
hand. `raven-kn listen`, `raven-kn add` and `raven-kn fink-listen` share it.

| Subdirectory | Contents | Written when |
|---|---|---|
| `optical/` | one file per optical transient: position, photometry, host, attached GRB and neutrino, notes | `raven-kn add`, `raven-kn fink-listen`; updated when a notice is attached or detached |
| `gw/` | latest alert per superevent, with its multi-order sky map (and RAVEN's combined map if any) | each GW alert; deleted on retraction; used for 14 days |
| `grb/` | latest notice per trigger: time, position, radius, and the X-ray position when there is one | each GRB or X-ray notice |
| `nu/` | latest alert per IceCube event | each neutrino alert; deleted on retraction |
| `coinc/` | `<superevent>__<optical>.json`: window, joint FAR and its factors | each GW alert or new optical transient whose window contains the merger time |
| `circulars/` | one file per GCN Circular linked to an alert: subject, text, labels, links; `unlinked/` holds for 30 days those that name an event not linked yet | each circular |

Each coincidence result is also printed on standard output as one JSON line. The log (standard error)
records the topics subscribed to, the optical transients updated by each notice, Kafka errors, and the traceback
of any message that could not be handled; such a message is skipped and the listener goes on.

## Running it permanently

The listener is one long-lived process. `examples/listener_service/` holds the two scripts used for the
production instance:

| Script | Role |
|---|---|
| `run_listener.sh` | Reads the GCN credentials from `~/.config/gcn_kafka.env`, runs `raven-kn listen` in a loop and restarts it 60 s after any exit. Results go to `results.jsonl`, the log to `listener.log`. |
| `ensure_listener.sh` | Starts `run_listener.sh` in a `screen` session if it is not running, and does nothing otherwise. |

```bash
mkdir ~/raven-kn-run && cp examples/listener_service/*.sh ~/raven-kn-run/
crontab -e      # add:  */5 * * * * $HOME/raven-kn-run/ensure_listener.sh
```

The cron line brings the listener back within 5 minutes of a reboot or of a kill (tested by stopping it:
cron restarted it at the next tick). An `@reboot` entry is not enough: it runs once, and failed on two
consecutive boots of the production server, where the home directory is a network mount. With the fixed
consumer group a restarted listener resumes from its last committed offset, so the messages of the gap
are processed as long as the broker still retains them.

To update the code (installed with `pip install -e`): `git pull` in the clone, then
`pkill -f "raven-kn listen"`; the loop starts the new version 60 s later (tested).

## GCN Circulars

A notice gives the time and position of an alert. Whether anybody searched for its optical counterpart,
and what was found, is reported by the observers in GCN Circulars, which are free text. The listener
reads them (`raven_kn.circulars`) and does three things with each one.

**Links** to the alerts of the store, in this order:

| Alert | Linked when |
|---|---|
| GW | a superevent of the store is named in the subject |
| Neutrino | a stored IceCube event is named in the subject |
| Optical transient | its name is in the subject or the text |
| GRB or X-ray notice | its trigger number is in the subject or the text (Swift `trigger=1263718`, Fermi `trigger 751960088.32731 / 241030242`, SVOM `sb26111302`) |
| | else a time written in the text is within 1.5 s of its trigger time, on the date of the event name (`GRB 261113A`: 2026-11-13). Einstein Probe notices carry no trigger number |
| | else the event name, once a circular linked by one of the two rules above has carried it: most follow-up reports only say `GRB 261113A` |

A circular linked to nothing is kept 30 days in `circulars/unlinked/` when it names an event, and is
linked when a later circular gives the trigger number of that event.

**Labels**, from keywords of the subject: the bands (`optical`, `xray`, `radio`, `gamma`, `gw`, `nu`),
and for an optical report its result: `detection`, `upper_limit`, or `observation` when neither the
subject nor the text says which. A redshift announced in the subject or the text (`z = 3.580`) is
recorded with the circular; it is not applied to an optical transient. These are keyword rules, not a
reading of the text: the subject is kept beside the labels, and the list of optical facilities
(`BANDS` in `raven_kn/circulars.py`) is the one met in the circulars of October and November 2026.

**Journal and service**: a `circular` entry is written when a circular is linked, and
[`/api/sky` and `/api/followups`](logbook.md) give for each alert the digest of its circulars: how
many, in which bands, whether an optical search was reported and its best outcome, the time of the
first optical report, the latest redshift.

```bash
raven-kn circular 46261 46270-46290 --store kn_store/    # from the archive, by number or range
```

reads circulars from the archive (`https://gcn.nasa.gov/circulars/<number>.json`): those sent before
the listener subscribed, or all of them again after a change of the rules. A circular already stored is
replaced, and written in the journal again only if its links changed.

Checked on 2026-11-19 with the 441 circulars 45880–46320 (10 October to 18 November 2026) and the 17
notices of the production store (10 to 18 November): 72 circulars were linked, to 11 of the 17 notices.

| Notices | Circulars | Optical search reported |
|---|---|---|
| SVOM-ECLAIRs sb26111302 (GRB 261113A) | 34 | detection, 27 reports, z = 3.58 |
| SVOM-ECLAIRs sb26111401 (GRB 261114A) | 10 | detection, 7 reports, z = 3.716 |
| Einstein Probe WXT (EP261115a) | 6 | upper limits, 3 reports |
| Einstein Probe WXT (EP261114a) | 3 | upper limits, 1 report |
| Fermi-GBM 816291960 (GRB 261113B) | 6 | 1 report: an optical fading source that may be associated |
| Fermi-GBM 816089232, 816473458, 816645341 | 4, 3, 1 | upper limits, 1 report each (MASTER) |
| Fermi-GBM 816659441, 816700674 | 1 each | none (the Fermi GBM localisation only) |
| Fermi-GBM, 6 other triggers | 0 | none |

Of the 72 circulars, 22 were linked by a trigger number, 8 by a trigger time (among them the first
circular of each of the two Einstein Probe transients) and 42 by the event name alone.

## Optical transients from Fink

`raven-kn fink-listen` fills the same store from the Fink broker's Kafka livestream (ZTF or Rubin):

```bash
pip install .[fink]                      # fink-client
finkctl auth register -survey ztf -username ... -groupid ... -servers ...   # values sent by Fink
raven-kn fink-listen --survey ztf --store kn_store/ --find-host --assess
```

Each alert becomes or updates an optical transient and goes through the same matching as `raven-kn add`
(stored GRB notices, neutrino alerts, GW alerts). `--find-host` looks the host galaxy up once per
optical transient, which the light-curve kilonova test needs; `--assess` re-runs the search decision on every
update and logs when a search becomes worth it. For ZTF the default topics are Fink's three kilonova
streams; Fink has no Rubin kilonova stream yet, so Rubin topics are chosen with `--topics`
(`finkctl topic list -survey lsst`). The closest Rubin topics are
`fink_fast_transient_gvom_lsst` (new, brighter than 21, extragalactic, fast-evolving) and
`fink_extragalactic_new_candidate_lsst` (new, brighter than 24, rising or fading). Both surveys were
checked on live alerts in November 2026.

If an LVK alert carries its own GW + GRB coincidence (`external_coinc`), it is kept and copied into every
RAVEN-KN result, with the combined sky map at the optical transient.

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

Outbound port 9092 is blocked on some campus networks, so the listener may have to run elsewhere: on a
server as described [above](#running-it-permanently), or on a cluster as an HTCondor local-universe job.
