# Logbook and its application

The store keeps the latest state only: one alert per superevent, one decision per optical transient, alerts deleted
on retraction and after 14 days. The **logbook** keeps what happened, and a read-only service gives it to
the RAVEN-KN logbook application (Flutter, Android and Linux), which copies it and shows it offline.

## What is written

| What | Where in the store | Written by |
|---|---|---|
| Journal: one JSON line per step, with the values of that moment | `journal/<source>.jsonl` | `listen`, `fink-listen`, `add`, `assess --store` |
| Heartbeat of a listener, every 30 s | `status/<source>.json` | `listen`, `fink-listen` |
| State and journal of a search or PE run | `<label>_state.json`, `<label>_journal.jsonl` in its `--out` directory | `targeted`, `runlog` |
| Reduced sky map of every GW alert (credible level per pixel, nside 128) | `gw_versions/<superevent>/` | `listen` |

Journal entries (`kind`):

| Group | Kinds |
|---|---|
| Optical | `optical_new`, `optical_update`, `host_lookup`, `lookup_failed` |
| Alerts | `gw_alert`, `gw_retraction`, `grb_notice`, `xray_position`, `nu_notice`, `nu_retraction` |
| Reports | `circular`: a GCN Circular linked to an alert of the store, with its labels; `grb_distance`: the distance of a GRB from a redshift announced in a circular, and whether a GW search can reach it ([GCN Circulars](listener.md#gcn-circulars)) |
| Associations | `grb_attached`, `grb_detached`, `nu_attached`, `nu_detached`, `coincidence` |
| Decisions | `decision` (only when it changes), `search_request` |
| Search and PE | `run_submitted`, `run_start`, `run_stage`, `run_progress`, `search_result`, `pe_result`, `run_end`, `run_failed` |
| Pipeline | `listener_start` (with the time without heartbeat before it), `kafka_error`, `message_failed` |

Each writing process has its own file, so no lock is needed on a network file system. Mock superevents
are not written in the journal; they update the listener's `last_message`, which shows that the GCN
stream is alive outside observing runs. A Kafka error is counted every time and written when its text
changes, or every 10 minutes while it repeats.

Sky-map versions are never deleted: after a retraction the directory holds a `retracted` marker.

## Runs on a cluster

`raven-kn targeted` writes `<label>_state.json` (stage, sub-segment i of N, loudest so far, exit, error) and keeps
its heartbeat alive from a thread, so a job killed without a word is seen as stale. `--pe` adds the
checks of each summary (posterior on a prior bound, Bayes factor above what the SNR allows,
maximum-likelihood SNR far from the search value) and a thinned samples file for plots.

```bash
raven-kn runlog out/ --request REQ.json --label search_free --note "condor cluster 1234"   # at submission
raven-kn targeted --request REQ.json --label search_free --out out/ ...                    # the job itself
```

The run is named by `--label`, so several jobs can share one output directory. A run is `public` when
it reads public data, else `lvk` (`--visibility` overrides); data read from local files (`--files`)
count as `lvk` unless `--visibility public` is given.
`examples/logbook_service/sync_runs.sh` copies the directories listed in `runs.txt` to
`<store>/runs/<optical>/<name>/`, with the queue state from `condor_q` (used for a directory that
holds a single run). It has not been tested on CIT.

The SNR time series of the loudest template is still written by `examples/gw170817_reverse_snr.py`
only; the service gives `<label>_snr_timeseries.json` when it is in the run directory.

## The service

```bash
export RAVEN_KN_TOKEN=...          # shows public runs
export RAVEN_KN_LVK_TOKEN=...      # optional: also shows runs on LVK data
raven-kn serve --store kn_store/ --host 0.0.0.0 --port 8765
```

Standard library only, read-only. Every request needs `Authorization: Bearer <token>`; a refusal is a
JSON 401 with `"service": "raven-kn"`, and every answer carries the header `X-Raven-KN`, so that the
application can tell the service from a portal login page answering in its place. Entries and runs
marked `lvk` are given to the LVK token only. The service speaks plain HTTP: keep it inside the Virgo
network or put a TLS proxy in front of it.

| Request | Gives |
|---|---|
| `GET /api/health` | listeners (alive if the heartbeat is younger than 3 min), store counts, age of the run copy |
| `GET /api/events?cursor=&limit=` | journal entries after the cursor, oldest first; `more` when there are others |
| `GET /api/optical` | optical transients with their steps |
| `GET /api/optical/<name>` | optical transient, decision, coincidences, runs, steps; `reported_photometry`: the magnitudes reported in GCN Circulars about it or its burst, not verified |
| `GET /api/runs/<optical>/<run>` | state, results, SNR series, PE summaries and samples |
| `GET /api/gw` | superevents and their sky-map versions |
| `GET /api/gw/<id>/map/<n>/<gw\|combined>` | reduced map: gzip of one byte per pixel, nside 128 NESTED, byte = ceil(credible level × 250), 255 where empty |
| `GET /api/followups` | alerts grouped by the links the pipeline made, each named by its detection type (optical, GW, GRB, neutrino); whether each group is a candidate, and why; `notes` and `too_far` when a redshift puts its GRB beyond the range of a GW search |
| `GET /api/sky` | positions and error radii of optical transients, GRBs, X-ray positions, neutrinos; for each optical transient, the GRB, neutrino and superevents attached to it; for each superevent, the stored GRB of LVK's own GW × GRB coincidence, matched by notice number, else by mission and time; for each alert, `circulars`: the digest of the GCN Circulars linked to it (number, bands, optical result, first optical report, redshift), or null |
| `GET /api/circulars` | GCN Circulars linked to the alerts of the store, with their labels and links, without their text |
| `GET /api/circulars/<number>` | one of them, with its text |

The cursor is the number of lines already read in each journal file. It stays exact when an older entry
arrives late, as with a run journal copied from a cluster hours after it was written.

Steps of an optical transient, each `pending`, `running`, `done`, `held`, `warning`, `failed` or `stale`:
optical transient received, host and redshift, GRB / neutrino, search decision, GW alert coincidence, search
request, targeted search, parameter estimation.

`examples/logbook_service/` holds `run_serve.sh` and `ensure_serve.sh`, to run the service as the
listener is run.

## Test store

```bash
python examples/logbook_fixture/make_fixture.py --out /tmp/fixture_store
RAVEN_KN_TOKEN=test raven-kn serve --store /tmp/fixture_store
```

The store holds the GW170817 reverse replay (real: notices, photometry, decisions, search results),
the GW170817 alert reconstructed from its two public sky maps (in 2017 it was a circular, not a Kafka
message; `fetch_gw170817_maps.py` gives the source), and
a synthetic case, marked as such, for what the real one does not show: a neutrino, two superevents with
two map versions and a retraction, a failed run, a run on LVK data with PE summaries that fail the
checks, and a listener outage.

## Follow-ups and candidates

An alert is named by how it was detected: **optical** (a transient from Fink or added by hand; the
`OpticalTransient` of the code), **GW**, **GRB**, **neutrino**. A **follow-up** is every alert the pipeline
linked together: a GRB or a neutrino attached to an optical transient, a superevent in the merger-time
window of an optical transient, and the GW × GRB coincidence carried by an LVK alert.

A follow-up is a **candidate** when one of these holds:

- two detection types are linked;
- its optical transient is confirmed as a kilonova by the decision;
- its GW alert, not retracted, is likely to involve a neutron star: P(BNS) + P(NSBH) ≥ 0.5, or
  HasRemnant ≥ 0.5.

The pipeline does not link a neutrino to a GW alert directly, only through an optical transient. The
steps (decision, search, PE) exist only for a follow-up that has an optical transient.

## Names before November 2026

The optical transient was called a candidate in the code. Old data are still read; new data use the new
names.

| Before | Now |
|---|---|
| class `KNCandidate`, `CandidateStore`, `add_candidate`, `load_candidate`, `save_candidate` | `OpticalTransient`, `AlertStore`, `add_optical`, `load_optical`, `save_optical` (the old names remain as aliases) |
| store folder `candidates/` | `optical/` (files of an old `candidates/` are read until a newer copy exists) |
| key `candidate` in decisions, search requests, run states and journal entries | `optical` (the old key is read) |
| journal kinds `candidate_new`, `candidate_update` | `optical_new`, `optical_update` (old entries are given with the new kinds) |
| `GET /api/candidates` | `GET /api/optical` |
