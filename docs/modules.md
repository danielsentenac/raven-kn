# Modules and commands

## Modules

The package mirrors LVK RAVEN so that it can be merged into it later.

| Module | Role | RAVEN analogue |
|---|---|---|
| `raven_kn.connection` | `KNCandidate` and connectors: Fink (ZTF, Rubin), ALeRCE, ANTARES, ATLAS, ZTF forced photometry (IPAC), GraceDB, GCN JSON, NED / NED-D, GLADE+, DES Y6, TNS, local JSON | `gracedb_events.ExtTrig` / `SE` |
| `raven_kn.treatment.t0` | merger-time window: GRB, light-curve-fit posterior with margin, photometric bracket | fixed `tl, th` |
| `raven_kn.treatment.lightcurve_fit` | NMMA wrapper: candidate → NMMA files and command → t0, distance and inclination samples | — |
| `raven_kn.treatment.distance_inclination` | joint (D_L, cos ι) posterior and search priors | — |
| `raven_kn.treatment.detectability` | expected network SNR over the window; trigger decision | — |
| `raven_kn.treatment.far` | joint FARs, targeted-search p-values and trials | `search.coinc_far` |
| `raven_kn.skymap` | kilonova 2D/3D HEALPix map, GW map reader, sky and distance overlaps | `search.skymap_overlap_integral` |
| `raven_kn.listener` | GCN Kafka listener | `listen` / `handle_*` |
| `raven_kn.fink_stream` | Fink livestream (ZTF, Rubin) → candidates | — |
| `raven_kn.search` | `query`, `search`, `coinc_far`, `calc_signif_gracedb`, `assess` | `search.*` |
| `raven_kn.targeted` | light PyCBC search, public (GWOSC) and LVK data sources, PyGRB configuration, DQSegDB segments | — |

`KNCandidate.to_raven_event_dict()` gives a GraceDB external-event dictionary accepted by
`ligo.raven.gracedb_events.ExtTrig(..., event_dict=..., use_radec=True)`.

## Commands

| Command | Does |
|---|---|
| `raven-kn assess CAND.json --out DIR` | window, distance and inclination, detectability decision, search request |
| `raven-kn coinc CAND.json --superevent-far F --gw-skymap MAP` | joint FAR against one GW alert |
| `raven-kn nmma CAND.json --workdir DIR` | NMMA light-curve fit (`--dry-run`, `--import-only` for cluster runs) |
| `raven-kn host CAND.json` | host galaxy and redshift (NED, GLADE+, DES Y6; `--no-des`) |
| `raven-kn tns CAND.json --objname NAME` | TNS classification and redshift |
| `raven-kn ztffp submit\|fetch CAND.json` | ZTF forced photometry from IPAC (queued; fetch merges the light curve) |
| `raven-kn targeted --request REQ.json` | light targeted search (`--statistic network\|coherent`) |
| `raven-kn pygrb --request REQ.json` | PyGRB configuration per sub-segment |
| `raven-kn listen --store DIR` | real-time GCN listener |
| `raven-kn fink-listen --survey ztf\|lsst --store DIR` | kilonova candidates from Fink's livestream into the same store |
| `raven-kn add CAND.json --store DIR` | add a candidate to the listener's store |

## Main parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--model-sys-frac` | 0.75 | model-systematic margin on a light-curve t0 |
| `fit_days` | 3 | days of data fitted after the first detection |
| BNS ranges | H1 150, L1 160, V1 55 Mpc | approximate O4 values; pass the live ones |
| SNR threshold, P_min, distance cut | 6, 0.1, 500 Mpc | search trigger decision |
| `--chirp-mass` | BNS, 1.2188 M☉ | rescales ranges as Mc^(5/6), e.g. for sub-solar alerts |
