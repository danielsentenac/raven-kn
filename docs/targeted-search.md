# Targeted GW search

Without a GW alert, the window, sky position and distance define a targeted search of the archived strain
data. Knowing the sky position and a narrow time window lowers the detection threshold well below that of
an all-sky search. RAVEN-KN provides two routes.

## Light search (`raven-kn targeted`)

A PyCBC-based search that runs on a laptop or a single cluster node:

- **Template bank**: equal-mass, non-spinning TaylorF2 templates over a chirp-mass range set by the
  distance, spaced for a minimal match of 0.97 (181 templates for 1.18–1.22 M☉, about 2150 for 1.0–1.5,
  about 10 000 for sub-solar 0.4–0.9).
- **Data**: public GWOSC data or LVK frames ([below](#data-public-or-lvk-raven_kntargeteddatasource)), or
  local HDF5 or frame files, in sub-segments of 1024 s by default. Loud glitches are gated
  (zeroed ±0.25 s with a 0.25 s taper, as PyGRB).
- **Statistic**: either
    - *network*: each detector's SNR, shifted by its time delay for the known position, combined as
      √Σρ², counted only when two detectors exceed 4; or
    - *coherent* (`--statistic coherent`, as PyGRB): complex SNRs projected onto the network's response to
      the two polarisations at the known position; the remainder (null stream) is cut as in PyGRB.
- **Signal consistency**: a power chi-square where a detector's SNR exceeds 5.
    - *network*: 16 bins; each detector's SNR is replaced by its newSNR.
    - *coherent*: as PyGRB, one chi-square summed over the detectors reweights the coherent SNR, with
      PyGRB's number of bins (about 125 for a binary neutron star). `coherent_chisq="per_ifo"` reweights
      each detector before the projection instead.
    - *Gate-aware*: where a gate falls inside the template, the chi-square uses the template gated like
      the data, so a real signal partly removed by a gate is not taken for a glitch.
    - The chi-square is computed at the 500 loudest samples per template and detector; other samples over
      the threshold take the value of the nearest computed sample within 1 s, so the wings of a loud
      glitch are down-weighted too.
- **Background**: time slides with independent shifts per detector; the loudest on-source candidate is
  removed from the background only if loud in two detectors. Optionally, off-source trials
  (`--n-offsource N`) sample the data around the window, as PyGRB does.
- **Significance**: a p-value per sub-segment, combined over sub-segments with a trials factor.

```bash
raven-kn targeted --request out/<name>_targeted_search_request.json --out targeted_out
```

!!! note "Signals next to a gate"
    GW170817's L1 glitch is gated 1.1 s before the merger, with part of the chirp. The plain per-detector
    chi-square then lowered L1 from 24 to 16 (coherent 24.4). With the gate-aware network chi-square the
    reweighted coherent SNR is 29.8, as PyGRB's (see [the comparison](examples/gw170817-search.md#comparison-with-pygrb)).

## Data: public or LVK (`raven_kn.targeted.datasource`)

The strain can come from two places. Both give the search the same two things, the analysable time of each
detector and its strain, so the search does not depend on which one is used.

| | Public (GWOSC) | LVK |
|---|---|---|
| Credentials | none | IGWN SciToken (`htgettoken -a vault.ligo.org -i igwn`) |
| Available | catalogue events; whole runs after the proprietary period | at once |
| Segments | GWOSC timelines (`H1_CBC_CAT1`, `H1_CBC_CAT2`); no DQSegDB | DQSegDB: science flag, veto definer |
| Frames | `datafind.gwosc.org` (`gw_data_find -r datafind.gwosc.org --show-types`): CVMFS where mounted, else OSDF, else gwosc.org's archive | IGWN datafind; on a cluster, files read in place |
| Conditioning | already analysis-ready | calibrated h(t) frames: science segments, CAT1 removed, gating |

- **Analysable time** is science data passing CAT1. A detector that is not analysable over a sub-segment is
  left out of it. CAT2 time is kept, and the result says if the loudest candidate falls in it.
- **Routing** (`--source auto`, the default): public data if they cover the window in at least two
  detectors; else LVK data if configured and a SciToken is present; else the search stops and says that the
  data are not public yet. `--source public` or `--source lvk` forces one.
- **LVK frames**: the calibrated h(t) frames, never raw data. Frame types, channels and science flags are
  given per detector, or by a preset checked on real data (`--lvk-preset O2`: cleaned C02 for H1 and L1,
  `V1O2Repro2A`).
- **Provenance**: each result file records the source, frame types, channels, frame files and segments
  used (`data`). The same event gives different numbers on different frames.

```bash
# public data, no credentials
raven-kn targeted --request REQ.json --source public
# LVK data on a cluster with a local datafind server
raven-kn targeted --request REQ.json --source lvk --lvk-preset O2 --veto-definer VDF.xml
```

Checked live on GW170817: the public back-end from a laptop (gwosc.org archive, the OSDF caches being
unreachable there) and on CIT (OSDF), the LVK back-end on CIT (DQSegDB, frames from `/ceph`). For production
on LVK data the engine is PyGRB, whose workflow does the same conditioning
([GW170817 on LVK frames](examples/gw170817-search.md#pygrb-on-lvk-frames)).

## Production: PyGRB (`raven-kn pygrb`)

PyGRB, LVK's coherent targeted search for GRBs, is built for on-source windows of seconds to minutes. RAVEN-KN
splits a kilonova window into sub-segments and writes, for each, an event file and a
`pycbc_make_offline_grb_workflow` command, sharing one template bank. Base configurations (e.g. LVK's O2
GRB170817A test configuration) are reused; event options are replaced and durations passed as overrides.

PyGRB aborts if its on-source window contains vetoed time, so sub-segments are placed in veto-free stretches.
The vetoes can be given as a file or queried from DQSegDB:

```bash
raven-kn pygrb --request out/<name>_targeted_search_request.json --out pygrb_configs \
    --base-config analysis.ini executables_common.ini executables_condorpool.ini injections.ini \
                  postprocessing.ini cats.ini o2-data.ini \
    --event-template GRB170817A.ini --bank-path /path/on/cluster/<name>_bank.hdf \
    --veto-definer VDF.xml --veto-categories 1 2 --science-flag H1=H1:DMT-ANALYSIS_READY:1 ...
```

DQSegDB is at `https://segments.igwn.org` and needs an IGWN SciToken
(`htgettoken -a vault.ligo.org -i igwn`).

On LIGO's CIT cluster, HTCondor jobs submitted from the login nodes need `+flock_local = True` and
`MY.POOLS = "IGWN,CIT"` to run on local nodes. Jobs start without `HOME`, so set `HOME` and `NUMBA_CACHE_DIR`
to the scratch directory.

See [the GW170817 example](examples/gw170817-search.md) for both routes on real data.
