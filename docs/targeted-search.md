# Targeted GW search

Without a GW alert, the window, sky position and distance define a targeted search of the archived strain
data. Knowing the sky position and a narrow time window lowers the detection threshold well below that of
an all-sky search. RAVEN-KN provides two routes.

## Light search (`raven-kn targeted`)

A PyCBC-based search that runs on a laptop or a single cluster node:

- **Template bank**: equal-mass, non-spinning TaylorF2 templates over a chirp-mass range set by the
  distance, spaced for a minimal match of 0.97 (181 templates for 1.18–1.22 M☉, about 2150 for 1.0–1.5,
  about 10 000 for sub-solar 0.4–0.9).
- **Data**: GWOSC, local HDF5 or frame files, in sub-segments of 1024 s by default. Loud glitches are gated
  (zeroed ±0.25 s with a 0.25 s taper, as PyGRB).
- **Statistic**: either
    - *network*: each detector's SNR, shifted by its time delay for the known position, combined as
      √Σρ², counted only when two detectors exceed 4; or
    - *coherent* (`--statistic coherent`, as PyGRB): complex SNRs projected onto the network's response to
      the two polarisations at the known position; the remainder (null stream) is cut as in PyGRB.
- **Signal consistency**: a 16-bin power chi-square replaces the SNR with newSNR where a detector's SNR
  exceeds 5.
- **Background**: time slides with independent shifts per detector; the loudest on-source candidate is
  removed from the background only if loud in two detectors. Optionally, off-source trials
  (`--n-offsource N`) sample the data around the window, as PyGRB does.
- **Significance**: a p-value per sub-segment, combined over sub-segments with a trials factor.

```bash
raven-kn targeted --request out/<name>_targeted_search_request.json --out targeted_out
```

!!! note "Known limitation"
    The per-detector chi-square over-penalises a real signal partly removed by a gate. On GW170817 it lowers
    L1 from 24 to 16. A gate-aware chi-square is the next fix.

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
