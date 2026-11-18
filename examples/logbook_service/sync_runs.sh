#!/bin/bash
# Copy the state, journal and results of the search / PE runs from the cluster
# into the store, for the logbook. From cron on the listener host, e.g. every
# 10 minutes:
#   */10 * * * * $HOME/raven-kn-run/sync_runs.sh
#
# runs.txt, one directory per line:   <candidate> <name> <directory on the cluster>
#   SSS17a  em_priors  raven_kn_em_priors
# The directory is the --out of `raven-kn targeted`; it may hold several runs,
# one per --label. Mark a run when it is submitted with:
#   raven-kn runlog <out> --request <request.json> --label <label>
#
# The copy needs a working ssh to the cluster (ssh -fN cit when the connection
# has dropped). When it fails, .last_sync is not touched, and the application
# shows how old the last copy is. NOT TESTED on CIT (written 2026-11-18).
REMOTE="${RAVEN_KN_CLUSTER:-cit}"
RUN="$HOME/raven-kn-run"
STORE="$RUN/kn_store"
ok=1
while read -r cand run dir; do
    case "$cand" in ""|\#*) continue;; esac
    dest="$STORE/runs/$cand/$run"
    mkdir -p "$dest"
    rsync -a --timeout=60 -e "ssh -o BatchMode=yes -o ConnectTimeout=20" \
        --include='*_state.json' --include='*_journal.jsonl' --include='*snr_timeseries.json' \
        --include='*_targeted_result.json' --include='*_pe_summary.json' --include='*_pe_samples.json' \
        --exclude='*' "$REMOTE:$dir/" "$dest/" 2>> "$RUN/sync.log" || ok=0
done < "$RUN/runs.txt"
# queue state (idle, running, held with its reason), matched to runs by directory name
ssh -o BatchMode=yes -o ConnectTimeout=20 "$REMOTE" \
    'condor_q -json -attributes ClusterId,ProcId,JobStatus,HoldReason,Iwd' \
    > "$STORE/runs/condor.json.tmp" 2>> "$RUN/sync.log" \
    && { [ -s "$STORE/runs/condor.json.tmp" ] || echo "[]" > "$STORE/runs/condor.json.tmp"; \
         mv "$STORE/runs/condor.json.tmp" "$STORE/runs/condor.json"; } || ok=0
if [ "$ok" = 1 ]; then touch "$STORE/runs/.last_sync"; else echo "$(date -Is) sync failed" >> "$RUN/sync.log"; fi
