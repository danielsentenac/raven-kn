#!/bin/bash
# RAVEN-KN GCN listener with a restart loop. Run it inside a screen session
# (ensure_listener.sh does it); paths assume the code in ~/raven-kn and the
# run directory in ~/raven-kn-run.
#   ~/.config/gcn_kafka.env : GCN_CLIENT_ID=... and GCN_CLIENT_SECRET=... (mode 600)
#   kn_store/               : candidates, alerts, notices, coincidences
#   results.jsonl           : one line per coincidence result
#   listener.log            : starts, exits, Kafka errors, candidates updated
set -a; . ~/.config/gcn_kafka.env; set +a
cd ~/raven-kn-run
while true; do
    echo "$(date -Is) starting listener on $(hostname)" >> listener.log
    ~/raven-kn/.venv/bin/raven-kn listen --store kn_store/ >> results.jsonl 2>> listener.log
    rc=$?
    echo "$(date -Is) listener exited ($rc), restarting in 60 s" >> listener.log
    sleep 60
done
