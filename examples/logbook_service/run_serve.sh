#!/bin/bash
# Read-only service for the logbook application, with a restart loop. Run it in
# a screen session (ensure_serve.sh does it), beside the listener.
#   ~/.config/raven_kn_serve.env : RAVEN_KN_TOKEN=...  (public runs only)
#                                  RAVEN_KN_LVK_TOKEN=...  (also the runs on LVK data; optional)
#                                  RAVEN_KN_HOST=0.0.0.0  RAVEN_KN_PORT=8765  (mode 600)
set -a; . ~/.config/raven_kn_serve.env; set +a
cd ~/raven-kn-run
while true; do
    echo "$(date -Is) starting serve on $(hostname)" >> serve.log
    ~/raven-kn/.venv/bin/raven-kn serve --store kn_store/ \
        --host "${RAVEN_KN_HOST:-127.0.0.1}" --port "${RAVEN_KN_PORT:-8765}" >> serve.log 2>&1
    echo "$(date -Is) serve exited ($?), restarting in 60 s" >> serve.log
    sleep 60
done
