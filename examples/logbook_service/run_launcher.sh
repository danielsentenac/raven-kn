#!/bin/bash
# Launcher of the targeted searches asked for from the logbook application, with a
# restart loop. Run it in a screen session (ensure_launcher.sh does it), beside
# the listener and the service, on the machine that sees the frames of the site.
#   RAVEN_KN_LIBS : library directory to put first in LD_LIBRARY_PATH, when PyCBC
#                   does not load otherwise (the conda environment's lib/)
#   RAVEN_KN_FARM : nodes to run the searches on (farmn1,farmn2,...), reached by
#                   ssh without a password, with the same home and frames;
#                   without it the searches run on this machine
cd ~/raven-kn-run
[ -f ~/.config/raven_kn_launcher.env ] && { set -a; . ~/.config/raven_kn_launcher.env; set +a; }
[ -n "$RAVEN_KN_LIBS" ] && export LD_LIBRARY_PATH="$RAVEN_KN_LIBS:$LD_LIBRARY_PATH"
while true; do
    echo "$(date -Is) starting launcher on $(hostname)" >> launcher.log
    ~/raven-kn/.venv/bin/raven-kn launcher --store kn_store/ >> launcher.log 2>&1
    echo "$(date -Is) launcher exited ($?), restarting in 60 s" >> launcher.log
    sleep 60
done
