#!/bin/bash
# Run from cron every 5 minutes: start the GCN listener if its restart loop
# (run_listener.sh, in the screen session "gcn-listener") is not running,
# e.g. after a reboot or if the screen session was killed. Crontab line:
#   */5 * * * * $HOME/raven-kn-run/ensure_listener.sh
pgrep -u "$(id -u)" -f "raven-kn-run/run_listener\.sh" > /dev/null && exit 0
echo "$(date -Is) run_listener.sh not running: started by cron on $(hostname)" >> "$HOME/raven-kn-run/listener.log"
/usr/bin/screen -dmS gcn-listener "$HOME/raven-kn-run/run_listener.sh"
