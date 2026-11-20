#!/bin/bash
# From cron every 5 minutes, as ensure_listener.sh:
#   */5 * * * * $HOME/raven-kn-run/ensure_launcher.sh
pgrep -u "$(id -u)" -f "raven-kn-run/run_launcher\.sh" > /dev/null && exit 0
echo "$(date -Is) run_launcher.sh not running: started by cron on $(hostname)" >> "$HOME/raven-kn-run/launcher.log"
/usr/bin/screen -dmS raven-kn-launcher "$HOME/raven-kn-run/run_launcher.sh"
