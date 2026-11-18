#!/bin/bash
# From cron every 5 minutes, as ensure_listener.sh:
#   */5 * * * * $HOME/raven-kn-run/ensure_serve.sh
pgrep -u "$(id -u)" -f "raven-kn-run/run_serve\.sh" > /dev/null && exit 0
echo "$(date -Is) run_serve.sh not running: started by cron on $(hostname)" >> "$HOME/raven-kn-run/serve.log"
/usr/bin/screen -dmS raven-kn-serve "$HOME/raven-kn-run/run_serve.sh"
