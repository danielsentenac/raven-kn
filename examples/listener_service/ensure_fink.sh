#!/bin/bash
# From cron every 5 minutes, as ensure_listener.sh: start the Fink listener of a
# survey (run_fink.sh, in the screen session "fink-<survey>") if it is not running.
#   */5 * * * * $HOME/raven-kn-run/ensure_fink.sh ztf
#   */5 * * * * $HOME/raven-kn-run/ensure_fink.sh lsst fink_fast_transient_gvom_lsst
survey=$1
pgrep -u "$(id -u)" -f "raven-kn-run/run_fink\.sh $survey" > /dev/null && exit 0
echo "$(date -Is) run_fink.sh $survey not running: started by cron on $(hostname)" >> "$HOME/raven-kn-run/fink_$survey.log"
/usr/bin/screen -dmS "fink-$survey" "$HOME/raven-kn-run/run_fink.sh" "$@"
