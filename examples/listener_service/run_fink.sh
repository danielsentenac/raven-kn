#!/bin/bash
# RAVEN-KN listener of a Fink livestream, with a restart loop:
#   run_fink.sh ztf                              the three ZTF kilonova topics
#   run_fink.sh lsst fink_fast_transient_gvom_lsst [more topics]
# Run it inside a screen session (ensure_fink.sh does it). Each alert creates or
# updates an optical transient of the store shared with the GCN listener, looks
# its host up once and runs the search decision again; a transient that becomes
# ready for a targeted search is announced (journal, mail). No search is started.
#   ~/raven-kn/.venv-fink            : environment with fink-client (docs/listener.md)
#   ~/.finkclient/<survey>_credentials.yml : finkctl auth register -survey <survey> ...
#   ~/.config/raven_kn_brokers.env   : optional, LASAIR_TOKEN=...
#   ~/.config/raven_kn_notify.env    : optional, RAVEN_KN_MAIL_TO=address[,address]
#   RAVEN_KN_FINK_MAX_NEW            : new transients accepted per UTC day (default 50)
survey=$1; shift
for f in raven_kn_brokers raven_kn_notify; do
    [ -f ~/.config/$f.env ] && { set -a; . ~/.config/$f.env; set +a; }
done
cd ~/raven-kn-run
while true; do
    echo "$(date -Is) starting fink-$survey on $(hostname)" >> fink_$survey.log
    ~/raven-kn/.venv-fink/bin/raven-kn fink-listen --survey "$survey" --store kn_store/ --assess --find-host \
        --max-new-per-night "${RAVEN_KN_FINK_MAX_NEW:-50}" ${@:+--topics "$@"} \
        >> fink_$survey.jsonl 2>> fink_$survey.log
    echo "$(date -Is) fink-$survey exited ($?), restarting in 60 s" >> fink_$survey.log
    sleep 60
done
