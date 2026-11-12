#!/bin/bash
# RAVEN-KN on CIT, LVK data path: GW170817 with PyGRB on the collaboration's
# O2 frames (not the GWOSC ones). Run on the login node before submitting:
#   1. IGWN SciToken (non-interactive when a vault credential exists),
#   2. DQSegDB: science segments and CAT1/CAT2 vetoes of LVK's O2 CBC veto
#      definer -> segments.json,
#   3. datafind (CIT's local server): frame caches H1.cache, L1.cache, V1.cache.
# Needs ~/raven-kn (RAVEN-KN in a venv on top of the igwn env), the veto
# definer XML, and bank/ + sky_grid.h5 as in examples/pygrb_gw170817.
set -e
START=1187006900; END=1187010900; PAD=8
VETO_DEFINER=${VETO_DEFINER:-$HOME/raven_kn_pygrb_test/H1L1V1-CBC_VETO_DEFINER_CLEANED_C02_O2_1164556817-23176801.xml}
source /cvmfs/software.igwn.org/conda/etc/profile.d/conda.sh
conda activate igwn
htgettoken -a vault.ligo.org -i igwn --nooidc
~/raven-kn/.venv/bin/python lvk_segments.py --start $START --end $END --veto-definer $VETO_DEFINER \
    --science H1=H1:DCH-CLEAN_SCIENCE_C02:1 L1=L1:DCH-CLEAN_SCIENCE_C02:1 V1=V1:ITF_SCIENCE:2 \
    --out segments.json
export GWDATAFIND_SERVER=datafind.ldas.cit:80
for spec in H1:H1_CLEANED_HOFT_C02 L1:L1_CLEANED_HOFT_C02 V1:V1O2Repro2A; do
  ifo=${spec%%:*}; type=${spec##*:}
  gw_data_find -o ${ifo:0:1} -t $type -s $((START - PAD)) -e $((END + PAD)) -u file --lal-cache > $ifo.cache
  echo "$ifo $type: $(wc -l < $ifo.cache) frame files"
done
# frame paths for HTCondor file transfer
cat H1.cache L1.cache V1.cache | awk '{sub("file://localhost", "", $5); print $5}' | paste -sd, > frames.txt
