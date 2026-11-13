#!/bin/bash
# RAVEN-KN light search on LIGO RAW frames, step 1 (CIT login node). The RAW
# frames (H1_R, L1_R: 2 GB per 64 s, CIT archive only, not on OSDF) hold no
# h(t): the front-end channel CAL-DELTAL_EXTERNAL_DQ is read, converted to
# strain and written to small frame files, with the low-latency science
# segments and the CAT1/CAT2 vetoes from DQSegDB -> manifest_raw.json.
# For comparison the same two detectors on the final O2 h(t) ->
# manifest_hoft.json. No Virgo: its raw data are on tape at CC-IN2P3.
# Needs ~/raven-kn (RAVEN-KN in a venv on top of the igwn env) and LVK's O2
# CBC veto definer.
set -e
VETO_DEFINER=${VETO_DEFINER:-$HOME/raven_kn_pygrb_test/H1L1V1-CBC_VETO_DEFINER_CLEANED_C02_O2_1164556817-23176801.xml}
source /cvmfs/software.igwn.org/conda/etc/profile.d/conda.sh
conda activate igwn
htgettoken -a vault.ligo.org -i igwn --nooidc
export GWDATAFIND_SERVER=datafind.ldas.cit:80
for kind in raw hoft; do
  preset=O2; [ $kind = raw ] && preset=O2-raw
  ~/raven-kn/.venv/bin/raven-kn targeted --request SSS17a_targeted_search_request.json --ifos H1,L1 \
      --source lvk --lvk-preset $preset --veto-definer $VETO_DEFINER \
      --write-manifest manifest_$kind.json --extract-frames frames_$kind
done
rm -rf raven_kn && cp -r ~/raven-kn/raven_kn . && find raven_kn -name __pycache__ -prune -exec rm -r {} +
ls -la frames_raw frames_hoft
