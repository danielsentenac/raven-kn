#!/bin/bash
# RAVEN-KN light search on LVK frames, step 1 (login node, where the SciToken
# and datafind are): resolve the science segments, CAT1/CAT2 vetoes (DQSegDB)
# and the frame files of the collaboration's O2 h(t) for the search window,
# into manifest.json. The job (light_search.sh) then needs neither network
# nor credentials. Needs ~/raven-kn (RAVEN-KN in a venv on top of the igwn
# env) and LVK's O2 CBC veto definer.
set -e
VETO_DEFINER=${VETO_DEFINER:-$HOME/raven_kn_pygrb_test/H1L1V1-CBC_VETO_DEFINER_CLEANED_C02_O2_1164556817-23176801.xml}
source /cvmfs/software.igwn.org/conda/etc/profile.d/conda.sh
conda activate igwn
htgettoken -a vault.ligo.org -i igwn --nooidc
export GWDATAFIND_SERVER=datafind.ldas.cit:80
~/raven-kn/.venv/bin/raven-kn targeted --request SSS17a_targeted_search_request.json \
    --source lvk --lvk-preset O2 --veto-definer $VETO_DEFINER --write-manifest manifest.json
# frame paths for HTCondor file transfer
python -c "
import json
m = json.load(open('manifest.json'))
print(','.join(sorted({u.replace('file://localhost', '').replace('file://', '') for i in m['ifos'].values() for u in i.get('frames', [])})))" > frames.txt
rm -rf raven_kn && cp -r ~/raven-kn/raven_kn . && find raven_kn -name __pycache__ -prune -exec rm -r {} +
echo "frames: $(tr ',' '\n' < frames.txt | wc -l)"
