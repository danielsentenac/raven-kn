#!/bin/bash
# Virgo leg, at CC-IN2P3 (login node): Virgo's raw frames hold no strain (h(t)
# needs Hrec), so the reconstructed O2 h(t) archived on the Lyon tapes is
# used. The frame covering the window is copied from tape, the strain around
# the window is cut out into a small frame (V-V1_RAVENKN_STRAIN-*.gwf) to be
# carried to CIT, and the staged copy is removed.
# Needs the raven_kn package next to this script (PYTHONPATH) and the IGWN
# conda environment from CVMFS.
set -e
START=${1:-1187008624}; END=${2:-1187009141}
TAPE=/hpss/in2p3.fr/group/virgo/Run/O2/proc/V2/V-V1O2Repro2A-1187005000-5000.gwf
WORK=${WORK:-/sps/virgo/USERS/$USER/raven_kn_virgo_hoft}
HERE=$(cd $(dirname $0) && pwd)
mkdir -p $WORK && cd $WORK
date; rfcp $TAPE . | tail -1; date
source /cvmfs/software.igwn.org/conda/etc/profile.d/conda.sh
conda activate igwn
PYTHONPATH=$HERE python - <<XX
import json
from raven_kn.targeted.datasource import extract_strain
e = extract_strain(["$(basename $TAPE)"], "V1:Hrec_hoft_V1O2Repro2A_16384Hz", "V1", $START, $END, "$HERE/frames_virgo")
e["original_frames"] = ["$TAPE"]
json.dump(e, open("$HERE/virgo_entry.json", "w"), indent=1)
print(e)
XX
rm -f $(basename $TAPE)
ls -la $HERE/frames_virgo
