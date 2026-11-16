#!/bin/bash
# Step 2 (HTCondor job): light search around one event.
#   $1 event   $2 data kind (raw | hoft)   $3 bank (tmpl | equal)   $4 statistic (coherent | network)
#   $5 (optional) polarisation of the coherent statistic: circular | face_on | face_away
# tmpl : templates with the mass ratio and spins of the pipeline event, chirp mass +-1 %
# equal: equal-mass non-spinning templates (the pipeline default), chirp mass +-1 %
set -e
source /cvmfs/software.igwn.org/conda/etc/profile.d/conda.sh
conda activate igwn
export HDF5_USE_FILE_LOCKING=FALSE OMP_NUM_THREADS=1 PYTHONPATH=$PWD
export HOME=${_CONDOR_SCRATCH_DIR:-$PWD} NUMBA_CACHE_DIR=${_CONDOR_SCRATCH_DIR:-$PWD}/numba_cache
mkdir -p $NUMBA_CACHE_DIR
eval $(python - "$1" "$3" <<'PY'
import json, sys
e = json.load(open("events.json"))[sys.argv[1]]
q = e["m1"] / e["m2"] if sys.argv[2] == "tmpl" else 1.0
s1, s2 = (e["spin1z"], e["spin2z"]) if sys.argv[2] == "tmpl" else (0.0, 0.0)
print(f'RA={e["ra"]} DEC={e["dec"]} A={e["t0"] + e["window"][0]} B={e["t0"] + e["window"][1]} IFOS={e["ifos"]} '
      f'FLOW={e["f_low"]} MCA={e["mchirp"] * 0.99} MCB={e["mchirp"] * 1.01} Q={q} S1={s1} S2={s2}')
PY
)
python -m raven_kn.cli targeted --ra $RA --dec $DEC --gps-start $A --gps-end $B --ifos $IFOS \
    --mchirp-min $MCA --mchirp-max $MCB --mass-ratio $Q --spin1z=$S1 --spin2z=$S2 --f-low $FLOW \
    --manifest manifest_$1_$2.json --frame-dir frames_$1_$2 --statistic $4 --n-slides 0 --n-offsource 20 \
    --polarisation ${5:-free} \
    --label $1_$2_$3_$4${5:+_$5} --out .
