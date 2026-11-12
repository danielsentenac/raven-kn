#!/bin/bash
# RAVEN-KN light search on LVK frames, step 2 (HTCondor job): the 6 s GRB
# window of the GW170817 reverse replay, coherent statistic, 1.0-1.5 Msun
# bank, on the collaboration's O2 frames resolved by prepare.sh. Inputs
# (raven_kn, request, manifest, frames) arrive by HTCondor file transfer.
set -e
source /cvmfs/software.igwn.org/conda/etc/profile.d/conda.sh
conda activate igwn
export HDF5_USE_FILE_LOCKING=FALSE OMP_NUM_THREADS=1 PYTHONPATH=$PWD
# HTCondor jobs start without HOME; numba must cache outside read-only CVMFS
export HOME=${_CONDOR_SCRATCH_DIR:-$PWD} NUMBA_CACHE_DIR=${_CONDOR_SCRATCH_DIR:-$PWD}/numba_cache
mkdir -p $NUMBA_CACHE_DIR
python -m raven_kn.cli targeted --request SSS17a_targeted_search_request.json \
    --manifest manifest.json --statistic coherent --n-slides 0 \
    --label SSS17a_GRB_window_coherent_lvk --out .
