#!/bin/bash
# RAVEN-KN light search on LIGO RAW frames, step 2 (HTCondor job): the 6 s
# GRB window of the GW170817 reverse replay, H1 L1, coherent statistic,
# 1.0-1.5 Msun bank. $1 = raw (strain from the RAW front-end channel) or
# hoft (final O2 h(t), for comparison); frames extracted by prepare.sh.
set -e
source /cvmfs/software.igwn.org/conda/etc/profile.d/conda.sh
conda activate igwn
export HDF5_USE_FILE_LOCKING=FALSE OMP_NUM_THREADS=1 PYTHONPATH=$PWD
# HTCondor jobs start without HOME; numba must cache outside read-only CVMFS
export HOME=${_CONDOR_SCRATCH_DIR:-$PWD} NUMBA_CACHE_DIR=${_CONDOR_SCRATCH_DIR:-$PWD}/numba_cache
mkdir -p $NUMBA_CACHE_DIR
IFOS=H1,L1; [ $1 = raw3 ] && IFOS=H1,L1,V1      # raw3: + Virgo h(t) from the Lyon tapes
python -m raven_kn.cli targeted --request SSS17a_targeted_search_request.json --ifos $IFOS \
    --manifest manifest_$1.json --frame-dir frames_$1 --statistic coherent --n-slides 0 \
    --label SSS17a_GRB_window_$1 --out .
