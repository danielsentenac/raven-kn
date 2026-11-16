#!/bin/bash
set -e
source /cvmfs/software.igwn.org/conda/etc/profile.d/conda.sh
conda activate igwn
export HDF5_USE_FILE_LOCKING=FALSE OMP_NUM_THREADS=1 PYTHONPATH=$PWD
export HOME=${_CONDOR_SCRATCH_DIR:-$PWD} NUMBA_CACHE_DIR=${_CONDOR_SCRATCH_DIR:-$PWD}/numba_cache
mkdir -p $NUMBA_CACHE_DIR
python run.py $1
