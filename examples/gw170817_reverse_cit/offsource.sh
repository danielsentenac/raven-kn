#!/bin/bash
# RAVEN-KN GW170817 reverse replay on CIT: off-source trials $1 .. $1+$2-1 of
# the 6 s GRB-window search (coherent statistic, 1.0-1.5 Msun bank), on the
# GWOSC GWTC-1 frames. Inputs arrive by HTCondor file transfer.
set -e
source /cvmfs/software.igwn.org/conda/etc/profile.d/conda.sh
conda activate igwn
export HDF5_USE_FILE_LOCKING=FALSE OMP_NUM_THREADS=1 PYTHONPATH=$PWD
# HTCondor jobs start without HOME; numba must cache outside read-only CVMFS
export HOME=${_CONDOR_SCRATCH_DIR:-$PWD} NUMBA_CACHE_DIR=${_CONDOR_SCRATCH_DIR:-$PWD}/numba_cache
mkdir -p $NUMBA_CACHE_DIR
python gw170817_reverse_offsource.py run --request SSS17a_targeted_search_request.json \
  --first $1 --n $2 --statistic coherent \
  --frames H1=H-H1_GWOSC_4KHZ_R1-1187006835-4096.gwf:H1:GWOSC-4KHZ_R1_STRAIN \
           L1=L-L1_GWOSC_4KHZ_R1-1187006835-4096.gwf:L1:GWOSC-4KHZ_R1_STRAIN \
           V1=V-V1_GWOSC_4KHZ_R1-1187006835-4096.gwf:V1:GWOSC-4KHZ_R1_STRAIN \
  --out chunk_$1.json
