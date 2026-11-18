#!/bin/bash
#   job.sh search <free|circular|distance|distance_samples>     job.sh pe <em|gw|em_samples>
# The *_samples modes read the request in samples/, next to the joint distance-inclination samples
# of the kilonova posterior, which then replace the intervals of the request.
set -e
source /cvmfs/software.igwn.org/conda/etc/profile.d/conda.sh
conda activate igwn
export HDF5_USE_FILE_LOCKING=FALSE OMP_NUM_THREADS=1 PYTHONPATH=$PWD
export HOME=${_CONDOR_SCRATCH_DIR:-$PWD} NUMBA_CACHE_DIR=${_CONDOR_SCRATCH_DIR:-$PWD}/numba_cache
mkdir -p $NUMBA_CACHE_DIR
G=GWOSC_4KHZ_R1-1187006835-4096.gwf
H=H1=H-H1_$G:H1:GWOSC-4KHZ_R1_STRAIN
V=V1=V-V1_$G:V1:GWOSC-4KHZ_R1_STRAIN
R=SSS17a_targeted_search_request.json; LABEL=SSS17a
case $2 in *_samples) R=samples/$R; LABEL=SSS17a_samples ;; esac
REQ="--request $R --statistic coherent --mchirp-min 1.19 --mchirp-max 1.205"
if [ "$1" = search ]; then
    case $2 in
        free)     OPT="--polarisation free" ;;
        circular) OPT="--polarisation circular" ;;
        distance|distance_samples) OPT="--polarisation circular --distance-prior" ;;
    esac
    python -m raven_kn.cli targeted $REQ --files $H,L1=L-L1_$G:L1:GWOSC-4KHZ_R1_STRAIN,$V \
        --n-slides 0 --n-offsource 100 $OPT --label search_$2 --out .
else
    # candidate of the unrestricted search (time and chirp mass); L1 with the glitch subtracted
    python -m raven_kn.cli targeted $REQ --files $H,L1=L-L1_LOSC_CLN_4_V1-1187007040-2048.gwf:L1:LOSC-STRAIN,$V \
        --pe ${2%_samples} --pe-result search_free_targeted_result.json --pe-npool ${3:-8} --label $LABEL --out .
fi
