#!/bin/bash
# RAVEN-KN on CIT: PyGRB's coherent search (pycbc_multi_inspiral, LVK O2
# GRB170817A test arguments, no injections) at the RAVEN-KN sky position on
# the collaboration's O2 frames (cleaned C02 h(t) for H1 and L1, V1O2Repro2A);
# $1 = bank piece (0-3). Inputs arrive by HTCondor file transfer in the job's
# working directory (the frames listed by prepare.sh in frames.txt), where
# the frame caches are rebuilt with local paths.
set -e
source /cvmfs/software.igwn.org/conda/etc/profile.d/conda.sh
conda activate igwn
export HDF5_USE_FILE_LOCKING=FALSE OMP_NUM_THREADS=1
# HTCondor jobs start without HOME; numba must cache outside read-only CVMFS
export HOME=${_CONDOR_SCRATCH_DIR:-$PWD} NUMBA_CACHE_DIR=${_CONDOR_SCRATCH_DIR:-$PWD}/numba_cache
mkdir -p $NUMBA_CACHE_DIR
for ifo in H1 L1 V1; do
  awk -v d="$PWD" '{n = split($5, p, "/"); print $1, $2, $3, $4, "file://localhost" d "/" p[n]}' $ifo.cache > $ifo.local.cache
done
pycbc_multi_inspiral \
  --frame-cache H1:H1.local.cache L1:L1.local.cache V1:V1.local.cache \
  --channel-name H1:DCH-CLEAN_STRAIN_C02 L1:DCH-CLEAN_STRAIN_C02 V1:Hrec_hoft_V1O2Repro2A_16384Hz \
  --verbose --strain-high-pass 25 --sample-rate 2048 --segment-length 256 \
  --segment-start-pad 113 --segment-end-pad 17 --psd-estimation median --psd-segment-length 32 \
  --psd-segment-stride 8 --psd-num-segments 29 --autogating-threshold 100 --autogating-cluster 0.5 \
  --autogating-width 0.25 --autogating-taper 0.25 --autogating-pad 0 --low-frequency-cutoff 30 \
  --approximant IMRPhenomD --order -1 --cluster-method window --cluster-window 0.1 \
  --chisq-bins "0.72*get_freq('fIMRPhenomDPeak',params.mass1,params.mass2,params.spin1z,params.spin2z)**0.7" \
  --sngl-snr-threshold 4. --coinc-threshold 6.0 --chisq-index 6.0 --chisq-nhigh 2.0 --do-null-cut \
  --null-min 5.25 --null-grad 0.2 --null-step 20. --processing-scheme mkl --pad-data 8 \
  --projection left+right --instruments H1 L1 V1 \
  --gps-start-time 1187006900 --gps-end-time 1187010900 \
  --trig-start-time 1187007013 --trig-end-time 1187010883 --trigger-time 1187008794 \
  --bank-file BANK${1}.hdf --sky-grid sky_grid.h5 \
  --output INSPIRAL_BANK${1}.hdf
