#!/bin/bash
# $1 = raw | hoft : LIGO frames; Virgo h(t) frame the same in both
k=$1; R=$SCRATCH/rawrun; L=$R/frames_raw3; [ $k = hoft ] && L=$R/frames_hoft
export OMP_NUM_THREADS=2 HDF5_USE_FILE_LOCKING=FALSE
pycbc_multi_inspiral \
  --frame-files H1:$L/H-H1_RAVENKN_STRAIN-1187008624-517.gwf L1:$L/L-L1_RAVENKN_STRAIN-1187008624-517.gwf V1:$R/frames_raw3/V-V1_RAVENKN_STRAIN-1187008624-517.gwf \
  --channel-name H1:RAVENKN_STRAIN L1:RAVENKN_STRAIN V1:RAVENKN_STRAIN \
  --verbose --strain-high-pass 25 --sample-rate 2048 --segment-length 256 \
  --segment-start-pad 113 --segment-end-pad 17 --psd-estimation median --psd-segment-length 32 \
  --psd-segment-stride 8 --psd-num-segments 29 --autogating-threshold 100 --autogating-cluster 0.5 \
  --autogating-width 0.25 --autogating-taper 0.25 --autogating-pad 0 --low-frequency-cutoff 30 \
  --approximant IMRPhenomD --order -1 --cluster-method window --cluster-window 0.1 \
  --chisq-bins "0.72*get_freq('fIMRPhenomDPeak',params.mass1,params.mass2,params.spin1z,params.spin2z)**0.7" \
  --sngl-snr-threshold 4. --coinc-threshold 6.0 --chisq-index 6.0 --chisq-nhigh 2.0 --do-null-cut \
  --null-min 5.25 --null-grad 0.2 --null-step 20. --processing-scheme cpu --pad-data 8 \
  --projection left+right --instruments H1 L1 V1 \
  --gps-start-time 1187008700 --gps-end-time 1187008956 \
  --trig-start-time 1187008813 --trig-end-time 1187008939 --trigger-time 1187008794 \
  --bank-file BANK1.hdf --sky-grid sky_grid.h5 \
  --output $k/INSPIRAL_BANK1.hdf
