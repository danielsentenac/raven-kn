#!/bin/bash
# Run once in the submit directory on CIT: GWOSC frames (already used by the reverse replay) and
# the glitch-subtracted L1 frame of GW170817 (public, LIGO-P1700349).
set -e
ln -sfn ${GWOSC_FRAMES:-$HOME/raven_kn_pygrb/frames} frames
mkdir -p logs clean
[ -s clean/L-L1_LOSC_CLN_4_V1-1187007040-2048.gwf ] || curl -sSL -o clean/L-L1_LOSC_CLN_4_V1-1187007040-2048.gwf \
    https://dcc.ligo.org/public/0146/P1700349/001/L-L1_LOSC_CLN_4_V1-1187007040-2048.gwf
ls -la frames/ clean/
