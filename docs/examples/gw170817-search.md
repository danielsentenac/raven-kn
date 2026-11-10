# GW170817 targeted search

GW170817 is searched for as if only its kilonova, AT2017gfo, were known: the sky position of NGC 4993, a
window around the merger, and a chirp-mass range of 1.18–1.22 M☉ (181 templates). Strain data from GWOSC.

## Light search

Network statistic, 1000 time slides, with H1 and V1 from GWOSC and L1 from the glitch-subtracted release
(LIGO-P1700349):

| Window | Network SNR | H1 / L1 / V1 | Time, chirp mass | p-value | Background max / median |
|---|---|---|---|---|---|
| containing the merger | **30.9** (published 32.4) | 18.6 / 24.5 / 3.2 | 1187008882.426, 1.1976 M☉ | < 1/1001 | 8.5 / 7.5 |
| 15 min earlier (noise only) | 7.4 | 4.0 / 5.4 / 3.1 | — | 0.62 | 8.5 / 7.5 |

On the frames with the L1 glitch still in (GWTC-1 release), the glitch 1.1 s before the merger is gated:

![GW170817 against the time-slide background](../img/gw170817_background.png)

*Loudest event in the window containing GW170817 (white line) against the loudest event of each of 1000
time slides, for the network and coherent statistics, with PyGRB's gate on the same frames.*

GW170817 stands far above every slide with both statistics. The coherent statistic lowers the background
by about one unit (7.6 against 8.5 at most) because of its null-stream cut.

## Comparison with PyGRB

PyGRB (`pycbc_multi_inspiral` with LVK's O2 GRB170817A arguments) ran on the same frames, with the
RAVEN-KN template bank split in 4 HTCondor jobs on CIT (4 min each). Both gate the L1 glitch at the same
instant.

| | Coherent SNR | H1 / L1 / V1 | Null SNR |
|---|---|---|---|
| Light search, raw SNR | **30.4** | 18.4 / 24.0 / 3.3 | 1.0 |
| Light search, newSNR (per-detector chi-square) | 24.4 | 18.4 / 15.9 / 3.2 | 2.5 |
| PyGRB, raw / reweighted | **30.9** / 29.8 | 19.0 / 24.5 / 3.6 | 4.4 |

- Both find the same chirp mass (1.1976 M☉) and the same time to 1–4 ms.
- Raw SNRs agree to 2 %.
- The light search's per-detector chi-square cuts the gated L1 signal from 24 to 16, where PyGRB's coherent
  chi-square costs 3.5 %. The light search therefore under-ranks signals near a gate (still p < 1/1001
  here). A gate-aware chi-square is the next fix.

## Data-quality vetoes

DQSegDB queried live from CIT (O2 CBC veto definer, ±600 s around the merger): H1, L1 and V1 in science
throughout, no category-1 time, one 2 s category-2 interval 444 s after the merger, outside the on-source
window. That interval stopped a first PyGRB attempt with a longer window, which is why sub-segments are
now placed in veto-free stretches.

Results: `examples/targeted_gw170817/`, `examples/pygrb_gw170817/pygrb_loudest.json`.
