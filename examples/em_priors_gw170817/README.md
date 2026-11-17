# GW170817: what the EM priors bring to the GW search

HTCondor jobs (CIT) on public data, from the request of the reverse replay
(`SSS17a_targeted_search_request.json`: AT2017gfo position, 30-54 Mpc, viewing angle below 58 deg).

- `prepare.sh`: links the GWOSC 4 kHz frames and downloads the L1 frame with the glitch subtracted
  (LIGO-P1700349), used for the parameter estimation.
- `job.sh search <free|circular|distance>`: coherent search of the GRB window, 100 off-source trials:
  unrestricted, circular polarisation, circular polarisation + distance prior (`--distance-prior`).
- `job.sh pe <em|gw>`: Bilby on the candidate of the `free` search, sky position fixed, with the EM
  distance and inclination priors (`em`) or without (`gw`).
- `jobs.sub`: the three searches; `pe.sub`: the two parameter estimations (8 cores each), to submit
  once the `free` search result is back.
