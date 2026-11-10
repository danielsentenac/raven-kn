# Host galaxies and redshifts

`raven-kn host` finds the host galaxy and its redshift: NED first (best spectroscopic completeness beyond
100 Mpc), then GLADE+, then DES Y6 Gold in its footprint. For hosts at z < 0.03 it adds NED-D's best
redshift-independent distance.

| Event | Host found | Truth |
|---|---|---|
| AT2017gfo | NGC 4993 (NED), z_cmb 0.0109 spectroscopic, 2.2 kpc projected; SBF 40.7 ± 1.5 Mpc (NED-D) | 40.7 Mpc (SBF) |
| GRB 211211A | WISEA J140910.48+275320.6 (NED), z 0.0763 spectroscopic | z 0.0763 (not in GLADE+) |
| AT2025ulz | GLADE+ 15877920, 0.9″, **photometric** z 0.102 ± 0.036 (no host redshift in NED yet) | DESI / Keck z 0.0848 |
| ESO 243-49, transient placed 8″ off | DES Y6 galaxy at 7.5″, photo-z 0.022 ± 0.013, chance-coincidence probability 8 × 10⁻⁴ | z 0.0224 |

## DES Y6 photometric redshifts

DES Y6 Gold covers 5000 deg² of the southern sky with about 450 million galaxies, 20 times more than GLADE+
there. Against 109 NED spectroscopic galaxies at 0.01 < z < 0.15 in 7 DES fields:

| Metric | Value |
|---|---|
| scatter, NMAD of Δz / (1 + z) | 0.013 |
| bias | < 0.007 |
| outliers | 1 % |

The catalogue's own error (`dnf_zsigma`) is 1.2–1.5 times too small, so RAVEN-KN uses 1.5 × `dnf_zsigma`,
with a floor of 0.01.

## In practice

For a new event, the spectroscopic host redshift usually arrives hours to days later, through TNS or GCN
circulars, and replaces a photometric one. AT2025ulz shows why it matters: the photometric z 0.102 ± 0.036 puts
the host at about 490 Mpc with a wide spread, close to the 500 Mpc cut, while the spectroscopic z 0.0848 gives
400 Mpc.
