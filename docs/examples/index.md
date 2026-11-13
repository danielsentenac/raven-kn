# Examples

Results on real events and on synthetic kilonovae. Every number comes from files committed under
[`examples/`](https://github.com/danielsentenac/raven-kn/tree/main/examples); the figures are remade with
`python docs/scripts/make_figures.py`.

| Example | Data | Shows |
|---|---|---|
| [AT2017gfo and GRB 211211A light-curve fits](kilonova-fits.md) | real kilonova photometry, NMMA | the merger-time window on real kilonovae, and why it needs a margin |
| [Injection calibration](injections.md) | 86 synthetic kilonovae, 201 fits, ZTF / Rubin / both | how often the window contains the merger, and how wide it is per survey |
| [GW170817 targeted search](gw170817-search.md) | GWOSC strain | the light search and PyGRB recover GW170817 at the kilonova position |
| [GW170817 from the GRB and the kilonova](gw170817-reverse.md) | GBM notice, real AT2017gfo light curve, GWOSC strain | GRB stored, kilonova confirmed by its light curve on the second night, then the 6 s search finds GW170817: reweighted coherent SNR 30.3, off-source maximum 7.6, p < 0.01 |
| [Host galaxies and redshifts](hosts.md) | NED, GLADE+, DES Y6 | host and distance for AT2017gfo, GRB 211211A, AT2025ulz |
| [S250818k / AT2025ulz](s250818k.md) | public LVK sky maps | joint FAR and distance tension for a real candidate pair |

## Synthetic candidates (`examples/make_examples.py`)

`raven-kn assess` on three made-up candidates with synthetic photometry and O4-like ranges:

| Candidate | Window | D_L (90 %) | θ_obs | Expected network SNR | Search |
|---|---|---|---|---|---|
| AT2017gfo-like, non-detection 3 h before merger | 19.4 h (bracket) | 44.6 [37, 52] Mpc | 21° [13, 29] | ~51 | yes |
| same + GRB | 6 s | same | same | ~68 | yes |
| z = 0.15 host | 7 d (no constraining limit) | 737 Mpc | unconstrained | ~2 | no: beyond 500 Mpc, P(SNR ≥ 6) = 0 |
