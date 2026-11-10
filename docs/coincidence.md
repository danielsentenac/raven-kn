# Distance, rate and joint FAR

## Distance and inclination

The distance decides whether a search is worth it and sets the chirp-mass range and distance prior of the
GW search. Sources, from most to least reliable:

1. **Spectroscopic host redshift**: from TNS, a GCN circular, or NED's nearest galaxy with a spectroscopic
   redshift within 50 kpc projected. Converted to the CMB frame, with 300 km/s peculiar-velocity error, at
   Planck18 cosmology.
2. **Redshift-independent host distance** (NED-D: SBF, TRGB, Cepheids, SNe Ia, ...), added for hosts at
   z < 0.03. For NGC 4993 this gives 40.7 ± 1.5 Mpc (SBF).
3. **Photometric redshift**: GLADE+ (LVK's galaxy catalogue), then DES Y6 Gold in its southern footprint.
   The host is the galaxy with the lowest chance-coincidence probability.
4. **Light-curve fit** distance: unreliable. For GRB 211211A it gives 165 Mpc (90 %: 115–274) against 357 Mpc
   from the host.

The **inclination** comes from external measurements (afterglow modelling, VLBI centroid motion,
light-curve fits), passed as `ViewingAngleMeasurement(mean, std, method)`. A prompt GRB adds a soft
θ_obs ≲ 30° constraint. Spectra are not used for it.

`raven_kn.treatment.distance_inclination` combines all of these into a joint (D_L, cos ι) posterior.

## Kilonova candidate rate

R_KN plays the role of RAVEN's GRB rate: how often a candidate this convincing turns up by chance.
RAVEN-KN applies no classifier itself; it takes what upstream filters give:

| Candidate | R_KN |
|---|---|
| spectroscopically classified kilonova | 5 / yr |
| with a classifier score and the filter's score distribution | 30 / month × fraction of candidates scoring as high (floor 1 / yr) |
| otherwise | 30 / month (everything the filter passes) |

Both numbers are placeholders to tune from the actual filter and TNS statistics. A spectroscopic
classification lowers the joint FAR by a factor of 72.

## Joint FAR

For a GW alert with FAR `FAR_GW` whose merger time falls inside the window, the joint FAR follows RAVEN's
form, extended with a distance term:

$$
\mathrm{FAR}_\mathrm{joint} = \frac{(t_h - t_l)\; R_\mathrm{KN}\; \mathrm{FAR}_\mathrm{GW}}{I_\mathrm{sky}\; I_\mathrm{dist}}
$$

- **I_sky**: the overlap of the GW sky map with the kilonova position (RAVEN's own overlap function is
  used if installed);
- **I_dist**: the overlap of the GW distance distribution along the candidate's direction with the EM
  distance, relative to the network horizon. Valid only for GW maps with distance layers, and only for
  hosts inside the reference horizon.

The temporal, sky and sky + distance FARs are all reported. `raven_kn.search.calc_signif_gracedb` can
upload the result to GraceDB (log prefixed `RAVEN-KN:`, tag `ext_coinc`, file `coincidence_far_kn.json`),
only after vetting.
