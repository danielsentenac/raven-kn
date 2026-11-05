"""Joint (luminosity distance, inclination) posterior from EM information.

Each piece of EM information is a likelihood on (D_L, cos iota). They are
multiplied over a prior uniform in volume (D_L^2) and isotropic in
orientation (cos iota uniform), by importance sampling. The result

* tells whether a GW search is worth it (``treatment.detectability``),
* gives the distance/inclination priors of the targeted search and of a
  later joint parameter estimation (the D_L--iota degeneracy of the GW
  data is what an EM inclination breaks; GW170817: Hotokezaka+2019),
* gives the distance term of the joint FAR (``skymap.distance_overlap``).

Conventions: iota in [0, pi] is the angle between the orbital angular
momentum and the line of sight; the EM viewing angle is
theta_obs = min(iota, pi - iota) because the EM emission is symmetric.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
from astropy import units as u
from astropy.cosmology import Planck18, z_at_value

from ..connection import KNCandidate

C_KMS = 299792.458
LogLike = Callable[[np.ndarray, np.ndarray], np.ndarray]   # (d_mpc, cos_iota) -> log L


# --------------------------------------------------------------------------
# Likelihood terms
# --------------------------------------------------------------------------

@dataclass
class Constraint:
    name: str
    loglike: LogLike


def _z_of_d(cosmo=Planck18, d_max=3000.0, n=2000):
    """Interpolator D_L [Mpc] -> z (cosmological)."""
    zmax = float(z_at_value(cosmo.luminosity_distance, d_max * u.Mpc))
    zz = np.linspace(0, zmax, n)
    dd = cosmo.luminosity_distance(zz).to_value(u.Mpc)
    return lambda d: np.interp(d, dd, zz)


def redshift_constraint(z: float, z_err: float = 1e-4, v_pec_err_kms: float = 300.0,
                        cosmo=Planck18, kind: str = "spectroscopic") -> Constraint:
    """Spectroscopic (host or transient) redshift. The peculiar velocity
    (~300 km/s) dominates below ~100 Mpc. The cosmology (H0) is fixed: for
    an H0 measurement use the GW distance instead."""
    z_of_d = _z_of_d(cosmo)
    sig = np.hypot(z_err, v_pec_err_kms / C_KMS)

    def ll(d, c):
        return -0.5 * ((z_of_d(d) - z) / sig) ** 2
    return Constraint(f"redshift({kind}) z={z}", ll)


def distance_constraint(mean_mpc: float, std_mpc: float, method: str = "") -> Constraint:
    """Redshift-independent host distance, or a distance from a light-curve
    fit summarised as mean/std."""
    def ll(d, c):
        return -0.5 * ((d - mean_mpc) / std_mpc) ** 2
    return Constraint(f"distance({method}) {mean_mpc}+-{std_mpc} Mpc", ll)


def viewing_angle_constraint(mean_deg: float, std_deg: float, method: str = "") -> Constraint:
    """Gaussian on theta_obs (radio/VLBI jet, afterglow, kilonova colour)."""
    m, s = np.radians(mean_deg), np.radians(std_deg)

    def ll(d, c):
        theta = np.arccos(np.abs(c))
        return -0.5 * ((theta - m) / s) ** 2
    return Constraint(f"viewing_angle({method}) {mean_deg}+-{std_deg} deg", ll)


def prompt_grb_constraint(theta_max_deg: float = 30.0, width_deg: float = 5.0) -> Constraint:
    """A detected prompt short GRB favours small viewing angles. Soft cut
    (logistic) at ``theta_max_deg``; GRB 170817A was seen at ~20 deg, so
    the default is permissive. Assumption to revisit with a jet-structure
    model."""
    tm, w = np.radians(theta_max_deg), np.radians(width_deg)

    def ll(d, c):
        theta = np.arccos(np.abs(c))
        return -np.logaddexp(0.0, (theta - tm) / w)
    return Constraint(f"prompt GRB theta<~{theta_max_deg} deg", ll)


def samples_constraint(d_samples: np.ndarray, cos_iota_samples: Optional[np.ndarray] = None,
                       name: str = "lc_fit") -> Constraint:
    """KDE likelihood from posterior samples of a light-curve fit (NMMA,
    Redback). Assumes the fit used broad priors, so posterior ~ likelihood;
    divide out the fit's prior first if it was informative."""
    from scipy.stats import gaussian_kde
    if cos_iota_samples is None:
        kde = gaussian_kde(np.asarray(d_samples))

        def ll(d, c):
            return np.log(kde(d) + 1e-300)
    else:
        # theta_obs symmetric: fold to |cos iota|
        kde = gaussian_kde(np.vstack([d_samples, np.abs(cos_iota_samples)]))

        def ll(d, c):
            return np.log(kde(np.vstack([d, np.abs(c)])) + 1e-300)
    return Constraint(name, ll)


def constraints_from_candidate(kn: KNCandidate) -> list[Constraint]:
    cons = []
    if kn.redshift is not None:
        cons.append(redshift_constraint(kn.redshift.z, kn.redshift.z_err,
                                        kind=kn.redshift.kind))
    if kn.distance is not None:
        cons.append(distance_constraint(kn.distance.mean_mpc, kn.distance.std_mpc,
                                        kn.distance.method))
    if kn.viewing_angle is not None:
        cons.append(viewing_angle_constraint(kn.viewing_angle.mean_deg,
                                             kn.viewing_angle.std_deg,
                                             kn.viewing_angle.method))
    if kn.grb is not None:
        cons.append(prompt_grb_constraint())
    if (kn.lc_fit_samples and "luminosity_distance" in kn.lc_fit_samples
            and not kn.lc_fit_informed):
        s = kn.lc_fit_samples
        ci = None
        if "inclination" in s:
            ci = np.cos(np.asarray(s["inclination"]))
        elif "cos_iota" in s:
            ci = np.asarray(s["cos_iota"])
        cons.append(samples_constraint(np.asarray(s["luminosity_distance"]), ci))
    return cons


# --------------------------------------------------------------------------
# Posterior
# --------------------------------------------------------------------------

@dataclass
class DistanceInclinationPosterior:
    d_mpc: np.ndarray
    cos_iota: np.ndarray
    constraints: list[str] = field(default_factory=list)
    n_eff: float = 0.0
    d_max_prior: float = 0.0

    @property
    def theta_obs_deg(self) -> np.ndarray:
        return np.degrees(np.arccos(np.abs(self.cos_iota)))

    def summary(self, cl: float = 0.9) -> dict:
        a = (1 - cl) / 2

        def q(x):
            lo, med, hi = np.quantile(x, [a, 0.5, 1 - a])
            return {"median": float(med), "lower": float(lo), "upper": float(hi)}
        return {"credible_level": cl,
                "luminosity_distance_mpc": q(self.d_mpc),
                "viewing_angle_deg": q(self.theta_obs_deg),
                "n_eff": float(self.n_eff),
                "constraints": self.constraints,
                "distance_prior": f"uniform in volume, D_L < {self.d_max_prior} Mpc",
                "inclination_prior": "isotropic (cos iota uniform)"}

    def search_priors(self, pad: float = 0.2, cl: float = 0.99) -> dict:
        """Bounded priors for the targeted search / joint PE: the ``cl``
        interval of D_L widened by ``pad``, and the viewing-angle interval.
        Both orientations (iota and pi - iota) are kept: EM data cannot tell
        them apart."""
        a = (1 - cl) / 2
        dlo, dhi = np.quantile(self.d_mpc, [a, 1 - a])
        tlo, thi = np.quantile(self.theta_obs_deg, [a, 1 - a])
        return {"luminosity_distance_mpc": [float(max(dlo * (1 - pad), 0.1)),
                                            float(dhi * (1 + pad))],
                "viewing_angle_deg": [float(max(tlo - 5, 0)), float(min(thi + 5, 90))],
                "inclination_ranges_rad": [[float(np.radians(max(tlo - 5, 0))),
                                            float(np.radians(min(thi + 5, 90)))],
                                           [float(np.pi - np.radians(min(thi + 5, 90))),
                                            float(np.pi - np.radians(max(tlo - 5, 0)))]],
                "credible_level": cl, "pad": pad}


def _grid_post(constraints, d_lo, d_hi, n_d, n_c):
    d = np.linspace(d_lo, d_hi, n_d + 1)
    d = 0.5 * (d[1:] + d[:-1])
    c = np.linspace(-1, 1, n_c + 1)
    c = 0.5 * (c[1:] + c[:-1])
    D, C = np.meshgrid(d, c, indexing="ij")
    logp = 2 * np.log(D)                         # uniform-in-volume prior; cos iota flat
    for con in constraints:
        logp = logp + con.loglike(D.ravel(), C.ravel()).reshape(D.shape)
    p = np.exp(logp - logp.max())
    return d, c, p / p.sum()


def combine(constraints: list[Constraint], d_max: float = 1500.0, n_out: int = 20_000,
            n_d: int = 1500, n_c: int = 300, seed: Optional[int] = 0) -> DistanceInclinationPosterior:
    """Posterior samples of (D_L, cos iota) from the product of constraints.

    Two-stage grid (the problem is 2D): a coarse grid over [0, d_max] finds
    the support, a fine grid over it gives the samples (jittered within
    cells). With no constraint the result is the prior: the distance is
    then unknown, which the detectability step reports.
    """
    rng = np.random.default_rng(seed)
    names = [x.name for x in constraints]
    d, c, p = _grid_post(constraints, 0.0, d_max, n_d, 50)
    pd = p.sum(axis=1)
    cdf = np.cumsum(pd)
    step = d[1] - d[0]
    lo = max(d[np.searchsorted(cdf, 1e-7)] - 3 * step, 1e-3)
    hi = min(d[min(np.searchsorted(cdf, 1 - 1e-7), len(d) - 1)] + 3 * step, d_max)
    d, c, p = _grid_post(constraints, lo, hi, n_d, n_c)
    flat = p.ravel()
    idx = rng.choice(flat.size, size=n_out, p=flat)
    i, j = np.unravel_index(idx, p.shape)
    dd, dc = d[1] - d[0], c[1] - c[0]
    ds = d[i] + rng.uniform(-0.5, 0.5, n_out) * dd
    cs = np.clip(c[j] + rng.uniform(-0.5, 0.5, n_out) * dc, -1, 1)
    n_eff = 1.0 / np.sum(flat ** 2)              # effective number of grid cells
    return DistanceInclinationPosterior(ds, cs, names, n_eff, d_max)


def posterior_for(kn: KNCandidate, **kw) -> DistanceInclinationPosterior:
    return combine(constraints_from_candidate(kn), **kw)
