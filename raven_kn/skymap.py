"""Sky map module: kilonova sky maps and overlaps with GW sky maps.

* The kilonova is localised to arcseconds, so its 2D map is essentially one
  HEALPix pixel (RAVEN's ``use_radec`` case, as for Swift GRBs). A Gaussian
  of width ``pos_err`` is used when that is larger than a pixel.
* The distance posterior is attached as LVK 3D sky map columns (PROB,
  DISTMU, DISTSIGMA, DISTNORM; Singer et al. 2016 ansatz), so the product is
  readable by ligo.skymap tools and usable as a targeted-search input.
* :func:`skymap_overlap_integral` uses ``ligo.raven``'s own implementation
  when installed, so values are identical to RAVEN's; a flat-map fallback is
  provided otherwise.
* :func:`distance_overlap` is the distance analogue of the sky overlap: how
  much more probable the GW distance in the kilonova's direction is under
  the kilonova distance than under the prior.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import healpy as hp
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import norm

from .treatment.distance_inclination import DistanceInclinationPosterior


# --------------------------------------------------------------------------
# Kilonova sky map
# --------------------------------------------------------------------------

def kn_prob_map(ra_deg: float, dec_deg: float, pos_err_arcsec: float = 0.5,
                nside: int = 512) -> np.ndarray:
    """2D probability map (NESTED) for a point-like source."""
    npix = hp.nside2npix(nside)
    prob = np.zeros(npix)
    vec = hp.ang2vec(ra_deg, dec_deg, lonlat=True)
    sig = np.radians(pos_err_arcsec / 3600.)
    if sig < hp.nside2resol(nside) / 2:
        prob[hp.vec2pix(nside, *vec, nest=True)] = 1.0
        return prob
    ipix = hp.query_disc(nside, vec, 5 * sig, nest=True, inclusive=True)
    ang = np.arccos(np.clip(np.dot(np.array(hp.pix2vec(nside, ipix, nest=True)).T, vec), -1, 1))
    prob[ipix] = np.exp(-0.5 * (ang / sig) ** 2)
    return prob / prob.sum()


def ansatz_from_moments(mean: float, std: float) -> tuple[float, float, float]:
    """(DISTMU, DISTSIGMA, DISTNORM) such that p(r) = DISTNORM r^2 N(r; mu, sigma)
    (r > 0) has the given mean and standard deviation."""
    r = np.linspace(1e-3, mean + 12 * std, 6000)

    def moments(mu, sig):
        p = r ** 2 * norm.pdf(r, mu, sig)
        z = np.trapezoid(p, r)
        m = np.trapezoid(r * p, r) / z
        s = np.sqrt(np.trapezoid((r - m) ** 2 * p, r) / z)
        return m, s, z

    def resid(x):
        m, s, _ = moments(x[0], abs(x[1]))
        return [(m - mean) / std, (s - std) / std]

    sol = least_squares(resid, [mean - 2 * std ** 2 / mean, std])
    mu, sig = sol.x[0], abs(sol.x[1])
    return float(mu), float(sig), float(1 / moments(mu, sig)[2])


def kn_skymap_3d(ra_deg: float, dec_deg: float, post: Optional[DistanceInclinationPosterior],
                 pos_err_arcsec: float = 0.5, nside: int = 512) -> dict[str, np.ndarray]:
    prob = kn_prob_map(ra_deg, dec_deg, pos_err_arcsec, nside)
    out = {"PROB": prob}
    if post is not None:
        mu, sig, nrm = ansatz_from_moments(float(np.mean(post.d_mpc)),
                                           float(np.std(post.d_mpc)))
        on = prob > 0
        for k, v in (("DISTMU", mu), ("DISTSIGMA", sig), ("DISTNORM", nrm)):
            col = np.full(prob.shape, np.inf if k == "DISTMU" else (1.0 if k == "DISTSIGMA" else 0.0))
            col[on] = v
            out[k] = col
    return out


def write_skymap(path: str | Path, columns: dict[str, np.ndarray], object_name: str = "",
                 gps: Optional[float] = None, post: Optional[DistanceInclinationPosterior] = None,
                 instrument: str = "RAVEN-KN") -> Path:
    path = Path(path)
    names = [k for k in ("PROB", "DISTMU", "DISTSIGMA", "DISTNORM") if k in columns]
    units = {"PROB": "pix-1", "DISTMU": "Mpc", "DISTSIGMA": "Mpc", "DISTNORM": "Mpc-2"}
    hdr = [("OBJECT", object_name), ("INSTRUME", instrument), ("CREATOR", "raven_kn")]
    if gps is not None:
        hdr.append(("GPSTIME", gps))
    if post is not None:
        hdr += [("DISTMEAN", float(np.mean(post.d_mpc))), ("DISTSTD", float(np.std(post.d_mpc)))]
    hp.write_map(str(path), [columns[k] for k in names], nest=True, column_names=names,
                 column_units=[units[k] for k in names], extra_header=hdr,
                 overwrite=True, dtype=np.float64)
    return path


# --------------------------------------------------------------------------
# GW sky maps
# --------------------------------------------------------------------------

def read_gw_skymap(path: str | Path) -> dict:
    """Read a GW sky map as a flat NESTED map. Multi-order (MOC) files need
    ligo.skymap; flattened ``*.fits.gz`` maps are read with healpy."""
    try:
        from ligo.skymap.io import read_sky_map
        (prob, *dist), meta = read_sky_map(str(path), nest=True, distances=True)
        out = {"PROB": prob, "nested": True, "meta": meta}
        if dist:
            out.update(DISTMU=dist[0], DISTSIGMA=dist[1], DISTNORM=dist[2])
        return out
    except ImportError:
        pass
    maps, hdr = hp.read_map(str(path), field=None, nest=True, h=True)
    maps = np.atleast_2d(maps)
    out = {"PROB": maps[0], "nested": True, "meta": dict(hdr)}
    if maps.shape[0] >= 4:
        out.update(DISTMU=maps[1], DISTSIGMA=maps[2], DISTNORM=maps[3])
    return out


def skymap_overlap_integral(gw_prob: np.ndarray, ra: float, dec: float,
                            nested: bool = True) -> float:
    """RAVEN sky overlap for a point-like external position:
    p_GW(pixel) * Npix / sum(p_GW) (= 1 for an uninformative GW map)."""
    try:
        from ligo.raven.search import skymap_overlap_integral as raven_overlap
        return float(raven_overlap(gw_prob, ra=ra, dec=dec, se_nested=nested))
    except ImportError:
        nside = hp.npix2nside(len(gw_prob))
        ipix = hp.ang2pix(nside, ra, dec, lonlat=True, nest=nested)
        return float(gw_prob[ipix] * len(gw_prob) / np.sum(gw_prob))


def conditional_distance(gw: dict, ra: float, dec: float,
                         d_compare_mpc: Optional[float] = None) -> dict:
    """GW luminosity-distance posterior along one direction.

    The 3D map stores, per pixel, the *parameters* of the ansatz
    p(r) = DISTNORM r^2 N(r; DISTMU, DISTSIGMA): DISTMU/DISTSIGMA are not the
    mean and standard deviation of the distance. This returns the actual
    conditional mean, sd and 90 % interval, and optionally P(D > d_compare)
    and the distance likelihood ratio N(d_compare)/N(DISTMU) (the ansatz
    Gaussian is the GW distance likelihood once the r^2 prior is removed).
    """
    if "DISTMU" not in gw:
        raise KeyError("GW sky map has no distance layers")
    nside = hp.npix2nside(len(gw["PROB"]))
    ipix = hp.ang2pix(nside, ra, dec, lonlat=True, nest=gw.get("nested", True))
    mu, sig = float(gw["DISTMU"][ipix]), float(gw["DISTSIGMA"][ipix])
    out = {"distmu": mu, "distsigma": sig}
    if not np.isfinite(mu) or sig <= 0:
        return out
    r = np.linspace(1e-3, max(mu, 0) + 12 * sig, 20000)
    p = r ** 2 * norm.pdf(r, mu, sig)
    p /= np.trapezoid(p, r)
    cdf = np.concatenate([[0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(r))])
    mean = float(np.trapezoid(r * p, r))
    out.update(mean_mpc=mean, sd_mpc=float(np.sqrt(np.trapezoid((r - mean) ** 2 * p, r))),
               interval90_mpc=[float(x) for x in np.interp([0.05, 0.95], cdf, r)])
    if d_compare_mpc is not None:
        out["p_beyond"] = float(1 - np.interp(d_compare_mpc, r, cdf))
        out["likelihood_ratio_to_peak"] = float(np.exp(-0.5 * ((d_compare_mpc - mu) / sig) ** 2))
    return out


def distance_overlap_valid(post: DistanceInclinationPosterior, d_ref_mpc: float) -> bool:
    """The distance factor assumes the kilonova distance lies inside the
    reference volume (unrelated GW candidates uniform in volume up to
    ``d_ref_mpc``); False if the bulk (95th percentile) is beyond it."""
    return bool(np.quantile(post.d_mpc, 0.95) <= d_ref_mpc)


def distance_overlap(gw: dict, ra: float, dec: float, post: DistanceInclinationPosterior,
                     d_ref_mpc: Optional[float] = None) -> float:
    """I_D = E_{D ~ p_KN}[ p_GW(D | Omega_KN) / pi(D) ],  pi(D) = 3 D^2 / D_ref^3.

    > 1 when the GW distance along the kilonova's direction agrees with the
    kilonova distance better than a random distance would. ``d_ref_mpc`` is
    the reach of the GW search (distance distribution of unrelated GW
    candidates); use the network horizon, not an arbitrary prior bound,
    or the factor is inflated as (d_ref / D)^3. It must exceed the bulk of
    the kilonova posterior (check with :func:`distance_overlap_valid`); if it
    does not, the GW search could hardly have seen the source at the
    kilonova distance, which is itself evidence against the association.
    """
    d_ref = d_ref_mpc if d_ref_mpc is not None else post.d_max_prior
    if "DISTMU" not in gw:
        raise KeyError("GW sky map has no distance layers")
    nside = hp.npix2nside(len(gw["PROB"]))
    ipix = hp.ang2pix(nside, ra, dec, lonlat=True, nest=gw.get("nested", True))
    mu, sig, nrm = gw["DISTMU"][ipix], gw["DISTSIGMA"][ipix], gw["DISTNORM"][ipix]
    if not np.isfinite(mu) or sig <= 0:
        return 0.0
    d = post.d_mpc
    # p_GW / pi = nrm d^2 N(d) / (3 d^2 / Dmax^3)
    return float(np.mean(nrm * norm.pdf(d, mu, sig) * d_ref ** 3 / 3))
