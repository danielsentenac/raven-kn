"""Parameter estimation of a targeted-search candidate with the EM priors.

The light search gives a time and a chirp mass at the known sky position.
This module runs Bilby on the same data with

- the sky position fixed at the EM position;
- the luminosity distance and the inclination (theta_jn) restricted by the
  EM measurement: the marginals of the joint distance-inclination samples of
  the search request when present, else the request's intervals (distance
  uniform in volume, isotropic orientation inside the allowed ranges);
- narrow chirp-mass and time priors around the search candidate.

``em=False`` keeps the sky position but takes a distance uniform in volume up
to ``d_max_mpc`` and an isotropic orientation: comparing the two posteriors
shows what the EM distance and inclination bring (the distance-inclination
degeneracy of the GW amplitude is broken, hence the inclination measured).

Aligned-spin IMRPhenomD without tides by default, relative-binning
likelihood (fast for long signals), dynesty. A follow-up of a candidate, not
a replacement of the LVK parameter-estimation pipelines.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from .search import StrainFetcher, TargetedSearchConfig, _masses


def _interped(x: np.ndarray, lo: float, hi: float, name: str, label: str, bins: int = 80, weight=None):
    """Bilby prior following the histogram of samples ``x`` within [lo, hi]."""
    from bilby.core.prior import Interped
    h, e = np.histogram(x, bins=bins, range=(lo, hi), weights=weight)
    xx = 0.5 * (e[1:] + e[:-1])
    return Interped(np.concatenate([[lo], xx, [hi]]), np.concatenate([[h[0]], h, [h[-1]]]).astype(float),
                    minimum=lo, maximum=hi, name=name, latex_label=label)


def em_priors(cfg: TargetedSearchConfig, mchirp: float, gps: float, em: bool = True,
              mchirp_frac: float = 0.02, q_min: float = 0.125, chi_max: float = 0.05,
              dt: float = 0.1, d_max_mpc: Optional[float] = None):
    """Bilby priors for a candidate of the search configured by ``cfg`` (which
    carries the EM priors of the request)."""
    import bilby
    from bilby.core.prior import Constraint, PowerLaw, Sine, Uniform
    from bilby.gw.prior import AlignedSpin, UniformInComponentsChirpMass, UniformInComponentsMassRatio

    p = bilby.gw.prior.BBHPriorDict(dictionary={})
    for k in list(p.keys()):
        del p[k]
    p["chirp_mass"] = UniformInComponentsChirpMass(mchirp * (1 - mchirp_frac), mchirp * (1 + mchirp_frac),
                                                   name="chirp_mass")
    p["mass_ratio"] = UniformInComponentsMassRatio(q_min, 1.0, name="mass_ratio")
    p["mass_1"] = Constraint(0.1, 100.0, name="mass_1")
    p["mass_2"] = Constraint(0.1, 100.0, name="mass_2")
    for k in ("chi_1", "chi_2"):
        p[k] = AlignedSpin(a_prior=Uniform(0, chi_max), z_prior=Uniform(-1, 1), name=k)
    p["ra"], p["dec"] = math.radians(cfg.ra_deg), math.radians(cfg.dec_deg)
    p["psi"] = Uniform(0, math.pi, name="psi", boundary="periodic")
    p["phase"] = Uniform(0, 2 * math.pi, name="phase", boundary="periodic")
    p["geocent_time"] = Uniform(gps - dt, gps + dt, name="geocent_time")
    if not em:
        hi = d_max_mpc or (5 * cfg.distance_range_mpc[1] if cfg.distance_range_mpc else 500.0)
        p["luminosity_distance"] = PowerLaw(2, 1.0, hi, name="luminosity_distance", unit="Mpc")
        p["theta_jn"] = Sine(name="theta_jn")
        return p
    if not cfg.distance_range_mpc:
        raise ValueError("EM priors need distance_range_mpc (a search request)")
    lo, hi = cfg.distance_range_mpc
    ranges = cfg.inclination_ranges_rad or [(0.0, math.pi)]
    t_lo, t_hi = min(a for a, b in ranges), max(b for a, b in ranges)
    if cfg.prior_samples:
        f = np.load(cfg.prior_samples)
        p["luminosity_distance"] = _interped(np.asarray(f["luminosity_distance"]), lo, hi,
                                             "luminosity_distance", "$d_L$")
        p["theta_jn"] = _interped(np.arccos(np.clip(f["cos_iota"], -1, 1)), t_lo, t_hi,
                                  "theta_jn", r"$\theta_{JN}$")
    else:
        p["luminosity_distance"] = PowerLaw(2, lo, hi, name="luminosity_distance", unit="Mpc")
        t = np.linspace(t_lo, t_hi, 4001)
        inside = np.any([(t >= a) & (t <= b) for a, b in ranges], axis=0)
        p["theta_jn"] = _interped(t[inside], t_lo, t_hi, "theta_jn", r"$\theta_{JN}$",
                                  bins=200, weight=np.sin(t[inside]) + 1e-6)
    return p


def _interferometers(cfg, fetch: StrainFetcher, ifos: Sequence[str], gps: float, duration: float,
                     f_max: float, psd_duration: float):
    """Bilby interferometers with the strain of [gps + 2 - duration, gps + 2]
    and a median Welch PSD of the ``psd_duration`` s before it."""
    import bilby
    from gwpy.timeseries import TimeSeries
    end = math.floor(gps) + 2.0
    start = end - duration
    out = bilby.gw.detector.InterferometerList([])
    for name in ifos:
        ts = fetch(name, start - psd_duration, end)
        g = TimeSeries(ts.numpy().astype(float), t0=float(ts.start_time), dt=float(ts.delta_t))
        ifo = bilby.gw.detector.get_empty_interferometer(name)
        ifo.set_strain_data_from_gwpy_timeseries(g.crop(start, end))
        psd = g.crop(start - psd_duration, start).psd(fftlength=duration, overlap=0.0,
                                                      window=("tukey", 0.4 / duration), method="median")
        ifo.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
            frequency_array=psd.frequencies.value, psd_array=psd.value)
        ifo.minimum_frequency, ifo.maximum_frequency = cfg.f_low, f_max
        out.append(ifo)
    return out


def summarise(posterior, cosmology_z=None) -> dict:
    """Median and 90 % interval of the main parameters of a Bilby posterior."""
    def q(x):
        lo, med, hi = np.quantile(np.asarray(x, dtype=float), [0.05, 0.5, 0.95])
        return {"median": float(med), "lower": float(lo), "upper": float(hi)}
    out = {k: q(posterior[k]) for k in ("chirp_mass", "mass_ratio", "luminosity_distance", "theta_jn",
                                        "chi_1", "chi_2") if k in posterior}
    if "theta_jn" in posterior:
        out["viewing_angle_deg"] = q(np.degrees(np.minimum(posterior["theta_jn"],
                                                           math.pi - np.asarray(posterior["theta_jn"]))))
    if cosmology_z is not None and "luminosity_distance" in posterior:
        z = cosmology_z(np.asarray(posterior["luminosity_distance"]))
        out["chirp_mass_source"] = q(np.asarray(posterior["chirp_mass"]) / (1 + z))
    if {"chi_1", "chi_2", "mass_ratio"} <= set(posterior):
        qq = np.asarray(posterior["mass_ratio"])
        out["chi_eff"] = q((np.asarray(posterior["chi_1"]) + qq * np.asarray(posterior["chi_2"])) / (1 + qq))
    return out


def _z_of_distance(d_mpc: np.ndarray) -> np.ndarray:
    from astropy import units as u
    from astropy.cosmology import Planck18, z_at_value
    grid = np.linspace(max(d_mpc.min() * 0.9, 1e-3), d_mpc.max() * 1.1, 200)
    zg = np.array([float(z_at_value(Planck18.luminosity_distance, g * u.Mpc, zmax=10.0)) for g in grid])
    return np.interp(d_mpc, grid, zg)


def run_pe(cfg: TargetedSearchConfig, fetch: StrainFetcher, mchirp: float, gps: float,
           out_dir: str | Path, em: bool = True, ifos: Optional[Sequence[str]] = None,
           duration: float = 64.0, psd_duration: float = 1024.0, approximant: str = "IMRPhenomD",
           nlive: int = 500, npool: int = 1, label: Optional[str] = None, sampler_kw: Optional[dict] = None,
           **prior_kw) -> dict:
    """Bilby posterior of the candidate (``mchirp``, geocentric ``gps``) of a
    search configured by ``cfg``; ``em`` selects the EM priors. Writes the
    Bilby result and ``<label>_pe_summary.json`` to ``out_dir``."""
    import bilby
    ifos = list(ifos or cfg.ifos)
    label = label or f"{cfg.label}_pe_{'em' if em else 'gw'}"
    priors = em_priors(cfg, mchirp, gps, em=em, **prior_kw)
    f_max = 0.875 * cfg.sample_rate / 2
    ifo_list = _interferometers(cfg, fetch, ifos, gps, duration, f_max, psd_duration)
    m1, m2 = _masses(mchirp, 1.0)
    d0 = (float(np.sqrt(np.prod(cfg.distance_range_mpc))) if cfg.distance_range_mpc else 100.0)
    fiducial = dict(chirp_mass=mchirp, mass_ratio=0.9, chi_1=0.0, chi_2=0.0,
                    ra=math.radians(cfg.ra_deg), dec=math.radians(cfg.dec_deg), psi=1.0, phase=1.0,
                    geocent_time=gps, luminosity_distance=d0, theta_jn=0.5)
    wg = bilby.gw.WaveformGenerator(
        duration=duration, sampling_frequency=cfg.sample_rate,
        frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole_relative_binning,
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        waveform_arguments=dict(waveform_approximant=approximant, reference_frequency=cfg.f_low,
                                minimum_frequency=cfg.f_low))
    like = bilby.gw.likelihood.RelativeBinningGravitationalWaveTransient(
        interferometers=ifo_list, waveform_generator=wg, priors=priors, fiducial_parameters=fiducial,
        update_fiducial_parameters=True)
    res = bilby.run_sampler(like, priors, sampler="dynesty", nlive=nlive, npool=npool, outdir=str(out_dir),
                            label=label, check_point_delta_t=1800, **(sampler_kw or {}))
    summary = {"label": label, "em_priors": em, "ifos": ifos, "candidate": {"mchirp": mchirp, "gps": gps},
               "approximant": approximant, "duration_s": duration, "f_low": cfg.f_low, "f_max": f_max,
               "log_evidence": float(res.log_evidence), "log_bayes_factor": float(res.log_bayes_factor),
               "n_samples": int(len(res.posterior)),
               "posterior": summarise(res.posterior, _z_of_distance)}
    Path(out_dir, f"{label}_pe_summary.json").write_text(json.dumps(summary, indent=2))
    return summary
