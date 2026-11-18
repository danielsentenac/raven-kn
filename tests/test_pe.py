"""EM priors for the parameter estimation of a search candidate, and the
relative-binning likelihood on a synthetic face-on signal."""
import numpy as np
import pytest

bilby = pytest.importorskip("bilby")
pytest.importorskip("pycbc")

from raven_kn.targeted.pe import _interferometers, em_priors, summarise  # noqa: E402
from raven_kn.targeted.search import TargetedSearchConfig  # noqa: E402
from test_targeted import DEC, RA, T_INJ, synthetic_fetcher  # noqa: E402


def _cfg(**kw):
    return TargetedSearchConfig(RA, DEC, T_INJ - 3, T_INJ + 3, ifos=("H1", "L1"), f_low=40.0,
                                distance_range_mpc=(90.0, 150.0),
                                inclination_ranges_rad=[(0.0, 0.5), (np.pi - 0.5, np.pi)], **kw)


def test_em_priors(tmp_path):
    p = em_priors(_cfg(), 1.2, T_INJ)
    s = p.sample(2000)
    assert s["luminosity_distance"].min() >= 90 and s["luminosity_distance"].max() <= 150
    t = s["theta_jn"]
    assert np.mean((t < 0.52) | (t > np.pi - 0.52)) > 0.97          # only the two allowed ranges
    assert 0.3 < np.mean(t < 1.5) < 0.7                              # both orientations
    assert p["ra"].peak == pytest.approx(np.radians(RA)) and abs(s["chirp_mass"] - 1.2).max() < 0.025
    assert abs(s["geocent_time"] - T_INJ).max() <= 0.1
    gw = em_priors(_cfg(), 1.2, T_INJ, em=False).sample(2000)
    assert gw["luminosity_distance"].max() > 400 and np.mean(np.abs(gw["theta_jn"] - np.pi / 2) < 0.5) > 0.3
    np.savez(tmp_path / "s.npz", luminosity_distance=np.random.default_rng(0).normal(120, 5, 5000),
             cos_iota=np.random.default_rng(1).uniform(0.9, 1.0, 5000))
    j = em_priors(_cfg(prior_samples=str(tmp_path / "s.npz")), 1.2, T_INJ).sample(2000)
    assert j["luminosity_distance"].std() < 7 and j["theta_jn"].max() < 0.52
    with pytest.raises(ValueError):
        em_priors(TargetedSearchConfig(RA, DEC, 0, 1), 1.2, T_INJ)


def test_likelihood_peaks_at_the_injection():
    """Injection at 120 Mpc, face-on: the likelihood is higher at the true
    distance than 3 times closer or farther, and at face-on than edge-on."""
    cfg = _cfg()
    ifos = _interferometers(cfg, synthetic_fetcher(inject=True), ("H1", "L1"), T_INJ, 32.0, 800.0, 256.0)
    assert len(ifos) == 2 and ifos[0].strain_data.duration == 32.0
    wg = bilby.gw.WaveformGenerator(
        duration=32.0, sampling_frequency=2048,
        frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        waveform_arguments=dict(waveform_approximant="IMRPhenomD", reference_frequency=40.0,
                                minimum_frequency=40.0))
    like = bilby.gw.likelihood.GravitationalWaveTransient(ifos, wg)
    base = dict(chirp_mass=1.2, mass_ratio=1.0, chi_1=0.0, chi_2=0.0, ra=np.radians(RA), dec=np.radians(DEC),
                luminosity_distance=120.0, theta_jn=0.0)

    def best(**kw):                 # maximum over the unknown phase, polarisation and time offsets
        out = -np.inf
        for t in T_INJ + np.arange(-3, 4) * 1e-3:
            for ph in np.linspace(0, np.pi, 8, endpoint=False):
                like.parameters.update({**base, "psi": 0.3, "phase": ph, "geocent_time": t, **kw})
                out = max(out, like.log_likelihood_ratio())
        return out
    at = best()
    assert at > 50
    assert best(luminosity_distance=40.0) < at and best(luminosity_distance=360.0) < at
    assert best(theta_jn=np.pi / 2) < at


def test_summarise():
    post = {"chirp_mass": np.full(100, 1.2), "mass_ratio": np.full(100, 0.8), "chi_1": np.zeros(100),
            "chi_2": np.zeros(100), "luminosity_distance": np.linspace(30, 50, 100),
            "theta_jn": np.full(100, 2.6)}
    s = summarise(post, lambda d: d * 2.3e-4)
    assert s["viewing_angle_deg"]["median"] == pytest.approx(np.degrees(np.pi - 2.6))
    assert s["chirp_mass_source"]["median"] == pytest.approx(1.2 / (1 + 40 * 2.3e-4), rel=1e-3)
    assert s["chi_eff"]["median"] == 0.0


def test_maximise_recovers_the_injection():
    """Started 0.03 % off in chirp mass and 4 ms off in time, the exact
    maximisation recovers the injection (120 Mpc, face-on) and gives the
    relative-binning likelihood a reference point with the same value."""
    from raven_kn.targeted.pe import maximise
    cfg = _cfg()
    ifos = _interferometers(cfg, synthetic_fetcher(inject=True), ("H1", "L1"), T_INJ, 32.0, 800.0, 256.0)
    ml = maximise(cfg, ifos, 1.2 * 1.0003, T_INJ + 0.004, 32.0, maxiter=150)
    p = ml["parameters"]
    assert ml["snr"] > 10 and ml["log_likelihood_ratio_at_phase"] == pytest.approx(ml["log_likelihood_ratio"], abs=1.0)
    assert p["chirp_mass"] == pytest.approx(1.2, rel=5e-4) and abs(p["geocent_time"] - T_INJ) < 2e-3
    assert 60 < p["luminosity_distance"] < 200
    wg = bilby.gw.WaveformGenerator(
        duration=32.0, sampling_frequency=2048,
        frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole_relative_binning,
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        waveform_arguments=dict(waveform_approximant="IMRPhenomD", reference_frequency=40.0,
                                minimum_frequency=40.0))
    like = bilby.gw.likelihood.RelativeBinningGravitationalWaveTransient(
        interferometers=ifos, waveform_generator=wg, fiducial_parameters=dict(p),
        update_fiducial_parameters=False)
    like.parameters.update(p)
    assert like.log_likelihood_ratio() == pytest.approx(ml["log_likelihood_ratio_at_phase"], abs=0.5)
    like.parameters.update({**p, "luminosity_distance": p["luminosity_distance"] / 3})
    assert like.log_likelihood_ratio() < ml["log_likelihood_ratio_at_phase"]
