"""Targeted search on synthetic Gaussian noise (aLIGO design PSD) with and
without an injected TaylorT4 signal at a known sky position."""
import numpy as np
import pytest

pycbc = pytest.importorskip("pycbc")

from raven_kn.targeted.search import (TargetedSearchConfig, build_bank,  # noqa: E402
                                      run, template_duration)

RA, DEC, PSI = 197.45, -23.38, 0.3
T_INJ = 1000000500.0          # geocentric merger time
FS = 2048


def synthetic_fetcher(inject: bool, distance_mpc: float = 120.0, seed: int = 1):
    import pycbc.noise
    import pycbc.psd
    from pycbc.detector import Detector
    from pycbc.waveform import get_td_waveform

    cache = {}

    def fetch(ifo, start, end):
        key = (ifo, start, end)
        if key in cache:
            return cache[key]
        n = int((end - start) * FS)
        df = 1.0 / 16
        psd = pycbc.psd.aLIGOZeroDetHighPower(int(FS / 2 / df) + 1, df, 10.0)
        s = seed + (0 if ifo == "H1" else 7)
        ts = pycbc.noise.noise_from_psd(n, 1.0 / FS, psd, seed=s)
        ts.start_time = start
        if inject:
            m = 1.2 * 2 ** 0.2
            hp, hc = get_td_waveform(approximant="TaylorT4", mass1=m, mass2=m, delta_t=1.0 / FS,
                                     f_lower=35.0, distance=distance_mpc)
            hp.start_time += T_INJ
            hc.start_time += T_INJ
            sig = Detector(ifo).project_wave(hp, hc, RA, DEC, PSI)
            i0 = int(round((float(sig.start_time) - start) * FS))
            a = ts.numpy()
            a[i0:i0 + len(sig)] += sig.numpy()[: max(0, min(len(sig), n - i0))]
        cache[key] = ts
        return ts
    return fetch


def _cfg():
    return TargetedSearchConfig(RA, DEC, T_INJ - 40, T_INJ + 24, ifos=("H1", "L1"),
                                mchirp_min=1.19, mchirp_max=1.21, f_low=40.0, sample_rate=FS,
                                subsegment_s=64, n_slides=60, slide_step_s=1.0, gate_threshold=None)


def test_bank_covers_range():
    b = build_bank(1.19, 1.21, 40.0)
    assert b[0] == 1.19 and b[-1] == pytest.approx(1.21) and 3 < len(b) < 200
    assert template_duration(1.2, 40.0) == pytest.approx(25, rel=0.15)


def test_injection_recovered():
    res = run(_cfg(), synthetic_fetcher(inject=True), progress=False)
    s = res.subsegments[0]
    assert abs(s.loudest_gps - T_INJ) < 0.05
    assert abs(s.loudest_mchirp - 1.2) < 0.005
    assert s.loudest_rho_net > 10
    assert s.p_value == pytest.approx(1 / (len(s.background) + 1))
    assert max(s.background) < s.loudest_rho_net


def test_noise_only_not_significant():
    res = run(_cfg(), synthetic_fetcher(inject=False), progress=False)
    s = res.subsegments[0]
    assert s.loudest_rho_net < 9
    assert s.p_value > 0.05


def test_config_from_assess_request(tmp_path):
    from raven_kn.search import assess
    from test_pipeline import make_kn
    assess(make_kn(), tmp_path, chirp_mass_msun=1.2)
    cfg = TargetedSearchConfig.from_request(tmp_path / "test_targeted_search_request.json")
    assert cfg.mchirp_min == pytest.approx(1.08) and cfg.mchirp_max == pytest.approx(1.32)
    assert cfg.gps_start < cfg.gps_end and cfg.ra_deg == pytest.approx(197.450374)
