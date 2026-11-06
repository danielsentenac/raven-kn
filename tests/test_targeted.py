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


def _glitchy():
    base = synthetic_fetcher(inject=False, seed=3)

    def fetch(ifo, start, end):
        ts = base(ifo, start, end).copy()
        if ifo == "H1":
            t = np.arange(len(ts)) / FS + float(ts.start_time) - (T_INJ - 10)
            ts.data += 4e-21 * np.exp(-(t / 0.02) ** 2) * np.sin(2 * np.pi * 150 * t)
        return ts
    return fetch


def test_chisq_suppresses_glitch():
    """A loud sine-Gaussian in H1 only (gating off): chi-square down-weights it."""
    no_coinc = {**_cfg().__dict__, "coinc_threshold": 0.0}
    raw = run(TargetedSearchConfig(**{**no_coinc, "chisq_bins": None}), _glitchy(), progress=False)
    rw = run(TargetedSearchConfig(**no_coinc), _glitchy(), progress=False)
    a, b = raw.subsegments[0], rw.subsegments[0]
    # the template merger time lies ~0.7 s after its 150 Hz part, where the glitch is
    assert a.loudest_rho_net > 12 and abs(a.loudest_gps - (T_INJ - 10)) < 2
    assert b.loudest_rho_ifo["H1"] < 0.7 * a.loudest_rho_ifo["H1"]


def test_single_detector_glitch_not_a_candidate():
    """An H1-only glitch with marginal L1 noise may still be the loudest value,
    but it is not removed from the background (2nd detector too weak), so its
    slid copies keep it non-significant."""
    s = run(_cfg(), _glitchy(), progress=False).subsegments[0]
    assert "removed from background" not in s.notes
    assert s.p_value > 0.05 and s.p_value_extrapolated > 0.05


def test_chisq_keeps_real_signal():
    raw = run(TargetedSearchConfig(**{**_cfg().__dict__, "chisq_bins": None}),
              synthetic_fetcher(inject=True), progress=False).subsegments[0]
    rw = run(_cfg(), synthetic_fetcher(inject=True), progress=False).subsegments[0]
    # TaylorT4 signal vs TaylorF2 templates: a few % mismatch raises chi-square at SNR ~25
    # (excess ~ 2 x mismatch x SNR^2), so newSNR loses 10-20 %; GW170817 keeps L1 24.5/26.4
    assert rw.loudest_rho_net > 0.8 * raw.loudest_rho_net


def test_slide_shifts():
    from raven_kn.targeted.search import slide_shifts
    two = slide_shifts(2, 50, 1000)
    assert len(two) == 49 and all(s[1] != 0 for s in two)
    three = slide_shifts(3, 50, 1000)
    assert len(three) == 1000 and len(set(three)) == 1000
    assert all(s[1] != s[2] and s[1] and s[2] for s in three)


def test_pygrb_overrides_and_bank(tmp_path):
    import configparser
    from raven_kn.targeted.pygrb import write_pygrb_overrides
    from pycbc.waveform.bank import TemplateBank
    cfg = TargetedSearchConfig(197.45, -23.38, 1187008000.0, 1187010048.0, mchirp_min=1.19,
                               mchirp_max=1.21, f_low=30.0, subsegment_s=1024, label="GW170817")
    files = write_pygrb_overrides(cfg, tmp_path)
    assert len(files) == 2 and (tmp_path / "run_pygrb.sh").exists()
    cp = configparser.ConfigParser()
    cp.optionxform = str
    cp.read(files[0])
    w = cp["workflow"]
    assert w["trigger-time"] == "1187008512" and w["ra"] == "197.450000deg"
    e = cp["workflow-exttrig_segments"]
    assert int(e["on-before"]) + int(e["on-after"]) == 1024
    assert int(e["max-duration"]) > int(e["min-duration"]) > 1024
    assert cp["workflow-tmpltbank"]["tmpltbank-method"] == "PREGENERATED_BANK"
    bank = TemplateBank(cp["workflow-tmpltbank"]["tmpltbank-pregenerated-bank"])
    assert len(bank.table) > 3
    assert np.allclose(bank.table.mass1, bank.table.mass2)
