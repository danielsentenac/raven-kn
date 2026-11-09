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
    from pycbc.types import angle_as_radians
    from pycbc.waveform.bank import TemplateBank
    tpl = tmp_path / "GRB170817A.ini"
    tpl.write_text("[workflow]\ntrigger-name = 170817A\ntrigger-time = 1187008884\n"
                   "ra = 3.44 rad\ndec = -0.41 rad\nsky-error = 0.0 rad\n"
                   "[workflow-injections]\nnsns-config-files = nsns.ini, GRB170817A.ini\n"
                   "[prior-ra+dec]\nmean_ra = ${workflow|ra}\n")
    cfg = TargetedSearchConfig(197.45, -23.38, 1187008000.0, 1187010048.0, mchirp_min=1.19,
                               mchirp_max=1.21, f_low=30.0, subsegment_s=1024, label="GW170817")
    files = write_pygrb_overrides(cfg, tmp_path / "out", base_config=["analysis.ini", "data.ini"],
                                  event_template=tpl)
    assert len(files) == 2
    cp = configparser.ConfigParser(interpolation=None)
    cp.optionxform = str
    cp.read(files[0])
    w = cp["workflow"]
    assert w["trigger-time"] == "1187008512" and w["sky-error"] == "0.0 rad"
    assert angle_as_radians(w["ra"]) == pytest.approx(np.radians(197.45))
    assert angle_as_radians(w["dec"]) == pytest.approx(np.radians(-23.38))
    assert cp["workflow-injections"]["nsns-config-files"] == "nsns.ini, GW170817_seg000.ini"
    assert cp["prior-ra+dec"]["mean_ra"] == "${workflow|ra}"
    run = (tmp_path / "out" / "run_pygrb.sh").read_text()
    assert "--config-files analysis.ini data.ini GW170817_seg000.ini" in run
    assert "workflow-exttrig_segments:on-before:512" in run
    m = [x for x in run.split() if x.startswith("workflow-tmpltbank:tmpltbank-pregenerated-bank:")][0]
    bank = TemplateBank(m.split(":", 2)[2])
    assert len(bank.table) > 3 and np.allclose(bank.table.mass1, bank.table.mass2)


def test_veto_free_subsegments(tmp_path):
    from raven_kn.targeted.pygrb import veto_free_subsegments, write_pygrb_overrides
    v = [(1187009326, 1187009328), (1187009900, 1187009905)]
    segs = veto_free_subsegments(1187008000, 1187010048, 1024, v)
    assert all(b - a <= 1024 for a, b in segs)
    for a, b in segs:                      # no on-source overlaps a (guarded) veto
        assert all(b <= va - 1 or a >= vb + 1 for va, vb in v)
    assert sum(b - a for a, b in segs) == pytest.approx(2048 - 2 - 5 - 4)
    cfg = TargetedSearchConfig(197.45, -23.38, 1187008000.0, 1187010048.0, mchirp_min=1.19,
                               mchirp_max=1.21, f_low=30.0, subsegment_s=1024, label="GW170817")
    files = write_pygrb_overrides(cfg, tmp_path, bank=[1.19, 1.2, 1.21], vetoes=v)
    run = (tmp_path / "run_pygrb.sh").read_text()
    import configparser
    for f in files:
        cp = configparser.ConfigParser(interpolation=None)
        cp.read(f)
        t = int(cp["workflow"]["trigger-time"])
        name = cp["workflow"]["trigger-name"]
        half = int(run.split(f"--workflow-name {name}")[1].split("on-before:")[1].split()[0])
        assert all(t + half <= va or t - half >= vb for va, vb in v)
