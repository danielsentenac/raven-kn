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


def synthetic_fetcher(inject: bool, distance_mpc: float = 120.0, seed: int = 1, flip=()):
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
        s = seed + {"H1": 0, "L1": 7}.get(ifo, 13)
        ts = pycbc.noise.noise_from_psd(n, 1.0 / FS, psd, seed=s)
        ts.start_time = start
        if inject:
            m = 1.2 * 2 ** 0.2
            hp, hc = get_td_waveform(approximant="TaylorT4", mass1=m, mass2=m, delta_t=1.0 / FS,
                                     f_lower=35.0, distance=distance_mpc)
            hp.start_time += T_INJ
            hc.start_time += T_INJ
            sig = Detector(ifo).project_wave(hp, hc, np.radians(RA), np.radians(DEC), PSI)
            if ifo in flip:                    # not a signal from that direction
                sig *= -1
            i0 = int(round((float(sig.start_time) - start) * FS))
            a, h = ts.numpy(), sig.numpy()
            lo, hi = max(i0, 0), min(i0 + len(h), n)      # overlap only (no wrap-around)
            if hi > lo:
                a[lo:hi] += h[lo - i0:hi - i0]
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
    assess(make_kn(), tmp_path, chirp_mass_msun=1.2, require_kn_evidence=False)
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


def _coh_cfg(**kw):
    return TargetedSearchConfig(**{**_cfg().__dict__, "ifos": ("H1", "L1", "V1"),
                                   "statistic": "coherent", **kw})


def test_coherent_projector_and_null_cut():
    from raven_kn.targeted.search import coherent_projector, null_cut
    P2 = coherent_projector([1.0, 2.0], [0.3, -0.5], [0.6, 0.2])
    assert np.allclose(P2, np.eye(2))                   # 2 detectors: no null stream
    P3 = coherent_projector([1.0, 1.0, 0.5], [0.3, -0.5, 0.1], [0.6, 0.2, -0.7])
    assert np.allclose(P3 @ P3, P3) and np.trace(P3) == pytest.approx(2)
    assert list(null_cut(np.array([10.0, 10.0, 30.0]), np.array([5.0, 6.0, 7.0]))) == \
        [False, True, False]                            # 30: threshold 5.25 + 0.2 * 10


def test_coherent_recovers_injection():
    s = run(_coh_cfg(), synthetic_fetcher(inject=True), progress=False).subsegments[0]
    assert s.statistic == "coherent" and abs(s.loudest_gps - T_INJ) < 0.01
    assert abs(s.loudest_mchirp - 1.2) < 0.005
    assert s.loudest_rho_coh > 10 and s.loudest_rho_null < 3
    assert s.loudest_rho_coh == pytest.approx(s.loudest_rho_network, rel=0.05)
    assert s.p_value == pytest.approx(1 / (len(s.background) + 1))


def test_coherent_null_stream_rejects_inconsistent_phases():
    """V1's signal sign-flipped: loud in every detector, but no source at that
    sky position gives those relative phases, so the null stream vetoes it."""
    ok = run(_coh_cfg(), synthetic_fetcher(inject=True, distance_mpc=60), progress=False).subsegments[0]
    bad = run(_coh_cfg(), synthetic_fetcher(inject=True, distance_mpc=60, flip=("V1",)),
              progress=False).subsegments[0]
    assert ok.loudest_rho_coh > 20 and ok.loudest_rho_null < 4
    assert "failed the null cut" in bad.notes
    assert bad.loudest_rho_net < 0.6 * ok.loudest_rho_net


def test_coherent_noise_only():
    s = run(_coh_cfg(), synthetic_fetcher(inject=False), progress=False).subsegments[0]
    assert s.loudest_rho_net < 9 and s.p_value > 0.05


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


def test_offsource_trials():
    from raven_kn.targeted.search import offsource_windows
    cfg = TargetedSearchConfig(**{**_cfg().__dict__, "n_offsource": 4, "n_slides": 20})
    w = offsource_windows(cfg, 64, 4)
    assert w[0][1] == cfg.gps_start - 64 and w[1][0] == cfg.gps_end + 64
    assert all(b <= cfg.gps_start or a >= cfg.gps_end for a, b in w)
    s = run(cfg, synthetic_fetcher(inject=True), progress=False).subsegments[0]
    assert len(s.offsource) == 4 and max(s.offsource) < 9
    assert s.p_value_offsource == pytest.approx(1 / 5)


def test_short_window_uses_offsource_pvalue():
    from raven_kn.targeted.search import SubsegmentResult, subsegment_pvalue
    r = SubsegmentResult(0.0, 6.0, ["H1", "L1"], 9.0, 3.0, 1.2, {}, [5.0], 0.5, {}, "")
    assert subsegment_pvalue(r) == 0.5                       # no off-source trials: slides
    r.p_value_offsource = 0.01
    assert subsegment_pvalue(r) == 0.01                      # 1 slide only: off-source alone
    r.background = [5.0] * 50
    assert subsegment_pvalue(r) == 0.5                       # enough slides: the larger


def _spiked(fetch, ifo="L1", dt=-1.1, amp=5e-19):
    """A loud short glitch in one detector ``dt`` s from the merger, as in
    GW170817's L1 data: the search gates it, and the gate cuts the chirp."""
    def f(i, start, end):
        ts = fetch(i, start, end).copy()
        if i == ifo:
            k = int(round((T_INJ + dt - float(ts.start_time)) * FS))
            if 0 <= k < len(ts):
                ts.data[k - 2:k + 3] += amp * np.array([0.3, 0.7, 1.0, 0.7, 0.3])
        return ts
    return f


def test_gate_aware_chisq_keeps_gated_signal():
    """The injection overlaps a gate in L1. With the template gated like the
    data, L1's chi-square stays near 1 and newSNR keeps most of the gated SNR;
    the plain chi-square sees the zeroed stretch of the chirp as a power
    deficit and cuts L1 hard (GW170817: 24 -> 16)."""
    base = {**_cfg().__dict__, "gate_threshold": 50.0}
    fetch = _spiked(synthetic_fetcher(inject=True, distance_mpc=60))
    aware = run(TargetedSearchConfig(**base), fetch, progress=False).subsegments[0]
    plain = run(TargetedSearchConfig(**{**base, "gate_aware_chisq": False}), fetch,
                progress=False).subsegments[0]
    raw = run(TargetedSearchConfig(**{**base, "chisq_bins": None}), fetch, progress=False).subsegments[0]
    assert aware.gates["L1"] and abs(aware.gates["L1"][0] - (T_INJ - 1.1)) < 0.1
    assert abs(aware.loudest_gps - T_INJ) < 0.01
    l1_raw, l1_aware, l1_plain = (r.loudest_rho_ifo["L1"] for r in (raw, aware, plain))
    assert l1_aware > 0.85 * l1_raw
    assert l1_plain < 0.8 * l1_aware


def test_gate_aware_chisq_is_standard_away_from_gates():
    """No gate in the template span: identical to the standard chi-square."""
    a = run(TargetedSearchConfig(**{**_cfg().__dict__}), synthetic_fetcher(inject=True),
            progress=False).subsegments[0]
    b = run(TargetedSearchConfig(**{**_cfg().__dict__, "gate_aware_chisq": False}),
            synthetic_fetcher(inject=True), progress=False).subsegments[0]
    assert a.loudest_rho_net == pytest.approx(b.loudest_rho_net)


def test_coherent_network_chisq_keeps_gated_signal():
    """Coherent statistic, injection overlapping an L1 gate: the coherent SNR
    reweighted by the network chi-square (PyGRB) stays close to the raw one."""
    base = {**_coh_cfg().__dict__, "gate_threshold": 50.0}
    fetch = _spiked(synthetic_fetcher(inject=True, distance_mpc=60))
    rw = run(TargetedSearchConfig(**base), fetch, progress=False).subsegments[0]
    raw = run(TargetedSearchConfig(**{**base, "chisq_bins": None}), fetch, progress=False).subsegments[0]
    assert rw.gates["L1"] and abs(rw.loudest_gps - T_INJ) < 0.01
    assert "network chi-square" in rw.notes and rw.loudest_chisq_r is not None
    assert rw.loudest_rho_net > 0.85 * raw.loudest_rho_net


def test_coherent_network_chisq_suppresses_glitch():
    """A loud H1-only sine-Gaussian (no gating, no coincidence requirement):
    the network chi-square brings the coherent statistic down."""
    base = {**_coh_cfg().__dict__, "coinc_threshold": 0.0}
    raw = run(TargetedSearchConfig(**{**base, "chisq_bins": None}), _glitchy(), progress=False).subsegments[0]
    rw = run(TargetedSearchConfig(**base), _glitchy(), progress=False).subsegments[0]
    assert raw.loudest_rho_net > 10
    assert rw.loudest_rho_net < 0.7 * raw.loudest_rho_net
