import json

import numpy as np
import pytest

from raven_kn.connection import (FinkConnector, GraceDbConnector, KNCandidate,
                                 PhotometryPoint, Redshift, ViewingAngleMeasurement, gps_to_mjd)
from raven_kn.search import assess, calc_signif_gracedb, coinc_far, search
from raven_kn.skymap import (ansatz_from_moments, distance_overlap, kn_skymap_3d,
                             read_gw_skymap, skymap_overlap_integral, write_skymap)
from raven_kn.treatment import (combine_subsegments, posterior_for, targeted_far,
                                targeted_pvalue, window_from_bracket, window_from_grb)
from raven_kn.treatment.far import coinc_far_kn

T0 = 1187008882.43
H = 1 / 24


def make_kn(ul_hours_before=3.0, ul_mag=24.5, grb=False, z=0.009783, theta=20.0):
    t0 = gps_to_mjd(T0)
    phot = [PhotometryPoint(t0 - ul_hours_before * H, "g", ul_mag, upper_limit=True, survey="LSST"),
            PhotometryPoint(t0 + 10.9 * H, "i", 17.1, 0.03, survey="LSST"),
            PhotometryPoint(t0 + 1.4, "r", 17.6, 0.03, survey="LSST")]
    kn = KNCandidate("test", 197.450374, -23.381495, photometry=phot,
                     redshift=Redshift(z, 1e-5) if z else None,
                     viewing_angle=ViewingAngleMeasurement(theta, 5.0) if theta else None)
    if grb:
        from raven_kn.connection import GRBAssociation
        kn.grb = GRBAssociation(T0 + 1.74, "Fermi-GBM")
    return kn


# ---------------- t0 window ----------------

def test_bracket_contains_merger_and_is_short_with_deep_limit():
    w = window_from_bracket(make_kn(), d_mpc_upper=50)
    assert w.gps_start < T0 < w.gps_end
    assert w.duration < 1.0 * 86400


def test_shallow_limit_is_not_constraining():
    w = window_from_bracket(make_kn(ul_mag=17.0), d_mpc_upper=50)
    assert w.duration > 6 * 86400


def test_grb_window_is_seconds():
    w = window_from_grb(T0 + 1.74)
    assert w.gps_start < T0 < w.gps_end and w.duration == 6.0


# ---------------- distance / inclination ----------------

def test_posterior_recovers_ngc4993_distance_and_angle():
    post = posterior_for(make_kn())
    s = post.summary()
    assert 35 < s["luminosity_distance_mpc"]["median"] < 50
    assert 12 < s["viewing_angle_deg"]["median"] < 28
    pri = post.search_priors()
    assert pri["luminosity_distance_mpc"][0] < 40 < pri["luminosity_distance_mpc"][1]


# ---------------- decision ----------------

def test_nearby_is_worth_far_is_not(tmp_path):
    near = assess(make_kn(), tmp_path)
    assert near["trigger_targeted_search"]
    assert (tmp_path / "test_targeted_search_request.json").exists()
    far = assess(make_kn(z=0.15, theta=None))
    assert not far["trigger_targeted_search"]
    assert any("Mpc cut" in r for r in far["detectability"]["reasons"])


def test_detector_off_means_no_search():
    res = assess(make_kn(), segments={"H1": [], "L1": [], "V1": []})
    assert not res["trigger_targeted_search"]


# ---------------- sky maps & FAR ----------------

def test_ansatz_moments():
    mu, sig, nrm = ansatz_from_moments(40.0, 3.0)
    r = np.linspace(1e-3, 80, 20000)
    p = nrm * r ** 2 * np.exp(-0.5 * ((r - mu) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
    assert np.trapezoid(p, r) == pytest.approx(1, rel=1e-3)
    assert np.trapezoid(r * p, r) == pytest.approx(40, rel=1e-3)


def test_skymap_roundtrip_and_overlaps(tmp_path):
    kn = make_kn()
    post = posterior_for(kn)
    cols = kn_skymap_3d(kn.ra, kn.dec, post, nside=64)
    f = write_skymap(tmp_path / "kn.fits", cols, "test", T0, post)
    gw = read_gw_skymap(f)
    assert skymap_overlap_integral(gw["PROB"], kn.ra, kn.dec) == pytest.approx(len(gw["PROB"]))
    assert distance_overlap(gw, kn.ra, kn.dec, post) > 1
    # a GW map pointing to 400 Mpc disfavours the association
    gw_far = dict(gw, DISTMU=np.where(np.isfinite(gw["DISTMU"]), 400.0, np.inf))
    assert distance_overlap(gw_far, kn.ra, kn.dec, post) < 1e-3


def test_coinc_far_scales_with_window():
    a = coinc_far_kn(86400, 1e-7)["temporal_coinc_far"]
    b = coinc_far_kn(3600, 1e-7)["temporal_coinc_far"]
    assert a / b == pytest.approx(24)


def test_targeted_stats():
    assert targeted_pvalue(10, [1, 2, 3]) == pytest.approx(0.25)
    assert combine_subsegments([0.01, 0.5]) == pytest.approx(1 - 0.99 ** 2)
    assert targeted_far(1e4, 1e-6, 1e-6) > 0


# ---------------- connectors (mocked) ----------------

class MockGraceDb:
    def __init__(self, superevents):
        self._se, self.logs = superevents, []

    def superevents(self, query):
        lo, hi = (float(x) for x in query.split(" .. "))
        return [s for s in self._se if lo <= s["t_0"] <= hi]

    def writeLog(self, graceid, message, filename=None, filecontents=None, tag_name=()):
        self.logs.append((graceid, message, filename, filecontents))


def test_search_and_signif_with_mock_gracedb():
    mock = MockGraceDb([{"superevent_id": "S170817a", "t_0": T0, "far": 1e-7},
                        {"superevent_id": "S170820x", "t_0": T0 + 3 * 86400, "far": 1e-5}])
    gdb = GraceDbConnector(client=mock)
    kn = make_kn()
    win = window_from_bracket(kn, d_mpc_upper=50)
    found = search(kn, win, gdb)
    assert [s["superevent_id"] for s in found] == ["S170817a"]
    res = calc_signif_gracedb(kn, found[0], win, gdb)
    assert res["temporal_coinc_far"] < 1e-7
    assert any(f == "coincidence_far_kn.json" for _, _, f, _ in mock.logs)
    assert json.loads(mock.logs[-1][3])["superevent_id"] == "S170817a"


def test_fink_parser():
    rows = [{"i:jd": 2460000.5, "i:fid": 1, "d:tag": "upperlim", "i:diffmaglim": 20.5,
             "i:ra": 1.0, "i:dec": 2.0},
            {"i:jd": 2460001.5, "i:fid": 2, "d:tag": "valid", "i:magpsf": 19.0,
             "i:sigmapsf": 0.1, "i:ra": 1.0, "i:dec": 2.0},
            {"i:jd": 2460001.6, "i:fid": 2, "d:tag": "badquality", "i:magpsf": 19.0,
             "i:sigmapsf": 0.1, "i:ra": 1.0, "i:dec": 2.0}]

    class Sess:
        def post(self, url, json, timeout):
            class R:
                def raise_for_status(self): pass
                def json(self): return rows
            return R()
    kn = FinkConnector(session=Sess()).candidate("ZTF00test")
    assert len(kn.photometry) == 2 and kn.upper_limits[0].mag == 20.5
    assert kn.first_detection.band == "r"


def test_raven_event_dict_shape():
    d = make_kn(grb=True).to_raven_event_dict()
    assert d["extra_attributes"]["GRB"]["ra"] == pytest.approx(197.450374)
    assert d["gpstime"] == pytest.approx(T0 + 1.74)


# ---------------- model-systematic margin ----------------

def _biased_fit(kn, shift_h, sigma_h):
    rng = np.random.default_rng(0)
    kn.t0_samples_gps = list(T0 + shift_h * 3600 + rng.normal(0, sigma_h * 3600, 4000))
    return kn


def test_model_margin_recovers_biased_fit_without_limit():
    from raven_kn.treatment.t0 import estimate_window
    # case-A-like: no pre-merger limit, fit 41.5 h early, 90 % width ~38 h
    kn = _biased_fit(make_kn(ul_mag=17.0), -41.5, 11.7)
    stat = estimate_window(kn, d_mpc_upper=50, model_sys_frac=0.0)
    sys_ = estimate_window(kn, d_mpc_upper=50, model_sys_frac=0.5)
    assert not (stat.gps_start <= T0 <= stat.gps_end)
    assert sys_.gps_start <= T0 <= sys_.gps_end
    assert sys_.gps_end <= kn.first_detection.gps          # still capped by the data
    assert "model-systematic margin" in sys_.notes


def test_model_margin_small_when_anchored_by_limit():
    from raven_kn.treatment.t0 import estimate_window
    kn = _biased_fit(make_kn(), -1.9, 2.6)                  # case-B-like
    w = estimate_window(kn, d_mpc_upper=50, model_sys_frac=0.5)
    early = (kn.t0_samples_gps and np.quantile(kn.t0_samples_gps, 0.05)) - w.gps_start
    assert early <= 1.01 * 3600                             # early side: only the 1 h floor
    assert w.duration < 24 * 3600 and w.gps_start <= T0 <= w.gps_end


# ---------------- chirp-mass scaling ----------------

def test_scale_ranges_chirp_mass():
    from raven_kn.treatment.detectability import BNS_CHIRP_MASS, scale_ranges
    r = {"H1": 150.0, "L1": 160.0}
    assert scale_ranges(r, None) == r
    assert scale_ranges(r, BNS_CHIRP_MASS)["L1"] == pytest.approx(160.0)
    assert scale_ranges(r, 0.5)["L1"] == pytest.approx(160.0 * (0.5 / 1.2188) ** (5 / 6))
    with pytest.raises(ValueError):
        scale_ranges(r, 0.0)


def test_subsolar_chirp_mass_lowers_detectability():
    kn = make_kn(z=0.0848, theta=None)          # ~400 Mpc, S250818k-like host
    bns = assess(kn)
    ssm = assess(kn, chirp_mass_msun=0.5)
    assert ssm["detectability"]["snr_median"] < 0.6 * bns["detectability"]["snr_median"]
    assert not ssm["trigger_targeted_search"]
    assert any("chirp mass 0.5" in r for r in ssm["detectability"]["reasons"])


def test_conditional_distance_is_not_ansatz_parameters():
    from raven_kn.skymap import conditional_distance, distance_overlap_valid
    nside = 16
    npix = 12 * nside ** 2
    gw = {"PROB": np.full(npix, 1.0 / npix), "DISTMU": np.full(npix, 218.0),
          "DISTSIGMA": np.full(npix, 80.0), "DISTNORM": np.ones(npix), "nested": True}
    c = conditional_distance(gw, 10.0, 20.0, d_compare_mpc=400.0)
    assert c["distmu"] == 218.0
    assert 260 < c["mean_mpc"] < 280            # r^2 factor shifts the mean above DISTMU
    assert c["interval90_mpc"][1] < 400 and 0.02 < c["p_beyond"] < 0.06
    assert c["likelihood_ratio_to_peak"] == pytest.approx(np.exp(-0.5 * (182 / 80) ** 2))
    post = posterior_for(make_kn(z=0.0848, theta=None))      # ~400 Mpc
    assert not distance_overlap_valid(post, 362.0) and distance_overlap_valid(post, 500.0)
