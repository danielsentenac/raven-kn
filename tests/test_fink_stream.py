"""Fink livestream alerts -> candidates in the listener store."""

import numpy as np
import pytest
from astropy.time import Time

from raven_kn.connection import Redshift
from raven_kn.fink_stream import (alert_to_candidate, handle_alert, lsst_alert_to_candidate,
                                  stream, ztf_alert_to_candidate)
from raven_kn.listener import CandidateStore, handle_grb_notice

RA, DEC = 197.450374, -23.381495
T0_JD = Time("2017-08-17T12:41:04.4", scale="utc").jd       # GW170817


def ztf_alert(days, mags, fids, limit_days=(), score=0.9):
    """A ZTF-format Fink alert: last detection in 'candidate', the rest in
    'prv_candidates' (non-detections with magpsf None)."""
    rows = [{"jd": T0_JD + d, "fid": f, "magpsf": m, "sigmapsf": 0.05, "ra": RA, "dec": DEC,
             "diffmaglim": 21.0} for d, m, f in zip(days, mags, fids)]
    rows += [{"jd": T0_JD + d, "fid": 2, "magpsf": None, "sigmapsf": None, "diffmaglim": 20.5}
             for d in limit_days]
    rows.sort(key=lambda r: r["jd"])
    det = [r for r in rows if r["magpsf"] is not None]
    return {"objectId": "ZTF17aaaaaaa", "candidate": det[-1],
            "prv_candidates": [r for r in rows if r is not det[-1]], "rf_kn_vs_nonkn": score}


def test_ztf_alert_conversion():
    kn = ztf_alert_to_candidate(ztf_alert([0.5, 1.5], [17.3, 18.0], [1, 2], limit_days=[-1.0]))
    assert kn.name == "ZTF17aaaaaaa" and kn.kn_score == pytest.approx(0.9)
    assert len(kn.detections) == 2 and len(kn.upper_limits) == 1
    assert {p.band for p in kn.detections} == {"g", "r"}
    assert kn.detections[0].mjd == pytest.approx(T0_JD + 0.5 - 2400000.5)


def test_lsst_alert_conversion_dedups_forced():
    mjd_tai = Time(T0_JD + 0.6, format="jd").mjd + 37 / 86400
    src = {"diaObjectId": 42, "ra": RA, "dec": DEC, "midpointMjdTai": mjd_tai, "band": "i",
           "psfFlux": 10 ** ((31.4 - 17.1) / 2.5), "psfFluxErr": 1000.0, "visit": 7}
    forced_same = dict(src)                                       # forced phot. of the same visit
    forced_before = {"midpointMjdTai": mjd_tai - 1.0, "band": "r", "psfFlux": 50.0,
                     "psfFluxErr": 100.0, "visit": 6}             # S/N 0.5: a limit
    alert = {"diaObject": {"diaObjectId": 42, "ra": RA, "dec": DEC}, "diaSource": src,
             "prvDiaSources": [], "prvDiaForcedSources": [forced_same, forced_before]}
    kn = lsst_alert_to_candidate(alert)
    assert kn.name == "LSST42" and len(kn.photometry) == 2
    det = kn.detections[0]
    assert det.mag == pytest.approx(17.1, abs=1e-6) and det.band == "i"
    assert det.mjd == pytest.approx(Time(T0_JD + 0.6, format="jd").mjd, abs=1e-6)   # TAI -> UTC
    lim = kn.upper_limits[0]
    assert lim.mag == pytest.approx(31.4 - 2.5 * np.log10(500.0))


def test_stream_alerts_merge_attach_grb_and_trigger(tmp_path):
    store = CandidateStore(tmp_path)
    grb = {"trigger_time": Time(T0_JD + 1.74 / 86400, format="jd").isot + "Z",
           "instrument": "BAT", "id": ["x"], "ra": RA + 1, "dec": DEC, "ra_dec_error": 3.0}
    handle_grb_notice(grb, store)
    # first alert: one night; the host redshift is put in the store as find_host would
    a1 = ztf_alert([0.45, 0.47], [17.3, 17.2], [1, 2], score=None)
    r1 = handle_alert(a1, "ztf", store, assess_kw={})
    assert r1["new"] and not r1["trigger"]
    kn = store.candidates()[0]
    assert kn.grb is not None                                     # stored GRB attached
    kn.redshift = Redshift(0.009783, 1e-5, "spectroscopic", "NGC 4993")
    store.save(kn)
    # second alert, a night later: g and r fade by ~0.7 mag/day
    a2 = ztf_alert([0.45, 0.47, 1.45, 1.47], [17.3, 17.2, 18.0, 17.9], [1, 2, 1, 2], score=None)
    r2 = handle_alert(a2, "ztf", store, assess_kw={})
    assert not r2["new"] and r2["n_photometry"] == 4
    assert r2["trigger"] and r2["trigger_changed"] and "light curve" in r2["kilonova_evidence"]
    assert (tmp_path / "assess" / "ZTF17aaaaaaa_targeted_search_request.json").exists()
    assert store.candidates()[0].redshift is not None             # stored host info kept


def test_stream_loop_with_fake_consumer(tmp_path, monkeypatch):
    import fink_client.consumer as fc
    alerts = [("fink_kn_candidates_ztf", ztf_alert([0.5], [17.3], [1]), "k"), (None, None, None)]

    class Fake:
        def __init__(self, topics, config, survey):
            assert config["username"] == "u" and config["bootstrap.servers"] == "h:1"
            self.q = list(alerts)

        def __enter__(self):
            return self

        def __exit__(self, *a):
            pass

        def poll(self, timeout):
            return self.q.pop(0) if self.q else (None, None, None)
    monkeypatch.setattr(fc, "AlertConsumer", Fake)
    out = []
    n = stream(CandidateStore(tmp_path), "ztf", ["fink_kn_candidates_ztf"],
               credentials={"username": "u", "servers": "h:1", "group_id": "g"},
               on_result=out.append, max_alerts=1)
    assert n == 1 and out[0]["candidate"] == "ZTF17aaaaaaa" and out[0]["topic"] == "fink_kn_candidates_ztf"
