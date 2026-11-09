"""GCN listener handlers on synthetic igwn.gwalert and GRB unified-schema messages."""
import base64
import io

import healpy as hp
import numpy as np
import pytest
from astropy.table import Table
from astropy.time import Time

from raven_kn.listener import (CandidateStore, add_candidate, handle_grb_notice,
                               handle_gw_alert, moc_lookup, moc_to_flat, uniq_to_order_ipix)
from test_pipeline import T0, make_kn

RA, DEC = 197.450374, -23.381495


def make_moc(ra=RA, dec=DEC, mu=40.0, sigma=8.0, coarse=3, fine=6):
    """Multi-order map: order ``coarse`` everywhere, the 3x3-ish block around
    (ra, dec) refined to order ``fine`` and holding all the probability."""
    hot_coarse = hp.ang2pix(1 << coarse, ra, dec, lonlat=True, nest=True)
    k = 4 ** (fine - coarse)
    uniq, dens = [], []
    for p in range(hp.nside2npix(1 << coarse)):
        if p == hot_coarse:
            for q in range(p * k, (p + 1) * k):
                uniq.append(4 * 4 ** fine + q)
                dens.append(1.0)
        else:
            uniq.append(4 * 4 ** coarse + p)
            dens.append(0.0)
    dens = np.array(dens)
    dens /= np.sum(dens * 4 * np.pi / hp.nside2npix(1 << fine))
    n = len(uniq)
    # r^2-ansatz parameters for a ~Gaussian distance posterior
    return Table({"UNIQ": np.array(uniq, dtype=np.int64), "PROBDENSITY": dens,
                  "DISTMU": np.full(n, mu), "DISTSIGMA": np.full(n, sigma),
                  "DISTNORM": np.full(n, 1 / (mu ** 2 + sigma ** 2))})


def gw_alert(t_gps=T0, far=1e-9, alert_type="PRELIMINARY", moc=None, sid="S170817a"):
    ev = None
    if alert_type != "RETRACTION":
        buf = io.BytesIO()
        (moc if moc is not None else make_moc()).write(buf, format="fits")
        ev = {"time": Time(t_gps, format="gps").utc.isot + "Z", "far": far, "significant": True,
              "instruments": ["H1", "L1", "V1"], "group": "CBC", "pipeline": "gstlal",
              "search": "AllSky", "classification": {"BNS": 0.99},
              "skymap": base64.b64encode(buf.getvalue()).decode()}
    return {"superevent_id": sid, "alert_type": alert_type,
            "time_created": "2017-08-17T12:41:30Z", "event": ev}


def test_moc_decoding():
    moc = make_moc()
    order, ipix = uniq_to_order_ipix(moc["UNIQ"])
    assert set(order) == {3, 6}
    flat = moc_to_flat(moc, max_order=5)        # finer pixels merged
    assert flat["PROB"].sum() == pytest.approx(1, rel=1e-9)
    assert np.isfinite(flat["DISTMU"][flat["PROB"] > 0]).all()
    assert moc_lookup(moc, RA, DEC)["DISTMU"] == 40.0
    assert moc_lookup(moc, RA, DEC)["PROBDENSITY"] > 50     # 1 / order-3 pixel area


def test_gw_alert_in_window_gives_joint_far(tmp_path):
    store = CandidateStore(tmp_path)
    store.save(make_kn())
    far_kn = make_kn(); far_kn.name = "elsewhere"; far_kn.ra += 120
    store.save(far_kn)
    res = handle_gw_alert(gw_alert(), store)
    by = {r["external_event"]: r for r in res}
    assert set(by) == {"test", "elsewhere"}       # both windows contain T0
    near, off = by["test"], by["elsewhere"]
    assert near["skymap_overlap"] > 100 and off["skymap_overlap"] == 0
    assert near["spatiotemporal_distance_coinc_far"] < near["spatiotemporal_coinc_far"] < near["temporal_coinc_far"]
    assert near["distance_overlap_valid"] is True and near["gw_map_at_kn"]["DISTMU"] == 40.0
    assert 10 < near["dt_from_first_detection_h"] < 12
    assert (tmp_path / "coinc" / "S170817a__test.json").exists()


def test_gw_alert_outside_window_ignored_and_retraction(tmp_path):
    store = CandidateStore(tmp_path)
    store.save(make_kn())
    assert handle_gw_alert(gw_alert(t_gps=T0 - 5 * 86400, sid="S1"), store) == []
    assert (tmp_path / "gw" / "S1.json").exists()
    handle_gw_alert(gw_alert(alert_type="RETRACTION", sid="S1"), store)
    assert not list((tmp_path / "gw").glob("S1.*"))


def test_late_candidate_matched_to_stored_alert(tmp_path):
    store = CandidateStore(tmp_path)
    assert handle_gw_alert(gw_alert(), store) == []          # no candidate yet
    res = add_candidate(make_kn(), store, now_gps=T0 + 2 * 86400)
    assert len(res) == 1 and res[0]["skymap_overlap"] > 100
    assert add_candidate(make_kn(), store, now_gps=T0 + 30 * 86400) == []   # alert too old


def test_grb_notice_attached_and_narrows_window(tmp_path):
    store = CandidateStore(tmp_path)
    kn = make_kn(); store.save(kn)
    notice = {"trigger_time": Time(T0 + 1.74, format="gps").utc.isot + "Z",
              "instrument": "BAT", "id": ["x"], "ra": RA + 1.0, "dec": DEC, "ra_dec_error": 3.0}
    assert handle_grb_notice(notice, store) == ["test"]
    kn2 = store.candidates()[0]
    assert kn2.grb.gps == pytest.approx(T0 + 1.74, abs=1e-3)
    far = {**notice, "ra": RA + 60}
    kn3 = make_kn(); kn3.name = "other"; store.save(kn3)
    assert handle_grb_notice(far, store) == []
    res = handle_gw_alert(gw_alert(), store)
    w = {r["external_event"]: r["window"] for r in res}
    assert w["test"]["gps_end"] - w["test"]["gps_start"] == pytest.approx(6.0)
