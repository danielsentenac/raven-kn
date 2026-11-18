"""Read-only service: tokens, LVK visibility, steps, events, reduced maps."""
import gzip
import json
import threading
import urllib.error
import urllib.request

import pytest

from raven_kn.journal import RunLog
from raven_kn.listener import AlertStore, add_optical, handle_gw_alert
from raven_kn.serve import Store, make_server, match_external_grb
from test_listener import gw_alert
from test_pipeline import T0, make_kn


@pytest.fixture()
def service(tmp_path):
    store = AlertStore(tmp_path, source="t")
    add_optical(make_kn(), store, now_gps=T0 + 86400)
    handle_gw_alert(gw_alert(), store)
    with RunLog(tmp_path / "runs" / "test" / "pub", "test", "pub", visibility="public") as rl:
        rl.subsegment(1, 1, {"rho_net": 9.0})
        rl.state["p_value"] = 0.01
    with RunLog(tmp_path / "runs" / "test" / "priv", "test", "priv", visibility="lvk") as rl:
        rl.pe("gw", "done", summary={}, warnings=["distance posterior on the lower prior bound"])
    srv = make_server(tmp_path, "127.0.0.1", 0, token="pub-token", lvk_token="lvk-token")
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{srv.server_address[1]}", tmp_path
    srv.shutdown()
    srv.server_close()


def get(base, path, token="pub-token", raw=False):
    req = urllib.request.Request(base + path, headers={"Authorization": f"Bearer {token}"} if token else {})
    with urllib.request.urlopen(req) as r:
        body = r.read()
    return body if raw else json.loads(body)


def test_token_required_and_json_refusal(service):
    base, _ = service
    for tok in (None, "wrong"):
        with pytest.raises(urllib.error.HTTPError) as e:
            get(base, "/api/health", tok)
        assert e.value.code == 401 and json.loads(e.value.read())["service"] == "raven-kn"
    h = get(base, "/api/health")
    assert h["service"] == "raven-kn" and h["lvk"] is False and h["store"]["optical"] == 1
    assert get(base, "/api/health", "lvk-token")["lvk"] is True
    with pytest.raises(urllib.error.HTTPError) as e:
        get(base, "/api/optical/..%2F..%2Fetc")
    assert e.value.code == 404


def test_lvk_runs_hidden_without_the_lvk_token(service):
    base, _ = service
    d = get(base, "/api/optical/test")
    assert [r["run"] for r in d["runs"]] == ["pub"] and d["hidden_runs"] == 1
    steps = {s["key"]: s for s in d["steps"]}
    assert steps["search"]["status"] == "done" and "not shown" in steps["search"]["detail"]
    assert steps["pe"]["status"] == "pending"
    assert steps["coincidence"]["status"] == "done" and "S170817a" in steps["coincidence"]["detail"]
    ev = get(base, "/api/events")["events"]
    assert all(e.get("visibility") != "lvk" for e in ev) and any(e.get("run") == "pub" for e in ev)
    with pytest.raises(urllib.error.HTTPError):
        get(base, "/api/runs/test/priv")
    d = get(base, "/api/optical/test", "lvk-token")
    assert {r["run"] for r in d["runs"]} == {"pub", "priv"}
    assert {s["key"]: s for s in d["steps"]}["pe"]["status"] == "warning"
    assert any(e.get("run") == "priv" for e in get(base, "/api/events", "lvk-token")["events"])
    assert get(base, "/api/runs/test/priv", "lvk-token")["state"]["visibility"] == "lvk"


def test_events_paging_candidates_and_maps(service):
    base, _ = service
    seen, cur = [], ""
    while True:
        r = get(base, f"/api/events?limit=3&cursor={cur}", "lvk-token")
        seen += [e["id"] for e in r["events"]]
        cur = r["cursor"]
        if not r["more"]:
            break
    assert len(seen) == len(set(seen)) == len(get(base, "/api/events", "lvk-token")["events"])
    c = get(base, "/api/optical")["optical"][0]
    assert c["name"] == "test" and c["steps"][0]["status"] == "done"
    gw = get(base, "/api/gw")["superevents"][0]
    assert gw["superevent_id"] == "S170817a" and gw["versions"][0]["maps"]["gw"]["nside"] == 128
    m = gzip.decompress(get(base, "/api/gw/S170817a/map/1/gw", raw=True))
    assert len(m) == 12 * 128 ** 2
    sky = get(base, "/api/sky")
    assert sky["optical"][0]["name"] == "test" and sky["gw"][0]["versions"]
    assert sky["optical"][0]["superevents"] == ["S170817a"] and sky["optical"][0]["grb_id"] is None


def test_run_status_from_state_and_queue(tmp_path):
    d = tmp_path / "runs" / "c" / "dir1"
    RunLog(d, "c", "r1", visibility="public").submitted()
    f = d / "r1_state.json"
    s = Store(tmp_path)
    assert s.run_state(f, False)["status"] == "pending"
    (tmp_path / "runs" / "condor.json").write_text(json.dumps(
        [{"Iwd": "/home/u/runs/dir1", "JobStatus": 5, "HoldReason": "memory exceeded", "ClusterId": 7}]))
    st = s.run_state(f, False)
    assert st["status"] == "held" and "memory" in st["detail"] and st["run"] == "r1"
    # a running job whose heartbeat stopped
    state = json.loads(f.read_text())
    state.update(stage="searching", heartbeat="2026-01-01T00:00:00.000Z")
    f.write_text(json.dumps(state))
    (tmp_path / "runs" / "condor.json").unlink()
    assert s.run_state(f, False)["status"] == "stale"


def test_several_runs_in_one_directory(tmp_path):
    """Jobs sharing an output directory (--out .) keep separate states, journals and results."""
    d = tmp_path / "runs" / "c" / "shared"
    for label, p in (("search_free", 0.5), ("search_circular", 0.01)):
        with RunLog(d, "c", label, visibility="public") as rl:
            rl.state["p_value"] = p
        (d / f"{label}_targeted_result.json").write_text(json.dumps({"p_value": p, "loudest": {}}))
    s = Store(tmp_path)
    runs, hidden = s.runs("c", False)
    assert {r["run"]: r["p_value"] for r in runs} == {"search_free": 0.5, "search_circular": 0.01}
    assert s.run_detail("c", "search_circular", False)["results"][0]["p_value"] == 0.01
    from raven_kn.journal import read_events
    assert {e["run"] for e in read_events(tmp_path)[0]} == {"search_free", "search_circular"}


def test_lvk_gw_grb_coincidence_is_tied_to_the_stored_grb(tmp_path):
    from raven_kn.listener import handle_grb_notice
    from test_listener import external_coinc
    store = AlertStore(tmp_path, source="t")

    def grb(trigger, instrument, gps):
        from astropy.time import Time
        handle_grb_notice({"trigger_time": Time(gps, format="gps").utc.isot + "Z", "ra": 176.8, "dec": -39.8,
                           "ra_dec_error": 11.6, "instrument": instrument, "id": trigger}, store)
    grb("524666471", "Fermi-GBM", T0 + 1.74)
    grb("999", "Swift-BAT", T0 + 1.74)
    grb("777", "Fermi-GBM", T0 + 900)
    handle_gw_alert({**gw_alert(far=1e-6), "external_coinc": external_coinc()}, store)
    s = Store(tmp_path)
    se = s.superevents()[0]
    assert (se["grb_id"], se["grb_match"]) == ("524666471", "notice number")
    assert se["external_coinc"]["observatory"] == "Fermi"
    # after a retraction the alert file is gone; the link comes from the stored version
    handle_gw_alert(gw_alert(alert_type="RETRACTION"), store)
    assert Store(tmp_path).superevents()[0]["grb_id"] == "524666471"

    grbs = {"a": {"grb": {"graceid": "1", "instrument": "Fermi-GBM", "gps": 1000.0 + 1.7}},
            "b": {"grb": {"graceid": "2", "instrument": "Swift-BAT", "gps": 1000.0 + 1.7}},
            "c": {"grb": {"graceid": "3", "instrument": "Fermi-GBM", "gps": 1000.0 + 60}}}
    ec = {"gcn_notice_id": 424242, "observatory": "Fermi", "time_difference": -1.74}
    assert match_external_grb(ec, 1000.0, grbs) == ("a", "mission and time")     # number unknown: mission and time
    assert match_external_grb({**ec, "time_difference": 1.74}, 1000.0, grbs)[0] == "a"   # either sign
    assert match_external_grb({**ec, "observatory": "SVOM"}, 1000.0, grbs) is None
    assert match_external_grb({**ec, "time_difference": -30.0}, 1000.0, grbs) is None
    # the same number under another mission is not a match
    assert match_external_grb({"gcn_notice_id": 2, "observatory": "Fermi"}, 1000.0, grbs) is None
    assert match_external_grb({"gcn_notice_id": 2, "observatory": "Swift"}, 1000.0, grbs) == ("b", "notice number")
    assert match_external_grb(None, 1000.0, grbs) is None


def test_followups_name_alerts_by_type_and_label_candidates(tmp_path):
    from raven_kn.listener import handle_grb_notice
    from test_listener import external_coinc
    from astropy.time import Time
    store = AlertStore(tmp_path, source="t")

    def grb(trigger, instrument, gps):
        handle_grb_notice({"trigger_time": Time(gps, format="gps").utc.isot + "Z", "ra": 176.8, "dec": -39.8,
                           "ra_dec_error": 11.6, "instrument": instrument, "id": trigger}, store)
    grb("524666471", "Fermi-GBM", T0 + 1.74)
    grb("lonely", "Swift-BAT", T0 + 5e5)
    handle_gw_alert({**gw_alert(far=1e-6), "external_coinc": external_coinc()}, store)       # BNS 0.99
    handle_gw_alert(gw_alert(sid="S999bbh", t_gps=T0 + 9e5) | {"event": {
        **gw_alert()["event"], "time": Time(T0 + 9e5, format="gps").utc.isot + "Z",
        "classification": {"BBH": 0.99, "BNS": 0.0, "NSBH": 0.0}}}, store)
    add_optical(make_kn(), store, now_gps=T0 + 86400)
    f = {x["id"]: x for x in Store(tmp_path).followups(True)}
    main = f["optical:test"]
    assert {m["key"] for m in main["members"]} == {"optical:test", "gw:S170817a", "grb:524666471"}
    assert main["candidate"] and main["origin"] == "gw" and main["types"] == ["gw", "grb", "optical"]
    assert any("LVK" in r for r in main["reasons"]) and any("neutron star" in r for r in main["reasons"])
    assert main["steps"] and main["optical"] == "test"
    by = {m["type"]: m["instrument"] for m in main["members"]}
    assert by == {"gw": "H1 L1 V1", "grb": "Fermi-GBM", "optical": "LSST"}
    assert not f["grb:lonely"]["candidate"] and f["grb:lonely"]["reasons"] == []
    assert not f["gw:S999bbh"]["candidate"]                   # a binary black hole alone is not a candidate
    # a neutron-star alert alone is one; retracted, it is not
    store2 = AlertStore(tmp_path / "b", source="t")
    handle_gw_alert(gw_alert(), store2)
    assert Store(tmp_path / "b").followups(True)[0]["candidate"]
    handle_gw_alert(gw_alert(alert_type="RETRACTION"), store2)
    assert not Store(tmp_path / "b").followups(True)[0]["candidate"]
