"""Logbook: journal and cursor, heartbeat, job state, sky-map versions, handler entries."""
import gzip
import json
from datetime import datetime, timedelta, timezone

import numpy as np
import pytest

from raven_kn.journal import (CRED_SCALE, Heartbeat, Journal, RunLog, credible_levels, decode_cursor,
                              log_decision, pe_warnings, read_events, reduce_map, sky_versions)
from raven_kn.listener import AlertStore, add_optical, handle_grb_notice, handle_gw_alert
from test_listener import RA, DEC, gw_alert, make_moc
from test_pipeline import T0, make_kn


def test_journal_cursor_is_exact_across_files(tmp_path):
    a, b = Journal.for_store(tmp_path, "a"), Journal.for_store(tmp_path, "b")
    t = [datetime(2026, 1, 1, tzinfo=timezone.utc)]

    def clock():
        t[0] += timedelta(seconds=1)
        return t[0]
    a.clock = b.clock = clock
    for i in range(5):
        (a if i % 2 else b).log("k", f"entry {i}", n=i, x=np.float64(1.5), bad=float("nan"))
    ev, cur, more = read_events(tmp_path)
    assert [e["data"]["n"] for e in ev] == [0, 1, 2, 3, 4] and not more
    assert ev[0]["data"] == {"n": 0, "x": 1.5, "bad": None} and ev[0]["id"] == "b:0"
    assert read_events(tmp_path, cur)[0] == []
    b.log("k", "late", n=5)
    # an older entry arriving later in another file (a job journal copied from a cluster)
    late = tmp_path / "runs" / "c" / "r" / "x_journal.jsonl"
    late.parent.mkdir(parents=True)
    late.write_text(json.dumps({"time": "2020-01-01T00:00:00.000Z", "kind": "old", "summary": "s",
                                "visibility": "lvk", "data": {}}) + "\n" + '{"incomplete": ')
    ev2, cur2, _ = read_events(tmp_path, cur)
    assert [e["kind"] for e in ev2] == ["old", "k"] and ev2[0]["id"] == "run/c/r/x:0"
    assert read_events(tmp_path, cur, allow_lvk=False)[0][0]["kind"] == "k"
    assert read_events(tmp_path, cur2)[0] == []
    # paging: every entry exactly once
    seen, c = [], None
    while True:
        page, c, more = read_events(tmp_path, c, limit=2)
        seen += [e["id"] for e in page]
        if not more:
            break
    assert sorted(seen) == sorted(e["id"] for e in read_events(tmp_path)[0]) and len(set(seen)) == 7
    assert decode_cursor("garbage") == {}


def test_heartbeat_reports_the_outage(tmp_path):
    j = Journal.for_store(tmp_path, "listen")
    hb = Heartbeat(tmp_path, "listen", j, ["t"])
    hb.message("t")
    hb.kafka_error("down")
    hb.kafka_error("down")                         # same text: counted, not written again
    hb.write(force=True)                           # the file is rewritten every 30 s
    st = json.loads(hb.path.read_text())
    st["heartbeat"] = "2026-01-01T00:00:00.000Z"
    hb.path.write_text(json.dumps(st))
    hb2 = Heartbeat(tmp_path, "listen", j, ["t"])
    ev = read_events(tmp_path)[0]
    assert [e["kind"] for e in ev] == ["listener_start", "kafka_error", "listener_start"]
    assert ev[2]["level"] == "warning" and ev[2]["data"]["down_s"] > 86400
    assert hb2.state["restarts"] == 1 and hb.state["kafka_errors"] == 2
    assert "t" in hb2.state["last_message"]         # kept across the restart


def test_runlog_states(tmp_path):
    rl = RunLog(tmp_path / "ok", "cand", "ok", visibility="public", heartbeat_s=0.05)
    rl.submitted("note")
    with rl:
        rl.subsegment(1, 2, {"rho_net": 5.0})
        rl.subsegment(2, 2, {"rho_net": 4.0})
        rl.pe("gw", "done", summary={}, warnings=["w"])
    st = json.loads((tmp_path / "ok" / "ok_state.json").read_text())
    assert (st["exit"], st["stage"], st["subsegments_done"], st["loudest"]["rho_net"]) == ("ok", "done", 2, 5.0)
    assert st["submitted"] and st["pe"]["gw"]["warnings"] == ["w"]
    with pytest.raises(ValueError):
        with RunLog(tmp_path / "bad", "cand", "bad"):
            raise ValueError("no data")
    st = json.loads((tmp_path / "bad" / "bad_state.json").read_text())
    assert st["exit"] == "failed" and "no data" in st["error"]
    kinds = [json.loads(l)["kind"] for l in (tmp_path / "bad" / "bad_journal.jsonl").read_text().splitlines()]
    assert kinds == ["run_start", "run_failed"]


def test_pe_warnings():
    good = {"maximum_likelihood": {"snr": 30.0}, "log_bayes_factor": 440.0, "distance_at_prior_edge": None}
    assert pe_warnings(good, 30.3) == []
    bad = {"maximum_likelihood": {"snr": 32.0}, "log_bayes_factor": 13942.0, "distance_at_prior_edge": "lower"}
    assert len(pe_warnings(bad, 30.0)) == 2
    assert len(pe_warnings({"maximum_likelihood": {"snr": 12.0}, "log_bayes_factor": 50.0}, 30.0)) == 1


def test_reduced_map():
    c = credible_levels(np.array([0.1, 0.6, 0.3]))
    assert c == pytest.approx([1.0, 0.6, 0.9])
    data, s = reduce_map(make_moc())
    b = np.frombuffer(gzip.decompress(data), dtype=np.uint8)
    assert len(b) == 12 * 128 ** 2 and s["nside"] == 128
    import healpy as hp
    assert b[hp.ang2pix(128, RA, DEC, lonlat=True, nest=True)] <= CRED_SCALE
    assert (b == 255).mean() > 0.99                  # all the probability in one order-3 pixel
    assert s["area90_deg2"] == pytest.approx(0.9 * hp.nside2pixarea(8, degrees=True), rel=0.02)


def test_handlers_write_the_logbook(tmp_path):
    store = AlertStore(tmp_path, source="t")
    add_optical(make_kn(), store, now_gps=T0 + 86400)
    add_optical(make_kn(), store, now_gps=T0 + 86400)           # unchanged: no entry
    handle_gw_alert(gw_alert(), store)
    handle_gw_alert(gw_alert(alert_type="INITIAL", far=1e-10), store)
    handle_gw_alert(gw_alert(alert_type="RETRACTION"), store)
    ev = read_events(tmp_path)[0]
    assert [e["kind"] for e in ev] == ["optical_new", "gw_alert", "coincidence", "gw_alert",
                                       "coincidence", "gw_retraction"]
    co = ev[2]
    assert co["optical"] == "test" and co["superevent"] == "S170817a"
    assert co["data"]["result"]["skymap_overlap"] > 100 and co["data"]["joint_far"] > 0
    # the store forgot the superevent; the versions and the journal did not
    assert not (tmp_path / "gw" / "S170817a.json").exists()
    v = sky_versions(tmp_path, "S170817a")
    assert [x["alert_type"] for x in v] == ["PRELIMINARY", "INITIAL"]
    assert (tmp_path / "gw_versions" / "S170817a" / "02_gw.cred.gz").exists()
    assert (tmp_path / "gw_versions" / "S170817a" / "retracted").exists()


def test_decision_logged_only_when_it_changes(tmp_path):
    j = Journal.for_store(tmp_path, "t")
    d = {"optical": "c", "trigger_targeted_search": False,
         "kilonova_evidence": {"confirmed": False, "by": [], "required": True},
         "detectability": {"worth_search": True}, "t0_window": {"method": "GRB"}}
    assert log_decision(j, d, None)["data"]["first"] is True
    assert log_decision(j, d, d) is None
    d2 = {**d, "trigger_targeted_search": True, "kilonova_evidence": {"confirmed": True, "by": ["light curve"]}}
    e = log_decision(j, d2, d)
    assert e["data"]["trigger_changed"] and "ready for a targeted search" in e["summary"]


def test_data_written_before_the_rename_is_still_read(tmp_path):
    """Store folder ``candidates/``, journal key ``candidate`` and kinds
    ``candidate_*``, request key ``candidate``: read; new files use the new names."""
    import raven_kn.connection as cn
    import raven_kn.listener as ls
    from raven_kn.connection import save_optical
    from raven_kn.serve import Store
    assert cn.KNCandidate is cn.OpticalTransient and ls.CandidateStore is ls.AlertStore   # old names kept
    (tmp_path / "candidates").mkdir()
    old = make_kn(); old.name = "old"
    save_optical(old, tmp_path / "candidates" / "old.json")
    (tmp_path / "journal").mkdir()
    (tmp_path / "journal" / "listen.jsonl").write_text(json.dumps(
        {"time": "2026-11-01T00:00:00.000Z", "kind": "candidate_new", "summary": "old: new candidate",
         "candidate": "old", "data": {}}) + "\n")
    store = AlertStore(tmp_path, source="t")
    assert [k.name for k in store.transients()] == ["old"]
    add_optical(make_kn(), store, now_gps=T0 + 86400)              # a new one: written in optical/
    assert (tmp_path / "optical" / "test.json").exists() and not (tmp_path / "candidates" / "test.json").exists()
    old.notes = "updated"
    add_optical(old, store, now_gps=T0 + 86400)                    # an old one updated: its new copy wins
    assert (tmp_path / "optical" / "old.json").exists()
    assert {k.name: k.notes for k in store.transients()}["old"] == "updated"
    ev = read_events(tmp_path)[0]
    assert ev[0]["kind"] == "optical_new" and ev[0]["optical"] == "old" and "candidate" not in ev[0]
    assert Store(tmp_path).optical_names() == ["old", "test"]
    # a search request written before the rename still names its run
    from raven_kn.targeted.search import TargetedSearchConfig
    req = tmp_path / "req.json"
    req.write_text(json.dumps({"candidate": "SSS17a", "ra_deg": 1.0, "dec_deg": 2.0, "gps_start": 10.0, "gps_end": 20.0}))
    assert TargetedSearchConfig.from_request(req).label == "SSS17a"
    req.write_text(json.dumps({"optical": "AT2017gfo", "ra_deg": 1.0, "dec_deg": 2.0, "gps_start": 10.0, "gps_end": 20.0}))
    assert TargetedSearchConfig.from_request(req).label == "AT2017gfo"
