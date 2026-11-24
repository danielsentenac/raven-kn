"""A transient ready for a targeted search is announced once: journal and mail."""
import json

from raven_kn.journal import read_events
from raven_kn.listener import AlertStore
from raven_kn.notify import announced, on_decision, send_mail
from raven_kn.serve import Store


def decision(ready, name="KN-1"):
    return {"optical": name, "ra": 197.45, "dec": -23.38, "trigger_targeted_search": ready,
            "kilonova_evidence": {"confirmed": ready, "by": ["light curve"] if ready else [],
                                  "reasons": ["g fades by 0.70 mag/day"] if ready else ["one epoch only"]},
            "detectability": {"worth_search": True, "snr_median": 12.3, "p_detect": 0.9, "reasons": []},
            "distance_inclination": {"distance_informed": True,
                                     "luminosity_distance_mpc": {"median": 41.0, "lower": 38.0, "upper": 44.0}},
            "t0_window": {"duration_h": 12.0, "method": "bracket", "gps_start": 1187000000.0,
                          "gps_end": 1187043200.0}}


def test_announced_once_then_withdrawn(tmp_path):
    store, sent = AlertStore(tmp_path, source="t"), []

    def mail(subject, body):
        sent.append((subject, body))
        return True
    assert on_decision(store, decision(False), mail) is None            # never ready: nothing to say
    e = on_decision(store, decision(True), mail)
    assert e["kind"] == "announcement" and e["level"] == "warning" and e["data"]["mailed"]
    assert "ready for a targeted search" in e["summary"] and len(sent) == 1
    subject, body = sent[0]
    assert subject == "[RAVEN-KN] KN-1: ready for a targeted search (confirmed as a kilonova, within reach)"
    assert "No search has been started" in body and "Distance: 41 Mpc (38 to 44)" in body
    assert "Kilonova evidence: light curve" in body and "43200 s from bracket" in body
    assert announced(tmp_path, "KN-1")["state"] == "ready"
    assert on_decision(store, decision(True), mail) is None and len(sent) == 1     # the same: not again
    w = on_decision(store, decision(False), mail)
    assert "no longer ready" in w["summary"] and len(sent) == 2 and "one epoch only" in sent[1][1]
    assert on_decision(store, decision(False), mail) is None
    assert on_decision(store, decision(True), mail) is not None and len(sent) == 3   # ready again
    rec = announced(tmp_path, "KN-1")
    assert [h["state"] for h in rec["history"]] == ["ready", "withdrawn", "ready"]
    kinds = [x["kind"] for x in read_events(tmp_path, None, 100)[0]]
    assert kinds.count("announcement") == 3


def test_no_address_no_mail_but_still_in_the_journal(tmp_path, monkeypatch):
    monkeypatch.delenv("RAVEN_KN_MAIL_TO", raising=False)
    store = AlertStore(tmp_path, source="t")
    e = on_decision(store, decision(True))
    assert e["data"]["mailed"] is False and announced(tmp_path, "KN-1")["mailed"] is False
    # the address is a setting of the machine: it is in no file of the store
    monkeypatch.setenv("RAVEN_KN_MAIL_TO", "someone@example.org, other@example.org")
    calls = []

    class Done:
        returncode = 0

    def run(cmd, input=None, **kw):
        calls.append((cmd, input.decode()))
        return Done()
    assert send_mail("subject line", "text", run=run)
    assert calls[0][0][1:] == ["-t", "-oi"] and "To: someone@example.org, other@example.org" in calls[0][1]
    assert "Subject: subject line" in calls[0][1]
    on_decision(store, decision(True, "KN-2"), lambda s, b: True)
    assert not any("example.org" in p.read_text() for p in tmp_path.rglob("*") if p.is_file())

    def broken(cmd, **kw):
        raise OSError("no sendmail")
    assert send_mail("s", "t", run=broken) is False                      # never raises


def test_the_service_says_which_transients_wait_for_a_person(tmp_path):
    store = AlertStore(tmp_path, source="t")
    (tmp_path / "optical").mkdir(exist_ok=True)
    (tmp_path / "optical" / "KN-1.json").write_text(json.dumps({"name": "KN-1", "ra": 197.45, "dec": -23.38}))
    (tmp_path / "assess").mkdir()
    (tmp_path / "assess" / "KN-1_decision.json").write_text(json.dumps(decision(True)))
    on_decision(store, decision(True), lambda s, b: True)
    view = Store(tmp_path)
    s = view.optical_summary("KN-1", False)
    assert s["ready"] and s["announced"]["state"] == "ready"
    f = next(f for f in view.followups(False) if f["optical"] == "KN-1")
    assert f["ready"] and f["candidate"]
    step = {x["key"]: x for x in s["steps"]}
    assert step["decision"]["detail"].startswith("ready for a targeted search")
    assert "a person starts it" in step["search"]["detail"]
    # a search asked for: it no longer waits
    (tmp_path / "launch").mkdir()
    (tmp_path / "launch" / "x.json").write_text(json.dumps({"optical": "KN-1", "status": "running",
                                                            "requested": "2026-11-24T10:00:00.000Z"}))
    assert not view.optical_summary("KN-1", False)["ready"]
