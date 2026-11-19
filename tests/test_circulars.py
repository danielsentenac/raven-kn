"""GCN Circulars: links to the alerts of the store and labels, on real
circulars of November 2026 (``tests/data/circulars``, e-mail addresses
removed) and the notices they are about (times and trigger numbers of the
production store)."""
import json
import threading
import urllib.request
from pathlib import Path

import pytest

from raven_kn.circulars import (classify, event_date, expand_ids, handle_circular, parse_circular,
                                stored_circulars, summary, times_in)
from raven_kn.connection import GRBAssociation, NeutrinoAssociation
from raven_kn.journal import read_events
from raven_kn.listener import AlertStore
from raven_kn.serve import Store, make_server

DATA = Path(__file__).parent / "data" / "circulars"


def circ(number):
    return json.loads((DATA / f"{number}.json").read_text())


@pytest.fixture()
def store(tmp_path):
    s = AlertStore(tmp_path, source="t")
    for instrument, trigger, gps in (("Fermi-GBM", "816089232", 1478431645.95),      # GRB 261111A
                                     ("Fermi-GBM", "816645341", 1478987754.47),
                                     ("SVOM-ECLAIRs", "sb26111302", 1478615477.667),  # GRB 261113A
                                     ("WXT", None, 1478720453.649)):                  # EP261114a
        s.save_grb({"grb": GRBAssociation(gps=gps, instrument=instrument, graceid=trigger),
                    "ra": 10.0, "dec": 20.0, "err_deg": 1.0})
    s.save_nu({"nu": NeutrinoAssociation(gps=1478600000.0, instrument="IceCube",
                                         event_name="IceCube-261113A", ra=1.0, dec=2.0, err_deg=1.0)})
    return s


def test_parse_keeps_no_email():
    c = parse_circular(circ(46244))
    assert c["id"] == 46244 and c["event"] == "GRB 261113A"
    assert c["time"].startswith("2026-11-13T") and c["url"].endswith("/circulars/46244")
    assert "@" not in c["submitter"] and "<" not in c["submitter"]


@pytest.mark.parametrize("subject, bands, result", [
    ("GRB 261113A: JinShan optical afterglow detection", ["optical"], "detection"),
    ("EP261114a: NOT Optical Upper Limits", ["optical"], "upper_limit"),
    ("EP261115a: FRAM-Auger optical limits", ["optical"], "upper_limit"),
    # "Akeno optical afterglow" holds the letters of "no optical afterglow"
    ("GRB 261113A : MITSuME Akeno optical afterglow detection", ["optical"], "detection"),
    ("GRB 261113A : GRANDMA and KNC optical detections and upper limits", ["optical"], "detection"),
    ("IceCube-260919A: No transient candidates from the Zwicky Transient Facility", ["optical"],
     "upper_limit"),
    ("GRB 261113A: NOT optical observations", ["optical"], "observation"),
    ("GRB 261010A: Enhanced Swift-XRT position", ["xray"], None),
    ("EP261006a: 10 GHz VLA observations", ["radio"], None),
    ("GRB 261111A: Fermi GBM Final Real-time Localization", ["gamma"], None),
    ("GRB 261113A: SVOM detection of a bright burst", ["gamma"], None),
    ("LIGO/Virgo/KAGRA S251112cm: Identification of a GW compact binary merger candidate", ["gw"], None),
    ("IceCube-261113A - IceCube observation of a high-energy neutrino candidate track-like event",
     ["nu"], None),
    ("GRB 260903A", [], None),
])
def test_classify_subject(subject, bands, result):
    c = classify(subject)
    assert (c["bands"], c["result"]) == (bands, result)


def test_classify_from_the_text_and_redshift():
    # the subject says "observations report": the text gives the upper limits
    assert classify(circ(46246)["subject"], circ(46246)["body"])["result"] == "upper_limit"
    c = classify(circ(46261)["subject"], circ(46261)["body"])
    assert c["result"] == "detection" and c["redshift"] == 3.58
    assert classify("GRB 261113A: optical observations", "photometric redshift z = 1.2")["redshift"] is None


def test_times():
    day = event_date("GRB 261113A")
    assert day.isoformat() == "2026-11-13T00:00:00+00:00"
    assert event_date("EP261114a").day == 14 and event_date("no date") is None
    got = times_in("At 05:48:03 UT the burst... began at 2026-11-13T07:09:35 (UTC); 1:2 ratio", day)
    assert [t.strftime("%d %H:%M:%S") for t in got] == ["13 07:09:35", "13 05:48:03"]
    assert expand_ids(["46261", "46270-46272"]) == [46261, 46270, 46271, 46272]


def test_links(store):
    # trigger number in the text (Fermi, SVOM) and in the subject
    assert handle_circular(circ(46224), store) == [
        {"type": "grb", "id": "816089232", "how": "trigger number"}]
    assert handle_circular(circ(46244), store)[0] == {
        "type": "grb", "id": "sb26111302", "how": "trigger number"}
    assert handle_circular(circ(46315), store)[0]["id"] == "816645341"
    # the Einstein Probe notice has no trigger number: the time of the text
    assert handle_circular(circ(46284), store) == [
        {"type": "grb", "id": "WXT_1478720453.65", "how": "trigger time"}]
    # follow-up reports only give the name learnt from the circulars above
    for n, burst in ((46248, "sb26111302"), (46261, "sb26111302"), (46293, "WXT_1478720453.65")):
        assert handle_circular(circ(n), store) == [{"type": "grb", "id": burst, "how": "event name"}]
    assert handle_circular(circ(46275), store) == [
        {"type": "nu", "id": "IceCube-261113A", "how": "event name"}]
    events = [e for e in read_events(store.root)[0] if e["kind"] == "circular"]
    assert len(events) == 8
    e = next(e for e in events if e["data"]["circular"] == 46261)
    assert e["data"]["redshift"] == 3.58 and e["data"]["bands"] == ["optical"]
    assert "GCN 46261" in e["summary"] and "grb sb26111302" in e["summary"]
    # the same circular again (a correction): stored, not written twice
    handle_circular(circ(46261), store)
    assert len([e for e in read_events(store.root)[0] if e["kind"] == "circular"]) == 8


def test_report_before_the_circular_that_names_the_burst(store):
    # the MASTER report on "GRB 261111A" read before the Fermi circular with the trigger number
    assert handle_circular(circ(46229), store) == []
    assert [c["id"] for c in stored_circulars(store.root, linked=False)] == [46229]
    assert not [e for e in read_events(store.root)[0] if e["kind"] == "circular"]
    handle_circular(circ(46224), store)
    linked = {c["id"]: c for c in stored_circulars(store.root)}
    assert linked[46229]["links"] == [{"type": "grb", "id": "816089232", "how": "event name"}]
    assert stored_circulars(store.root, linked=False) == []
    # a circular about an event that is not in the store, and one without event name
    other = {**circ(46248), "circularId": 1, "eventId": "GRB 990101A",
             "subject": "GRB 990101A: optical observations", "body": "nothing"}
    assert handle_circular(other, store) == []
    assert handle_circular({**other, "circularId": 2, "eventId": None}, store) == []
    assert [c["id"] for c in stored_circulars(store.root, linked=False)] == [1]


def test_short_identifiers_and_longer_numbers_are_not_links(store):
    store.save_grb({"grb": GRBAssociation(gps=1.0e9, instrument="X", graceid="12345"),
                    "ra": 0.0, "dec": 0.0, "err_deg": 1.0})
    msg = {"circularId": 5, "subject": "note", "createdOn": 0, "submitter": "a",
           "body": "trigger 12345 and trigger 8160892321 and 1816089232"}
    assert handle_circular(msg, store) == []


def test_summary_and_service(store):
    for n in (46244, 46246, 46248, 46261, 46304, 46284, 46293):
        handle_circular(circ(n), store)
    s = summary(c for c in stored_circulars(store.root)
                if c["links"][0]["id"] == "sb26111302")
    assert s["n"] == 5 and s["bands"] == {"gamma": 1, "optical": 4} and s["name"] == "GRB 261113A"
    assert s["optical"] == "detection" and s["n_optical"] == 4 and s["redshift"] == 3.58
    assert s["first_optical"].startswith("2026-11-13T")
    view = Store(store.root)
    sky = {g["id"]: g for g in view.sky()["grb"]}
    assert sky["sb26111302"]["circulars"]["optical"] == "detection"
    assert sky["WXT_1478720453.65"]["circulars"]["optical"] == "upper_limit"
    assert sky["816089232"]["circulars"] is None                 # nothing read about it
    member = next(m for f in view.followups(False) for m in f["members"] if m["key"] == "grb:sb26111302")
    assert member["circulars"]["n"] == 5
    assert view.health(False)["store"]["circulars"] == 7

    srv = make_server(store.root, "127.0.0.1", 0, token="tok")
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        def get(path):
            req = urllib.request.Request(f"http://127.0.0.1:{srv.server_address[1]}{path}",
                                         headers={"Authorization": "Bearer tok"})
            with urllib.request.urlopen(req) as r:
                return json.loads(r.read())
        listed = get("/api/circulars")["circulars"]
        assert [c["id"] for c in listed] == [46244, 46246, 46248, 46261, 46284, 46293, 46304]
        assert all("body" not in c for c in listed)
        assert "X-shooter" in get("/api/circulars/46261")["body"]
        with pytest.raises(urllib.error.HTTPError):
            get("/api/circulars/1")
    finally:
        srv.shutdown()
        srv.server_close()
