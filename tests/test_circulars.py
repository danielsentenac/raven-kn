"""GCN Circulars: links to the alerts of the store and labels, on real
circulars of November 2026 (``tests/data/circulars``, e-mail addresses
removed) and the notices they are about (times and trigger numbers of the
production store)."""
import json
import threading
import urllib.request
from pathlib import Path

import pytest

from raven_kn.circulars import (accept_photometry, classify, distance_gate, event_date, expand_ids,
                                handle_circular, parse_circular, photometry_in, position_in,
                                reported_photometry, stored_circulars, summary, times_in)
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


def test_distance_gate(store):
    far = distance_gate(3.716, 46290)
    assert not far["within_reach"] and 33000 < far["d_mpc"] < 34500
    assert far["text"] == ("redshift 3.716 (GCN 46290), luminosity distance 33.6 Gpc: "
                           "too far for a GW search (limit 500 Mpc)")
    near = distance_gate(0.0098)                                  # NGC 4993
    assert near["within_reach"] and 40 < near["d_mpc"] < 45
    assert near["text"].endswith("within the range of a GW search (500 Mpc)")
    # the redshift of GRB 261113A is written in the journal, and its follow-up is marked
    for n in (46244, 46261):
        handle_circular(circ(n), store)
    e = [e for e in read_events(store.root)[0] if e["kind"] == "grb_distance"]
    assert len(e) == 1 and e[0]["data"]["triggers"] == ["sb26111302"]
    handle_circular(circ(46261), store)                           # read again: not written twice
    assert len([e for e in read_events(store.root)[0] if e["kind"] == "grb_distance"]) == 1
    assert e[0]["summary"].startswith("GRB 261113A (sb26111302): redshift 3.58 (GCN 46261)")
    f = next(f for f in Store(store.root).followups(False) if f["id"] == "grb:sb26111302")
    assert f["too_far"] and not f["candidate"] and "too far for a GW search" in f["notes"][0]
    # a nearby burst is a candidate
    near_msg = {**circ(46261), "circularId": 99, "eventId": "GRB 261111A",
                "subject": "GRB 261111A: spectroscopic redshift z = 0.05", "body": "z = 0.05"}
    handle_circular(circ(46224), store)
    handle_circular(near_msg, store)
    f = next(f for f in Store(store.root).followups(False) if f["id"] == "grb:816089232")
    assert f["candidate"] and not f["too_far"]
    assert f["reasons"][0].startswith("GRB 261111A: redshift 0.05 (GCN 99), luminosity distance 2")
    assert [e["level"] for e in read_events(store.root)[0] if e["kind"] == "grb_distance"] == ["info", "warning"]


def test_position_in_the_text():
    p = position_in("candidate is R.A. (J2000) = 23h03m15.21s\nDec. (J2000) = -32d19m01.2s\n\n"
                    "with an uncertainty of 0.50 arcsec.")
    assert p == {"ra": 345.813375, "dec": -32.317, "err_arcsec": 0.5}
    assert position_in("localization (RA, Dec = 3.984, +51.5986 deg, J2000)") == {
        "ra": 3.984, "dec": 51.5986, "err_arcsec": 1.0}
    assert position_in("RA(J2000) = 12:30:49.42, Dec(J2000) = +12:23:28.0")["ra"] == pytest.approx(187.70592, abs=1e-5)
    assert position_in("RA = 187.7059 deg, Dec = -12.3911 deg")["dec"] == -12.3911
    assert position_in("no coordinates here, r = 21.2") is None


def test_photometry_in_the_text():
    t0 = event_date("GRB 261114A").replace(hour=9, minute=31, second=29)
    body = ("Observations started at 20:54:24 UT (about 11.5 hr after the SVOM trigger).\n\n"
            "In our stacked image, we measure the preliminary AB magnitude:\n\nr = 22.38 ± 0.13\n\n"
            "No source is detected in the g band, with an upper limit g > 23.5 (3 sigma). "
            "The afterglow at redshift z = 3.716 was J = 19.37 +/- 0.14 mag (Vega), at a mid-time of "
            "~6.17 hours post-trigger.\n\n| a | r = 18.1 | table |")
    pts = photometry_in(body, t0)
    assert [(p["band"], p["mag"], p["mag_err"], p["upper_limit"]) for p in pts] == [
        ("r", 22.38, 0.13, False), ("g", 23.5, None, True), ("J", 19.37, 0.14, False)]
    assert pts[0]["time"] == "2026-11-14T21:01:29.000Z" and pts[0]["time_from"] == "delay after the trigger"
    assert pts[2]["time"] == "2026-11-14T15:41:41.000Z" and pts[2]["system"] == "Vega"
    assert pts[0]["text"] == "r = 22.38 ± 0.13"
    # without the trigger time only a date of the text gives the time
    assert photometry_in("We measure r = 20.1 +/- 0.1, 3 hours after the trigger.")[0]["time"] is None
    assert photometry_in("On 2026-11-14T10:00:00 UT we measure r = 20.1.")[0]["time_from"] == "date in the text"
    assert photometry_in("The source was detected with L = 16.32 ± 0.01 AB mag (+0.5h).", t0)[0]["time"] \
        == "2026-11-14T10:01:29.000Z"


def test_magnitude_some_time_after_the_discovery():                       # 46276
    body = ("SVOM/GWAC detected an uncatalogued fading source G261113_D02737 at 2026-11-13T19:45:45 UTC. "
            "The source initially rose and then faded, reaching a peak magnitude of W~13.8 mag.\n\n"
            "Automatically follow-up observations by GWAC-F60 showed that it was fading continuously with the "
            "decay slope of -0.94 with a single power law, down to the R~18.5 mag at 3 ks after the GWAC "
            "discovery time.")
    r = next(p for p in photometry_in(body) if p["band"] == "R")
    assert r["time"] == "2026-11-13T20:35:45.000Z"                       # not the discovery time itself
    assert r["time_from"] == "delay after the discovery time of the text"


def test_optical_transient_of_a_reported_counterpart(store):
    from raven_kn.connection import Redshift
    sub = "GRB 261113A: optical afterglow detection"
    base = {**circ(46248), "eventId": "GRB 261113A"}
    handle_circular(circ(46244), store)                          # names the burst sb26111302 (10, 20; 1 deg)
    # a detection without a position, a stellar flare, a position outside the burst: no transient
    handle_circular({**base, "circularId": 1, "subject": sub, "body": "We detect it at r = 19.1 +/- 0.1."}, store)
    handle_circular({**base, "circularId": 2, "subject": sub, "body":
                     "RA, Dec = 10.2, 20.1. We detect a stellar flare of a known M-dwarf."}, store)
    handle_circular({**base, "circularId": 3, "subject": sub, "body": "We detect it at RA, Dec = 50.0, -20.0."}, store)
    assert store.transients() == []
    # a spectroscopic redshift read before the transient exists is applied when it is created
    handle_circular(circ(46261), store)
    handle_circular({**base, "circularId": 4, "subject": sub, "body":
                     "We detect the afterglow at RA, Dec = 10.2001, 20.1002 with an uncertainty of 0.3 arcsec, "
                     "with r = 19.14 +/- 0.05, 0.5 hours after the trigger.\n\nA\n\nB\n\nC\n\nLater, i > 21.0."}, store)
    kn, = store.transients()
    assert (kn.name, kn.ra, kn.dec, kn.pos_err_arcsec) == ("GRB261113A-OT", 10.2001, 20.1002, 0.3)
    assert kn.source == "GCN circular 4" and kn.grb.graceid == "sb26111302" and kn.photometry == []
    assert kn.redshift == Redshift(z=3.58, z_err=1e-3, kind="spectroscopic", source="GCN circular 46261")
    kinds = [e["kind"] for e in read_events(store.root)[0]]
    assert kinds.count("optical_new") == 1
    # the search decision was run: held, too far; its sky map is not kept
    dec = json.loads((store.root / "assess" / "GRB261113A-OT_decision.json").read_text())
    assert not dec["trigger_targeted_search"] and not dec["detectability"]["worth_search"]
    assert not (store.root / "assess" / "GRB261113A-OT_skymap.fits").exists()
    held = [e for e in read_events(store.root)[0] if e["kind"] == "decision"]
    assert len(held) == 1 and held[0]["summary"].startswith("GRB261113A-OT: held")
    steps = {s["key"]: s for s in Store(store.root).steps("GRB261113A-OT", False)}
    assert steps["decision"]["status"] == "held"
    # the optical reports about the burst are linked to its transient: the one that created it,
    # those read before (1, 46261), not the gamma-ray one (46244) nor the report on a flare star... (2 is optical)
    by = {c["id"]: [(ln["type"], ln["id"], ln["how"]) for ln in c["links"]] for c in stored_circulars(store.root)}
    link = ("optical", "GRB261113A-OT", "optical report about its burst")
    assert all(link in by[n] for n in (1, 2, 3, 4, 46261)) and link not in by[46244]
    dig = Store(store.root).circular_digests()["optical:GRB261113A-OT"]
    assert dig["n"] == 5 and dig["optical"] == "detection"
    # a second report does not create another one; the same redshift again changes nothing
    handle_circular({**base, "circularId": 5, "subject": sub, "body": "We detect it at RA, Dec = 10.2003, 20.1."}, store)
    handle_circular({**circ(46261), "circularId": 6}, store)
    assert len(store.transients()) == 1
    assert link in [(ln["type"], ln["id"], ln["how"]) for c in stored_circulars(store.root) if c["id"] == 5
                    for ln in c["links"]]
    n_entries = len([e for e in read_events(store.root)[0] if e["kind"] == "circular"])
    again = next(c for c in stored_circulars(store.root) if c["id"] == 4)
    handle_circular({**base, "circularId": 4, "subject": sub, "body": again["body"]}, store)
    # read again: the same links, nothing written
    assert len([e for e in read_events(store.root)[0] if e["kind"] == "circular"]) == n_entries
    assert [e["kind"] for e in read_events(store.root)[0]].count("optical_update") == 0
    # too far and no GW alert: two detection types, not a candidate
    f = next(f for f in Store(store.root).followups(False) if f["optical"] == "GRB261113A-OT")
    assert f["types"] == ["grb", "optical"] and f["too_far"] and not f["candidate"]
    # a revised redshift, now nearby: the decision is run again and says it is detectable
    handle_circular({**circ(46261), "circularId": 99998, "subject": "GRB 261113A: revised spectroscopic redshift z = 0.03",
                     "body": "A new spectrum gives z = 0.03."}, store)
    assert store.transients()[0].redshift.z == 0.03
    dec = json.loads((store.root / "assess" / "GRB261113A-OT_decision.json").read_text())
    assert dec["detectability"]["worth_search"] and not dec["trigger_targeted_search"]   # not confirmed as a kilonova
    assert len([e for e in read_events(store.root)[0] if e["kind"] == "decision"]) == 2
    f = next(f for f in Store(store.root).followups(False) if f["optical"] == "GRB261113A-OT")
    assert f["candidate"] and not f["too_far"]
    # a redshift measured otherwise is not replaced
    kn.redshift = Redshift(z=0.01, kind="spectroscopic", source="host NGC 1")
    store.save(kn)
    handle_circular({**circ(46261), "circularId": 7, "subject": "GRB 261113A: spectroscopic redshift z = 2.0",
                     "body": "From the spectrum, z = 2.0."}, store)
    assert store.transients()[0].redshift.z == 0.01

    # reported photometry: shown, not in the photometry of the transient until accepted
    kn = store.transients()[0]
    rep = reported_photometry(store.root, kn)
    assert [(p["circular"], p["band"], p["mag"], p["accepted"]) for p in rep if p["circular"] in (1, 4)] == [
        (1, "r", 19.1, False), (4, "r", 19.14, False), (4, "i", 21.0, False)]
    detail = Store(store.root).optical_detail("GRB261113A-OT", False)
    assert len(detail["reported_photometry"]) == len(rep) and detail["optical"]["photometry"] == []
    added = accept_photometry(store, "GRB261113A-OT", [1, 4])
    assert [(p["circular"], p["band"]) for p in added] == [(4, "r")]        # the others have no time
    pt, = store.transients()[0].photometry
    assert (pt.band, pt.mag, pt.mag_err, pt.survey, pt.upper_limit) == ("ps1::r", 19.14, 0.05, "GCN 4", False)
    assert pt.mjd == pytest.approx(61357.6253, abs=1e-3)
    assert [p["accepted"] for p in reported_photometry(store.root, store.transients()[0])
            if p["circular"] == 4] == [True, False]
    assert [p["circular"] for p in accept_photometry(store, "GRB261113A-OT")] == [6, 46261]   # all the others
    with pytest.raises(ValueError):
        accept_photometry(store, "unknown")


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
