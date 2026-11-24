"""Brokers at the position of a transient: answers in the shape the live
services gave on 2026-11-24 (ZTF21abasjty; long lists cut)."""
import json

import pytest

from raven_kn.brokers import (BROKERS, LasairConnector, _alerce, _fink_ztf, crossmatch,
                              crossmatch_transient, stored)
from raven_kn.connection import OpticalTransient
from raven_kn.serve import Store

OBJECT = {
    "objectId": "ZTF21abasjty",
    "objectData": {"ncand": 16, "ramean": 194.516899, "decmean": 47.980826, "discMjd": 60990.4940394,
                   "latestMjd": 61030.5355787, "peakMag": "19.96±0.20", "peakFilter": "g"},
    "candidates": [
        {"candid": 3276535574315015012, "mjd": 61030.5355787, "fid": 2, "magpsf": 20.356, "sigmapsf": 0.19,
         "isdiffpos": "t"},
        {"mjd": 61029.5, "fid": 1, "diffmaglim": 20.9, "magpsf": 20.9},
        {"candid": 1, "mjd": 61028.5, "fid": 1, "magpsf": 20.1, "sigmapsf": 0.1, "isdiffpos": "f"}],
    "sherlock": {"classification": "SN", "catalogue_table_name": "NED/SDSS", "catalogue_object_id": "NGC 4993",
                 "catalogue_object_type": "galaxy", "separationArcsec": 10.3, "physical_separation_kpc": 2.1,
                 "direct_distance": 0.0, "distance": 41.0, "z": 0.0098, "photoZ": 0.0, "photoZErr": 0.0,
                 "description": 'The transient is possibly associated with <em><a href="https://x">NGC 4993</a></em>; '
                                'a J=9.3 mag galaxy.'},
    "TNS": {"tns_name": "2017gfo", "tns_prefix": "AT", "type": "", "z": None},
}
SHERLOCK = {"classifications": {"transient_00000": ["SN", "The transient is <b>possibly</b> associated."]},
            "crossmatches": [{"catalogue_object_id": "far", "rank": 2, "separationArcsec": 9.0},
                             {"catalogue_object_id": "SDSSJ125758.11+475929.8", "rank": 1, "z": None,
                              "photoZ": 0.157, "photoZErr": 0.136, "separationArcsec": 4.26,
                              "catalogue_table_name": "NED/SDSS/PS1", "catalogue_object_type": "galaxy"}]}


class Resp:
    def __init__(self, data, status=200):
        self.data, self.status = data, status

    def raise_for_status(self):
        if self.status != 200:
            raise RuntimeError(f"{self.status} Server Error for url: https://x/?token=secret")

    def json(self):
        return self.data


class Session:
    """Answers by the end of the address; records what was asked."""

    def __init__(self, answers):
        self.answers, self.asked = answers, []

    def get(self, url, params=None, headers=None, timeout=None, json=None):
        self.asked.append((url, params or json, headers))
        return next(v for k, v in self.answers.items() if url.rstrip("/").endswith(k))

    post = get


def test_lasair_object_with_host_and_tns():
    s = Session({"cone": Resp([{"object": "ZTF2", "separation": 2.5}, {"object": "ZTF21abasjty", "separation": 0.01}]),
                 "object": Resp(OBJECT)})
    c = LasairConnector(token="tok", session=s)
    m = c.match(194.5169, 47.98083)
    assert (m["status"], m["object"], m["sep_arcsec"], m["others"]) == ("match", "ZTF21abasjty", 0.01, ["ZTF2"])
    assert m["classification"] == "SN" and "supernova" in m["classification_text"]
    assert m["host"] == {"name": "NGC 4993", "catalogue": "NED/SDSS", "type": "galaxy", "z": 0.0098,
                         "photo_z": None, "photo_z_err": None, "sep_arcsec": 10.3, "sep_kpc": 2.1,
                         "distance_mpc": 41.0}
    assert m["description"] == "The transient is possibly associated with NGC 4993; a J=9.3 mag galaxy."
    assert m["tns"] == {"name": "AT 2017gfo", "type": None, "z": None}
    assert m["first_detection"].startswith("2025-11-11T11:51") and m["url"].endswith("/objects/ZTF21abasjty/")
    assert s.asked[0][2] == {"Authorization": "Token tok"} and s.asked[0][1]["radius"] == 3.0
    pts = c.parse(OBJECT)                      # the negative difference is left out
    assert [(p.band, p.mag, p.upper_limit) for p in pts] == [("g", 20.9, True), ("r", 20.356, False)]


def test_lasair_without_object_gives_the_catalogue_context():
    s = Session({"cone": Resp([]), "sherlock/position": Resp(SHERLOCK)})
    m = LasairConnector(token="tok", session=s).match(1.0, 2.0)
    assert m["status"] == "none" and m["classification"] == "SN"
    assert m["host"]["name"] == "SDSSJ125758.11+475929.8" and m["host"]["photo_z"] == 0.157
    assert m["host"]["z"] is None and m["description"] == "The transient is possibly associated."
    assert LasairConnector(token="").match(1.0, 2.0)["status"] == "not_configured"


def test_fink_and_alerce():
    f = _fink_ztf(194.5169, 47.98083, 3.0, Session({"conesearch": Resp([
        {"d:classification": "QSO", "d:nalerthist": 16, "i:dec": 47.98089, "i:ra": 194.51697, "i:objectId": "ZTF21abasjty",
         "i:jd": 2461031.0355787, "i:jdstarthist": 2458654.69, "v:separation_degree": 0.0000746697}])}))
    assert (f["object"], f["classification"], f["n_detections"]) == ("ZTF21abasjty", "QSO", 16)
    assert f["sep_arcsec"] == pytest.approx(0.2688, abs=1e-3) and "first_detection" not in f
    a = _alerce(194.5169, 47.98083, 3.0, Session({"objects": Resp({"items": [
        {"oid": "ZTF21abasjty", "ndet": 19, "firstmjd": 60990.494, "lastmjd": 61030.536, "meanra": 194.51689,
         "meandec": 47.98083, "class": "SNIa", "classifier": "lc_classifier", "probability": 0.61}]})}))
    assert (a["object"], a["classification"], a["classification_text"]) == \
        ("ZTF21abasjty", "SNIa", "lc_classifier, probability 0.61")
    assert _alerce(1.0, 2.0, 3.0, Session({"objects": Resp({"items": []})})) == {"status": "none"}


def test_a_broker_that_fails_does_not_stop_the_others(tmp_path):
    def down(ra, dec, radius):
        raise RuntimeError("500 Server Error for url: https://x/?token=secret")

    brokers = {"a": ("A", lambda ra, dec, radius: {"status": "match", "object": "ZTF1", "sep_arcsec": 0.4}),
               "b": ("B", down)}
    r = crossmatch(10.0, 20.0, brokers=brokers)
    assert [(b["broker"], b["status"]) for b in r["brokers"]] == [("a", "match"), ("b", "error")]
    assert "secret" not in json.dumps(r) and r["brokers"][1]["note"] == "RuntimeError: 500 Server Error"
    assert set(BROKERS) == {"lasair", "fink_ztf", "fink_lsst", "alerce", "antares"}

    kn = OpticalTransient("KN-1", 10.0, 20.0)
    first = crossmatch_transient(kn, tmp_path, brokers=brokers)
    assert stored(tmp_path, "KN-1") == first and first["optical"] == "KN-1"
    assert crossmatch_transient(kn, tmp_path, max_age_h=12, brokers={"a": ("A", down)}) == first   # young: kept
    # the broker that answered before and fails now keeps its last answer, marked as old
    again = crossmatch_transient(kn, tmp_path, brokers={"a": ("A", down), "b": ("B", down)})
    a, b = again["brokers"]
    assert (a["status"], a["object"], a["stale_since"]) == ("match", "ZTF1", first["time"])
    assert b["status"] == "error" and "stale_since" not in b
    (tmp_path / "optical").mkdir()
    (tmp_path / "optical" / "KN-1.json").write_text(json.dumps({"name": "KN-1", "ra": 10.0, "dec": 20.0}))
    assert Store(tmp_path).optical_detail("KN-1", False)["brokers"]["brokers"][0]["object"] == "ZTF1"
