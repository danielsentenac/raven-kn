"""Host galaxy, redshift and distance connectors on recorded replies (no network)."""
import json
from pathlib import Path

import numpy as np
import pytest

from raven_kn.connection import (GladePlusConnector, KNCandidate, NEDConnector, TNSConnector,
                                 find_host, helio_to_cmb)
from raven_kn.treatment.distance_inclination import posterior_for

DATA = Path(__file__).parent / "data"


class Resp:
    def __init__(self, text=None, js=None):
        self.text, self._js = text, js

    def raise_for_status(self):
        pass

    def json(self):
        return self._js


class FakeSession:
    """Routes requests to recorded replies by URL substring."""

    def __init__(self, routes):
        self.routes, self.calls = routes, []

    def _get(self, url, **kw):
        self.calls.append((url, kw))
        for key, resp in self.routes.items():
            if key in url:
                return resp
        raise ConnectionError(url)

    get = post = _get


def ned_session():
    return FakeSession({"tap/sync": Resp((DATA / "ned_cone_ngc4993.csv").read_text()),
                        "nDistance": Resp((DATA / "nedd_ngc4993.html").read_text())})


def test_ned_host_and_nedd_distance():
    ned = NEDConnector(session=ned_session())
    hosts = ned.hosts(197.450374, -23.381495)
    g, n_dist = hosts[0]
    assert g.name == "NGC 4993" and g.z_kind == "spectroscopic" and n_dist > 0
    assert g.z == pytest.approx(helio_to_cmb(0.00987683, g.ra, g.dec), rel=1e-6)
    assert 0.0107 < g.z < 0.0111 and g.offset_kpc < 5
    rows = ned.parse_distances((DATA / "nedd_ngc4993.html").read_text())
    assert any(r["method"].startswith("SBF") and r["d_mpc"] == 40.7 for r in rows)
    d = ned.distance("NGC 4993")
    assert d.mean_mpc == 40.7 and "SBF" in d.method and 1 < d.std_mpc < 2
    assert "Grav" not in d.method


def test_find_host_ngc4993_gives_distance_and_posterior():
    kn = KNCandidate("AT2017gfo", 197.450374, -23.381495)
    log = find_host(kn, ned=NEDConnector(session=ned_session()))
    assert kn.redshift.kind == "spectroscopic" and "NGC 4993" in kn.redshift.source
    assert kn.distance.mean_mpc == 40.7
    s = posterior_for(kn).summary()["luminosity_distance_mpc"]
    assert 37 < s["median"] < 44
    assert any("NED-D" in x for x in log["steps"])


def test_glade_fallback_photometric_host():
    """AT2025ulz: no NED galaxy redshift; GLADE+ host 0.9" away, photometric z."""
    ned = NEDConnector(session=FakeSession({"tap/sync": Resp("prefname,prefphytype,ra,dec,z,zunc,"
                                                             "zflag,n_dist,sep\n")}))
    glade = GladePlusConnector(session=FakeSession(
        {"asu-tsv": Resp((DATA / "gladep_at2025ulz.tsv").read_text())}))
    kn = KNCandidate("AT2025ulz", 237.9757129, 30.9023146)
    find_host(kn, ned=ned, glade=glade)
    assert kn.redshift.kind == "photometric" and kn.redshift.z == pytest.approx(0.102, abs=1e-3)
    assert kn.redshift.z_err > 0.02 and kn.distance is None


def test_glade_mirror_fallback_on_database_error():
    bad = Resp("#INFO\tError=The database is not currently reachable.\t\n")
    good = Resp((DATA / "gladep_at2025ulz.tsv").read_text())
    sess = FakeSession({"cfa.harvard": bad, "idia.ac.za": good})
    g = GladePlusConnector(session=sess)
    assert len(g.parse(g.fetch(237.9757129, 30.9023146, 600))) > 10
    with pytest.raises(ConnectionError):
        GladePlusConnector(session=FakeSession({"cfa.harvard": bad, "idia": bad, "unistra": bad})
                           ).fetch(0, 0, 10)


def test_spectroscopic_redshift_kept():
    kn = KNCandidate("x", 237.9757129, 30.9023146)
    from raven_kn.connection import Redshift
    kn.redshift = Redshift(0.0848, 1e-4, source="Keck spectrum")
    ned = NEDConnector(session=FakeSession({"tap/sync": Resp("prefname,prefphytype,ra,dec,z,zunc,"
                                                             "zflag,n_dist,sep\n")}))
    glade = GladePlusConnector(session=FakeSession(
        {"asu-tsv": Resp((DATA / "gladep_at2025ulz.tsv").read_text())}))
    find_host(kn, ned=ned, glade=glade)
    assert kn.redshift.z == 0.0848 and kn.redshift.source == "Keck spectrum"


TNS_REPLY = {"id_code": 200, "data": {"reply": {
    "objname": "2025ulz", "name_prefix": "AT", "radeg": 237.9757129, "decdeg": 30.9023146,
    "object_type": {"name": "Kilonova", "id": 999}, "redshift": None,
    "hostname": "DESI J237.9757+30.9023", "host_redshift": 0.08484,
    "spectra": [{"jd": 2460906.8, "obsdate": "2025-08-20 07:00:00", "instrument": {"name": "Keck1 - LRIS"},
                 "asciifile": "https://www.wis-tns.org/system/files/tns_2025ulz.txt", "fitsfile": "",
                 "source_group": {"name": "GROWTH"}}]}}}


def test_tns_classification_and_host_redshift():
    sess = FakeSession({"api/get/object": Resp(js=TNS_REPLY)})
    tns = TNSConnector(123, "raven_kn_bot", "KEY", session=sess)
    kn = KNCandidate("AT2025ulz", 237.9757129, 30.9023146)
    info = tns.add_to(kn, "AT2025ulz")
    assert kn.classification == "Kilonova" and info["spectra"][0]["instrument"] == "Keck1 - LRIS"
    assert kn.redshift.kind == "spectroscopic"
    assert kn.redshift.z == pytest.approx(helio_to_cmb(0.08484, kn.ra, kn.dec))
    url, kw = sess.calls[0]
    assert '"tns_id": 123' in kw["headers"]["User-Agent"] and kw["data"]["api_key"] == "KEY"
    assert json.loads(kw["data"]["data"])["objname"] == "2025ulz"
