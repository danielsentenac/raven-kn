"""Host galaxy, redshift and distance connectors on recorded replies (no network)."""
import json
from pathlib import Path

import numpy as np
import pytest

from raven_kn.connection import (DESY6Connector, GladePlusConnector, KNCandidate, NEDConnector,
                                 TNSConnector, find_host, helio_to_cmb, p_chance)
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


ESO243 = (17.6154 + 8 / 3600 / 0.693, -46.0742)     # 8" east of ESO 243-49 (z_spec 0.0224)
EMPTY_NED = "prefname,prefphytype,ra,dec,z,zunc,zflag,n_dist,sep\n"


def des_session():
    return FakeSession({"datalab": Resp((DATA / "desy6_eso243-49.csv").read_text())})


def test_desy6_host_photoz():
    des = DESY6Connector(session=des_session())
    hosts = des.hosts(*ESO243)
    h = hosts[0]
    assert h.catalog == "DES Y6" and h.z_kind == "photometric"
    assert h.z == pytest.approx(0.0224, abs=0.005) and h.z_err >= 0.01
    assert h.sep_arcsec < 9 and h.offset_kpc < 5 and h.p_chance < 1e-2
    assert all(g.offset_kpc <= 50 for g in hosts)
    q = des.session.calls[0][1]["data"]["QUERY"]
    assert "q3c_radial_query" in q and "des_dr2.y6_gold" in q and "ext_mash >= 2" in q


def test_desy6_error_document_raises():
    err = Resp('<?xml version="1.0"?><VOTABLE><RESOURCE><INFO name="QUERY_STATUS" '
               'value="ERROR">ParseException</INFO></RESOURCE></VOTABLE>')
    with pytest.raises(ConnectionError, match="ParseException"):
        DESY6Connector(session=FakeSession({"datalab": err})).fetch(0, 0, 10)


def test_find_host_falls_back_to_des():
    """No NED or GLADE+ host: DES Y6 photo-z; outside the footprint: nothing."""
    ned = NEDConnector(session=FakeSession({"tap/sync": Resp(EMPTY_NED)}))
    glade = GladePlusConnector(session=FakeSession({"asu-tsv": Resp("#\n")}))
    kn = KNCandidate("x", *ESO243)
    log = find_host(kn, ned=ned, glade=glade, des=DESY6Connector(session=des_session()))
    assert kn.redshift.kind == "photometric" and "DES Y6" in kn.redshift.source
    assert any(s.startswith("DES Y6:") for s in log["steps"])
    kn2 = KNCandidate("y", 237.9757, 30.9023)
    log2 = find_host(kn2, ned=ned, glade=glade,
                     des=DESY6Connector(session=FakeSession({"datalab": Resp("ra\n")})))
    assert kn2.redshift is None and "outside the footprint" in " ".join(log2["steps"])
    kn3 = KNCandidate("z", *ESO243)
    find_host(kn3, ned=ned, glade=glade, use_des=False)
    assert kn3.redshift is None


def test_p_chance_half_light_radius():
    assert p_chance(14.0, 7.5, 20.0) > p_chance(14.0, 7.5)
    assert GladePlusConnector.p_chance(18.0, 5.0) == p_chance(18.0, 5.0)
    assert p_chance(None, 5.0) is None


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
