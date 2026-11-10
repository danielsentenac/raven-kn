import json
from pathlib import Path

import numpy as np
import pytest

from raven_kn.connection import (ALeRCEConnector, ANTARESConnector, ATLASForcedPhotometryConnector,
                                 KNCandidate, PhotometryPoint, last_nondetection)

DATA = Path(__file__).parent / "data"


class _Resp:
    def __init__(self, payload=None, text=""):
        self._p, self.text = payload, text

    def raise_for_status(self):
        pass

    def json(self):
        return self._p


def test_alerce_parse_real_lightcurve():
    """ALeRCE response for ZTF21abdpqpq (fetched 2026-11-06): 3 detections, 91 non-detections."""
    lc = json.loads((DATA / "alerce_ZTF21abdpqpq_lightcurve.json").read_text())

    class Sess:
        def get(self, url, timeout):
            assert url.endswith("/objects/ZTF21abdpqpq/lightcurve")
            return _Resp(lc)
    kn = ALeRCEConnector(session=Sess()).candidate("ZTF21abdpqpq")
    assert len(kn.detections) == 3 and len(kn.upper_limits) == 91
    assert kn.first_detection.mjd == pytest.approx(59369.3157, abs=1e-3)
    assert {p.band for p in kn.photometry} <= {"g", "r", "i"}
    last = last_nondetection(kn)
    assert last.mjd < kn.first_detection.mjd and last.survey == "ZTF"


def test_atlas_queue_poll_parse():
    table = ("###MJD m dm uJy duJy F err chi/N RA Dec x y maj min phi apfit mag5sig Sky Obs\n"
             "60000.10 19.80 0.05 44.0 2.0 o 0 1.0 1 2 3 4 2 2 0 -0.4 19.5 20 01a\n"
             "59999.10 20.50 0.40 3.0 4.0 c 0 1.0 1 2 3 4 2 2 0 -0.4 19.9 20 01a\n")
    calls = []

    class Sess:
        def post(self, url, headers, data, timeout):
            calls.append(("post", url, headers["Authorization"], data["mjd_min"]))
            return _Resp({"url": "https://x/queue/1/"})

        def get(self, url, headers, timeout):
            calls.append(("get", url))
            if url.endswith("/queue/1/"):
                return _Resp({"finishtimestamp": "now", "result_url": "https://x/result.txt"})
            return _Resp(text=table)
    kn = KNCandidate("t", 10.0, 20.0, photometry=[PhotometryPoint(60000.5, "g", 19.0, 0.1)])
    ATLASForcedPhotometryConnector("TOKEN", session=Sess(), poll_s=0).add_to(kn, days_before=10)
    assert calls[0][2] == "Token TOKEN" and calls[0][3] == pytest.approx(59990.5)
    atl = [p for p in kn.photometry if p.survey == "ATLAS"]
    det = [p for p in atl if not p.upper_limit]
    lim = [p for p in atl if p.upper_limit]
    assert len(det) == 1 and det[0].band == "o" and det[0].mag == pytest.approx(23.9 - 2.5 * 1.6435, abs=0.01)
    assert len(lim) == 1 and lim[0].band == "c" and lim[0].mag == 19.9
    assert last_nondetection(kn).survey == "ATLAS"


def test_atlas_no_data_is_empty():
    """Live reply for AT2017gfo, Aug 2017 (2026-11-09): the task ends with
    error_msg "No data returned" (field too close to the Sun for ATLAS)."""
    class Sess:
        def post(self, url, headers, data, timeout):
            return _Resp({"url": "https://x/queue/2/"})

        def get(self, url, headers, timeout):
            return _Resp({"finishtimestamp": "now", "result_url": None,
                          "error_msg": "No data returned"})
    a = ATLASForcedPhotometryConnector("TOKEN", session=Sess(), poll_s=0)
    assert a.request(197.45, -23.38, 57975.0, 57993.0) == "" and a.parse("") == []


def test_atlas_parse_real_reply():
    """Live reply for SN 2023ixf, MJD 60075-60090 (2026-11-09): limits up to
    hours before discovery (MJD ~60082.7), then o ~ 10.9-11.5 detections."""
    a = ATLASForcedPhotometryConnector("TOKEN")
    pts = a.parse((DATA / "atlas_sn2023ixf.txt").read_text())
    det = [p for p in pts if not p.upper_limit]
    lim = [p for p in pts if p.upper_limit]
    assert len(pts) == 19 and len(det) == 4 and {p.band for p in pts} == {"c", "o"}
    assert all(p.mjd > 60086 and 10.5 < p.mag < 13 for p in det)
    # S/N 3.8 point on the 11-mag supernova: a non-detection, but not a deep one
    bad = [p for p in lim if abs(p.mjd - 60086.3642) < 1e-3][0]
    assert bad.mag < 11
    # just before the first detection ATLAS has a non-detection only to ~9.8 (bad error);
    # the last limit deep enough to matter is hours before discovery
    kn = KNCandidate("SN2023ixf", 210.910674, 54.3117, photometry=pts)
    assert last_nondetection(kn).mag < 11
    last = last_nondetection(kn, depth_mag=19.0)
    assert last.band == "o" and 60082.4 < last.mjd < 60082.5 and 19.1 < last.mag < 19.5


def test_last_nondetection_depth():
    kn = KNCandidate("t", 0, 0, photometry=[PhotometryPoint(10.0, "r", 19.0, 0.1),
                                            PhotometryPoint(9.5, "r", 18.0, upper_limit=True),
                                            PhotometryPoint(9.0, "r", 21.0, upper_limit=True)])
    assert last_nondetection(kn).mjd == 9.5
    assert last_nondetection(kn, depth_mag=20.0).mjd == 9.0


def test_fink_lsst_flux_to_mag_and_limits():
    import json
    from pathlib import Path
    from raven_kn.connection import FinkLSSTConnector
    d = Path(__file__).parent / "data"
    src = json.loads((d / "fink_lsst_sources_170591538829590749.json").read_text())
    fp = json.loads((d / "fink_lsst_fp_170591538829590749.json").read_text())

    class S:
        def get(self, url, params=None, timeout=None):
            class R:
                def raise_for_status(self):
                    pass

                def json(self):
                    return src if url.endswith("sources") else fp
            return R()
    c = FinkLSSTConnector(session=S())
    kn = c.candidate("170591538829590749")
    assert kn.ra == pytest.approx(src[0]["r:ra"]) and kn.photometry
    # one point per visit (forced photometry kept), sorted, TAI -> UTC
    assert len(kn.photometry) == len({r["r:visit"] for r in src} | {r["r:visit"] for r in fp})
    first = min(fp, key=lambda r: r["r:midpointMjdTai"])
    p0 = next(p for p in kn.photometry if abs(p.mjd - (first["r:midpointMjdTai"] - 37 / 86400)) < 1e-6)
    assert p0.mjd == pytest.approx(first["r:midpointMjdTai"] - 37 / 86400)
    assert p0.mag == pytest.approx(31.4 - 2.5 * np.log10(first["r:psfFlux"]))
    assert all(p.survey == "LSST" for p in kn.photometry)
    # a faint forced point becomes a 5-sigma limit
    pt = c._point({"r:midpointMjdTai": 61000.0, "r:band": "g", "r:psfFlux": 300.0, "r:psfFluxErr": 200.0})
    assert pt.upper_limit and pt.mag == pytest.approx(31.4 - 2.5 * np.log10(1000.0))


class _AntaresSess:
    """Recorded ANTARES replies for ZTF21abdpqpq (2026-11-09); the alerts page
    carries a 'next' link although it already holds all 94 alerts."""

    def __init__(self):
        self.calls = []

    def get(self, url, params=None, timeout=None):
        self.calls.append((url, params))
        if url.endswith("/loci"):
            return _Resp(json.loads((DATA / "antares_cone_ZTF21abdpqpq.json").read_text()))
        return _Resp(json.loads((DATA / "antares_alerts_ANT2021owtei_p1.json").read_text()))


def test_antares_matches_alerce():
    sess = _AntaresSess()
    kn = ANTARESConnector(session=sess).candidate(241.3264025, 13.9156065)
    assert kn.name == "ZTF21abdpqpq" and kn.source == "ANTARES"
    dets = [p for p in kn.photometry if not p.upper_limit]
    uls = [p for p in kn.photometry if p.upper_limit]
    assert len(dets) == 3 and len(uls) == 91 and all(p.survey == "ZTF" for p in kn.photometry)
    assert len(sess.calls) == 2                      # stopped at meta.count, no 'next' loop
    q = json.loads(sess.calls[0][1]["elasticsearch_query[locus_listing]"])
    assert "sky_distance" in q["query"]["bool"]["filter"]
    # same points as ALeRCE: merging the two adds nothing
    alerce = ALeRCEConnector.parse(json.loads((DATA / "alerce_ZTF21abdpqpq_lightcurve.json").read_text()))
    a = sorted((round(p.mjd, 4), p.band, p.upper_limit, p.mag) for p in alerce)
    b = sorted((round(p.mjd, 4), p.band, p.upper_limit, p.mag) for p in kn.photometry)
    assert [x[:3] for x in a] == [x[:3] for x in b]
    assert max(abs(x[3] - y[3]) for x, y in zip(a, b)) < 0.02     # float32 storage in ANTARES
    n = len(kn.photometry)
    kn.add_photometry(alerce)
    assert len(kn.photometry) == n


def test_antares_cone_rejects_far_loci():
    """If the server ignored the filter, unrelated loci must not be taken."""
    class Far(_AntaresSess):
        def get(self, url, params=None, timeout=None):
            return _Resp({"data": [{"id": "ANTx", "attributes": {"ra": 37.28, "dec": 9.26,
                                                                  "properties": {}}}]})
    assert ANTARESConnector(session=Far()).cone(241.3264025, 13.9156065) == []
