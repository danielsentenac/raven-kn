import json
from pathlib import Path

import pytest

from raven_kn.connection import (ALeRCEConnector, ATLASForcedPhotometryConnector, KNCandidate,
                                 PhotometryPoint, last_nondetection)

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


def test_last_nondetection_depth():
    kn = KNCandidate("t", 0, 0, photometry=[PhotometryPoint(10.0, "r", 19.0, 0.1),
                                            PhotometryPoint(9.5, "r", 18.0, upper_limit=True),
                                            PhotometryPoint(9.0, "r", 21.0, upper_limit=True)])
    assert last_nondetection(kn).mjd == 9.5
    assert last_nondetection(kn, depth_mag=20.0).mjd == 9.0
