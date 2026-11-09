import json
from pathlib import Path

import numpy as np
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
