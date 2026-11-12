"""Data sources of the targeted search: public (GWOSC) and LVK back-ends, routing (offline)."""
import numpy as np
import pytest

pytest.importorskip("pycbc")

from raven_kn.targeted import datasource as D  # noqa: E402
from raven_kn.targeted.datasource import (DataNotAvailable, LVKSource, PublicSource,  # noqa: E402
                                          select_source)

T0 = 1000000000


class FakePublic(PublicSource):
    """GWOSC timelines replaced by fixed segments ("passes" flags, as GWOSC gives them)."""
    def __init__(self, timelines, **kw):
        super().__init__(**kw)
        self.timelines = timelines

    def _timeline(self, flag, start, end):
        return D._clip(self.timelines.get(flag, []), start, end)


def test_segment_arithmetic():
    assert D._subtract([(0, 10)], [(2, 3), (5, 6)]) == [(0, 2), (3, 5), (6, 10)]
    assert D._merge([(5, 7), (0, 2), (2, 3)]) == [(0, 3), (5, 7)]
    assert D._covers([(0, 5), (5, 9)], 1, 8) and not D._covers([(0, 5), (6, 9)], 1, 8)


def test_public_segments_from_timelines():
    p = FakePublic({"H1_CBC_CAT1": [(T0, T0 + 500), (T0 + 600, T0 + 1000)],
                    "H1_CBC_CAT2": [(T0, T0 + 100), (T0 + 102, T0 + 500), (T0 + 600, T0 + 1000)]})
    assert p.analysable("H1", T0, T0 + 1000) == [(T0, T0 + 500), (T0 + 600, T0 + 1000)]
    assert p.vetoed("H1", T0, T0 + 1000) == [(T0 + 100, T0 + 102)]      # CAT2: analysable but flagged
    assert p.available("H1", T0 + 10, T0 + 400) and not p.available("H1", T0 + 400, T0 + 700)
    with pytest.raises(DataNotAvailable):                                # gap: the search drops H1
        p.fetch("H1", T0 + 400, T0 + 700)


def test_public_prefers_cvmfs_then_https(monkeypatch):
    urls = ["pelican://osg-htc.org/gwdata/O2/x/H-H1_GWOSC_O2_4KHZ_R1-1-4096.gwf",
            "osdf:///gwdata/O2/x/H-H1_GWOSC_O2_4KHZ_R1-1-4096.gwf",
            "https://osdf-director.osg-htc.org/gwdata/O2/x/H-H1_GWOSC_O2_4KHZ_R1-1-4096.gwf"]
    assert PublicSource._best(urls).startswith("https://")
    monkeypatch.setattr(D.os.path, "exists", lambda p: p.startswith(D.GWOSC_CVMFS))
    assert PublicSource._best(urls) == "file://" + D.GWOSC_CVMFS + "/gwdata/O2/x/H-H1_GWOSC_O2_4KHZ_R1-1-4096.gwf"


def test_public_frame_type_for_the_run(monkeypatch):
    import gwdatafind
    monkeypatch.setattr(gwdatafind, "find_types", lambda site, host: [
        "H1_AUX_AR1", "H1_GWOSC_O2_16KHZ_R1", "H1_GWOSC_O2_4KHZ_R1", "H1_GWOSC_O3a_4KHZ_R1"])
    monkeypatch.setattr(gwdatafind, "find_times", lambda site, ft, a, b, host: (
        [(T0 - 100, T0 + 5000)] if "O3a" in ft else []))
    assert PublicSource().frame_type("H1", T0, T0 + 64) == "H1_GWOSC_O3a_4KHZ_R1"
    assert PublicSource().frame_type("H1", T0 + 6000, T0 + 6064) is None
    assert PublicSource().channel("H1") == "H1:GWOSC-4KHZ_R1_STRAIN"


ACTIVE = {"H1:DCH-CLEAN_SCIENCE_C02:1": [(T0, T0 + 900)], "L1:DCH-CLEAN_SCIENCE_C02:1": [(T0, T0 + 1000)],
          "V1:ITF_SCIENCE:2": [], "H1:DCH-GLITCH:1": [(T0 + 100, T0 + 102)],
          "L1:DCH-BAD:1": [(T0 + 500, T0 + 510)]}


def test_lvk_segments_science_minus_cat1(tmp_path):
    pytest.importorskip("igwn_ligolw.lsctables")
    from test_segments import write_vdf          # H1 CAT2 (padded 1 s), L1 CAT1
    s = LVKSource.preset("O2", veto_definer=write_vdf(tmp_path / "vdf.xml"), active=ACTIVE)
    assert s.analysable("H1", T0, T0 + 1000) == [(T0, T0 + 900)]
    assert s.vetoed("H1", T0, T0 + 1000) == [(T0 + 99, T0 + 103)]
    assert s.analysable("L1", T0, T0 + 1000) == [(T0, T0 + 500), (T0 + 510, T0 + 1000)]
    assert s.analysable("V1", T0, T0 + 1000) == []
    assert s.channel("V1") == "V1:Hrec_hoft_V1O2Repro2A_16384Hz"
    p = s.provenance()
    assert p["source"] == "lvk" and p["veto_definer"] == "vdf.xml"


def test_routing(monkeypatch):
    lvk = LVKSource.preset("O2", active=ACTIVE)
    pub = FakePublic({"H1_CBC_CAT1": [(T0, T0 + 1000)], "L1_CBC_CAT1": [(T0, T0 + 1000)]})
    none = FakePublic({})
    assert select_source("auto", ["H1", "L1", "V1"], T0, T0 + 1000, pub, lvk) is pub
    monkeypatch.setattr(D, "has_lvk_credentials", lambda: True)
    assert select_source("auto", ["H1", "L1"], T0, T0 + 1000, none, lvk) is lvk
    assert select_source("lvk", ["H1", "L1"], T0, T0 + 1000, pub, lvk) is lvk
    monkeypatch.setattr(D, "has_lvk_credentials", lambda: False)
    with pytest.raises(DataNotAvailable, match="SciToken"):
        select_source("auto", ["H1", "L1"], T0, T0 + 1000, none, lvk)
    with pytest.raises(DataNotAvailable, match="not released yet"):
        select_source("auto", ["H1", "L1"], T0, T0 + 1000, none)
    with pytest.raises(DataNotAvailable):
        select_source("lvk", ["H1", "L1"], T0, T0 + 1000, pub)


def test_search_records_provenance_and_drops_unanalysable_detector():
    """A source wraps the synthetic data: V1 has a science gap, so the search
    runs on H1 L1 and the result says which data and segments were used."""
    from test_targeted import T_INJ, _cfg, synthetic_fetcher
    from raven_kn.targeted.search import TargetedSearchConfig, run
    base = synthetic_fetcher(inject=True)

    class Src(FakePublic):
        def _read(self, ifo, start, end):
            return base(ifo, start, end)

    w = (T_INJ - 400, T_INJ + 400)
    src = Src({"H1_CBC_CAT1": [w], "L1_CBC_CAT1": [w], "V1_CBC_CAT1": [(w[0], T_INJ - 5)],
               "H1_CBC_CAT2": [(w[0], T_INJ - 0.5), (T_INJ + 0.5, w[1])], "L1_CBC_CAT2": [w]},
              sample_rate=_cfg().sample_rate)
    cfg = TargetedSearchConfig(**{**_cfg().__dict__, "ifos": ("H1", "L1", "V1"), "n_slides": 20})
    res = run(cfg, src, progress=False)
    s = res.subsegments[0]
    assert s.ifos == ["H1", "L1"] and "DataNotAvailable" in s.gates["V1"][0]
    assert abs(s.loudest_gps - T_INJ) < 0.05
    assert res.data["source"] == "public" and res.data["loudest_in_cat2"] == ["H1"]
    assert res.data["window"]["V1"]["analysable"] == [(cfg.gps_start, T_INJ - 5)]
    assert np.isfinite(res.p_value)
