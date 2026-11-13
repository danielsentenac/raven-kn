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


def test_manifest_resolved_beforehand(tmp_path):
    """Segments and frames resolved once (login node), read back without network (job)."""
    from raven_kn.targeted.datasource import ManifestSource, write_manifest

    class Src(FakePublic):
        def frame_urls(self, ifo, start, end):
            return [f"file://localhost/data/{ifo[0]}-{ifo}_X-{T0 - 4096}-4096.gwf",
                    f"file://localhost/data/{ifo[0]}-{ifo}_X-{T0}-4096.gwf"]

    src = Src({"H1_CBC_CAT1": [(T0 - 500, T0 + 500)], "H1_CBC_CAT2": [(T0 - 500, T0 + 10), (T0 + 12, T0 + 500)]})
    man = write_manifest(src, ["H1", "L1"], T0, T0 + 100, tmp_path / "m.json", margin=200)
    assert man["span"] == [T0 - 200, T0 + 300] and len(man["ifos"]["H1"]["frames"]) == 2
    (tmp_path / f"H-H1_X-{T0}-4096.gwf").write_bytes(b"")           # shipped with the job
    m = ManifestSource.load(tmp_path / "m.json", frame_dir=str(tmp_path))
    assert m.name == "public" and m.channel("H1") == "H1:GWOSC-4KHZ_R1_STRAIN"
    assert m.analysable("H1", T0, T0 + 100) == [(T0, T0 + 100)]
    assert m.vetoed("H1", T0, T0 + 100) == [(T0 + 10, T0 + 12)]
    assert m.frame_urls("H1", T0 + 1, T0 + 100) == [f"file://{tmp_path}/H-H1_X-{T0}-4096.gwf"]   # local copy
    assert m.frame_urls("H1", T0 - 10, T0 - 1) == [f"file://localhost/data/H-H1_X-{T0 - 4096}-4096.gwf"]
    assert m.analysable("L1", T0, T0 + 100) == [] and not m.available("L1", T0, T0 + 100)
    with pytest.raises(DataNotAvailable, match="outside the manifest"):
        m.analysable("H1", T0 - 1000, T0)


def test_public_rate_choice():
    assert PublicSource(2048, rate_khz=16).channel("L1") == "L1:GWOSC-16KHZ_R1_STRAIN"
    assert PublicSource(16384).channel("L1") == "L1:GWOSC-16KHZ_R1_STRAIN"


def test_lvk_frames_in_place_or_over_osdf(monkeypatch):
    """A cluster's datafind gives files read in place; datafind.igwn.org gives
    OSDF URLs, downloaded with the SciToken."""
    import gwdatafind
    osdf = ["pelican://osg-htc.org/igwn/ligo/frames/O2/x/H-H1_CLEANED_HOFT_C02-1-4096.gwf",
            "https://osdf-director.osg-htc.org/igwn/ligo/frames/O2/x/H-H1_CLEANED_HOFT_C02-1-4096.gwf"]
    monkeypatch.setattr(gwdatafind, "find_urls", lambda *a, **k: osdf)
    s = LVKSource.preset("O2", datafind_host="datafind.igwn.org")
    assert s.frame_urls("H1", T0, T0 + 10) == [osdf[1]]
    monkeypatch.setattr(gwdatafind, "find_urls", lambda *a, **k: osdf + ["file://localhost/ceph/x/H-H1_CLEANED_HOFT_C02-1-4096.gwf"])
    assert s.frame_urls("H1", T0, T0 + 10) == ["file://localhost/ceph/x/H-H1_CLEANED_HOFT_C02-1-4096.gwf"]
    assert s._local_path("file://localhost/ceph/x/a.gwf") == "/ceph/x/a.gwf"
    monkeypatch.setenv("BEARER_TOKEN", "abc")
    assert D.bearer_token() == "abc" and D.has_lvk_credentials() and s._token() == "abc"
    assert PublicSource()._token() is None


def test_strain_from_raw_front_end_channel():
    """The RAW channel is the arm length whitened by 6 zeros at 0.3 Hz and 6
    poles at 30 Hz; de-whitening recovers the strain in band."""
    from pycbc.types import TimeSeries
    from raven_kn.targeted.datasource import strain_from_deltal
    fs, n = 4096, 4096 * 16
    t = np.arange(n) / fs
    h = 1e-21 * (np.sin(2 * np.pi * 50 * t) + 0.5 * np.sin(2 * np.pi * 300 * t + 1.0))
    f = np.fft.rfftfreq(n, 1 / fs)
    W = ((1 + 1j * f / 0.3) / (1 + 1j * f / 30.0)) ** 6
    raw = np.fft.irfft(np.fft.rfft(h) * W * 3994.5, n)          # what the front end records
    out = strain_from_deltal(TimeSeries(raw, delta_t=1 / fs, epoch=T0)).numpy()
    mid = slice(2 * fs, n - 2 * fs)                               # away from the tapered ends
    assert np.max(np.abs(out[mid] - h[mid])) < 1e-3 * 1e-21
    assert float(strain_from_deltal(TimeSeries(raw, delta_t=1 / fs, epoch=T0)).start_time) == T0
    raw_src = LVKSource.preset("O2-raw")
    assert raw_src.channel("L1") == "L1:CAL-DELTAL_EXTERNAL_DQ" and "V1" not in raw_src.frame_types


def test_manifest_extracts_strain_into_small_frames(tmp_path):
    from pycbc.types import TimeSeries
    from raven_kn.targeted.datasource import ManifestSource, write_manifest

    class Src(FakePublic):
        def frame_urls(self, ifo, start, end):
            return ["file://localhost/archive/huge.gwf"]

        def _read(self, ifo, start, end):
            n = int(round((end - start) * 256))
            return TimeSeries(np.arange(n, dtype=np.float32), delta_t=1 / 256, epoch=start)

    src = Src({"H1_CBC_CAT1": [(T0 - 50.5, T0 + 300)]})
    man = write_manifest(src, ["H1"], T0, T0 + 10, tmp_path / "m.json", margin=100,
                         extract_dir=str(tmp_path / "frames"))
    m = man["ifos"]["H1"]
    assert m["channel"] == "H1:RAVENKN_STRAIN" and m["original_frames"] == ["file://localhost/archive/huge.gwf"]
    assert m["analysable"] == [(T0 - 50, T0 + 110)]
    ms = ManifestSource.load(tmp_path / "m.json", sample_rate=256)
    ts = ms.fetch("H1", T0, T0 + 10)
    assert len(ts) == 2560 and ts.dtype == "float64" and float(ts[0]) == 50 * 256


def test_extract_elsewhere_and_add_to_manifest(tmp_path):
    """A detector whose archive is at another site: its strain is cut out there
    and added to the manifest made where the search runs."""
    from pycbc.frame import write_frame
    from pycbc.types import TimeSeries
    from raven_kn.targeted.datasource import ManifestSource, add_to_manifest, extract_strain, write_manifest
    big = tmp_path / "V-V1O2Repro2A-1000000000-64.gwf"
    write_frame(str(big), "V1:Hrec_hoft", TimeSeries(np.arange(64 * 256, dtype=np.float32), delta_t=1 / 256, epoch=T0))
    entry = extract_strain([str(big)], "V1:Hrec_hoft", "V1", T0 + 10.5, T0 + 50, tmp_path / "out")
    assert entry["extracted"] == [T0 + 11, T0 + 50] and entry["original_frames"] == [big.name]
    write_manifest(FakePublic({}), [], T0 + 20, T0 + 30, tmp_path / "m.json", margin=20)
    add_to_manifest(tmp_path / "m.json", "V1", entry, [(T0, T0 + 45)], [(T0 + 40, T0 + 41)], note="h(t) from tape")
    m = ManifestSource.load(tmp_path / "m.json", sample_rate=256)
    assert m.analysable("V1", T0, T0 + 40) == [(T0 + 11, T0 + 40)] and m.vetoed("V1", T0 + 30, T0 + 40) == []
    ts = m.fetch("V1", T0 + 20, T0 + 30)
    assert len(ts) == 2560 and float(ts[0]) == 20 * 256 and ts.dtype == "float64"
    assert m.provenance()["original"]["V1"]["original_channel"] == "V1:Hrec_hoft"
    assert m.provenance()["notes"]["V1"] == "h(t) from tape"
