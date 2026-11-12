"""Where the strain of a targeted search comes from: two back-ends behind one
interface.

``public`` (:class:`PublicSource`, no credentials)
    Data released by GWOSC: catalogue events and whole runs after the
    proprietary period. Frames are found with datafind
    (``datafind.gwosc.org``, the server behind ``gw_data_find -r
    datafind.gwosc.org --show-types``) and read from CVMFS where mounted, else
    downloaded over OSDF, else from gwosc.org's archive; segments come from
    the GWOSC timelines
    (``H1_DATA``, ``H1_CBC_CAT1``, ``H1_CBC_CAT2``). DQSegDB is not needed.

``lvk`` (:class:`LVKSource`, IGWN SciToken)
    The collaboration's calibrated h(t) frames, available at once. Science
    segments and veto-definer categories come from DQSegDB
    (:mod:`raven_kn.targeted.segments`), frames from an IGWN datafind server
    (on a cluster: its local one, files read in place).

Both give the search the same two things: the analysable time of a detector
(science data passing CAT1) and its strain. A source is a
:data:`~raven_kn.targeted.search.StrainFetcher` (callable ``(ifo, start,
end)``); it refuses stretches that are not analysable, so the search drops
that detector for the sub-segment. :func:`select_source` picks the back-end:
public data when it covers the window, else LVK data when credentials are
present, else :class:`DataNotAvailable`.

Every result records what was used (:meth:`DataSource.provenance`): the same
event gives different numbers on different frames (GW170817: L1 from GWOSC
or from the glitch-cleaned frame).
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

GWOSC_DATAFIND = "datafind.gwosc.org"
GWOSC_CVMFS = "/cvmfs/gwosc.osgstorage.org"
DQSEGDB = "https://segments.igwn.org"

Segs = list[tuple[float, float]]


class DataNotAvailable(RuntimeError):
    """No analysable data from this source for the requested detector and time."""


def _covers(segs: Segs, start: float, end: float) -> bool:
    return any(a <= start and b >= end for a, b in _merge(segs))


def _merge(segs: Segs) -> Segs:
    out: Segs = []
    for a, b in sorted((float(a), float(b)) for a, b in segs):
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def _subtract(segs: Segs, minus: Segs) -> Segs:
    out = []
    for a, b in _merge(segs):
        for c, d in _merge(minus):
            if d <= a or c >= b:
                continue
            if c > a:
                out.append((a, c))
            a = max(a, d)
        if a < b:
            out.append((a, b))
    return out


def _clip(segs: Segs, start: float, end: float) -> Segs:
    return [(max(a, start), min(b, end)) for a, b in _merge(segs) if b > start and a < end]


@dataclass
class DataSource:
    """Common part of the two back-ends. Sub-classes give
    :meth:`analysable`, :meth:`vetoed`, :meth:`frame_urls`, :meth:`channel`."""
    sample_rate: int = 2048
    name: str = ""
    used: dict = field(default_factory=dict, repr=False)   # ifo -> what was read (provenance)

    # -- to be provided -----------------------------------------------------
    def analysable(self, ifo: str, start: float, end: float) -> Segs:
        """Science time passing CAT1 in [start, end]."""
        raise NotImplementedError

    def vetoed(self, ifo: str, start: float, end: float) -> Segs:
        """CAT2 time in [start, end]: analysable, but triggers there are suspect."""
        raise NotImplementedError

    def frame_urls(self, ifo: str, start: float, end: float) -> list[str]:
        raise NotImplementedError

    def channel(self, ifo: str) -> str:
        raise NotImplementedError

    # -- shared -------------------------------------------------------------
    def available(self, ifo: str, start: float, end: float) -> bool:
        try:
            return _covers(self.analysable(ifo, start, end), start, end)
        except Exception:
            return False

    cache_dir: Optional[str] = None          # downloaded frames (default ~/.cache/raven-kn/frames)

    def _token(self) -> Optional[str]:
        return None

    def _local_path(self, url: str) -> str:
        if url.startswith("file://"):
            return url[len("file://"):].replace("localhost", "", 1)
        if url.startswith("https://"):       # OSDF: the director redirects to a cache (port 8443)
            import requests
            d = Path(self.cache_dir or Path.home() / ".cache" / "raven-kn" / "frames")
            d.mkdir(parents=True, exist_ok=True)
            dest = d / url.rsplit("/", 1)[-1]
            if not dest.exists():
                tmp = dest.with_suffix(".part")
                tok = self._token()
                try:
                    with requests.get(url, stream=True, timeout=(15, 600),
                                      headers={"Authorization": f"Bearer {tok}"} if tok else {}) as r:
                        r.raise_for_status()
                        with open(tmp, "wb") as f:
                            for chunk in r.iter_content(1 << 20):
                                f.write(chunk)
                except requests.RequestException as e:
                    tmp.unlink(missing_ok=True)
                    raise DataNotAvailable(f"OSDF download failed ({type(e).__name__}): {url}") from e
                tmp.rename(dest)
            return str(dest)
        raise DataNotAvailable(f"cannot read {url}")

    def fetch(self, ifo: str, start: float, end: float):
        """Strain of ``ifo`` over [start, end] as a PyCBC TimeSeries at
        ``sample_rate``. Raises :class:`DataNotAvailable` where the stretch
        is not entirely analysable."""
        ok = self.analysable(ifo, start, end)
        if not _covers(ok, start, end):
            raise DataNotAvailable(f"{ifo} not analysable over [{start:.0f}, {end:.0f}] "
                                   f"({self.name}): analysable {ok}")
        ts = self._read(ifo, start, end)
        if ts.dtype != "float64":            # Virgo's O2 h(t) is single precision
            ts = ts.astype("float64")
        from pycbc.filter import resample_to_delta_t
        if abs(ts.delta_t - 1.0 / self.sample_rate) > 1e-12:
            ts = resample_to_delta_t(ts, 1.0 / self.sample_rate)
        return ts

    def _read(self, ifo: str, start: float, end: float):
        from pycbc.frame import read_frame
        urls = self.frame_urls(ifo, start, end)
        if not urls:
            raise DataNotAvailable(f"no {self.name} frames for {ifo} over [{start:.0f}, {end:.0f}]")
        paths = [self._local_path(u) for u in urls]
        u = self.used.setdefault(ifo, {"channel": self.channel(ifo), "frames": []})
        u["frames"] = sorted(set(u["frames"]) | {os.path.basename(p) for p in paths})
        return read_frame(paths, self.channel(ifo), start_time=start, end_time=end)

    __call__ = fetch

    def provenance(self) -> dict:
        return {"source": self.name, "sample_rate": self.sample_rate, "read": self.used}


# --------------------------------------------------------------------------
# Public data (GWOSC)
# --------------------------------------------------------------------------

@dataclass
class PublicSource(DataSource):
    """GWOSC open data: datafind + OSDF/CVMFS frames, GWOSC timelines."""
    name: str = "public"
    host: str = GWOSC_DATAFIND
    ca_bundle: Optional[str] = None
    rate_khz: Optional[int] = None           # 4 or 16: GWOSC frame types <IFO>_GWOSC_<run>_4KHZ_R1 or
                                             # _16KHZ_R1 (default: 4, or 16 if the search samples faster)
    snippets: bool = True                    # fall back to gwosc.org's archive when frames cannot be read
    _types: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        if self.ca_bundle:
            os.environ.setdefault("REQUESTS_CA_BUNDLE", self.ca_bundle)
            os.environ.setdefault("SSL_CERT_FILE", self.ca_bundle)

    def _timeline(self, flag: str, start: float, end: float) -> Segs:
        from gwosc.timeline import get_segments
        return [(float(a), float(b)) for a, b in get_segments(flag, int(start), int(-(-end // 1)))]

    def analysable(self, ifo, start, end):
        return _clip(self._timeline(f"{ifo}_CBC_CAT1", start, end), start, end)

    def vetoed(self, ifo, start, end):
        return _subtract(self.analysable(ifo, start, end),
                         self._timeline(f"{ifo}_CBC_CAT2", start, end))

    def frame_type(self, ifo: str, start: float, end: float) -> Optional[str]:
        """The run's GWOSC frame type covering [start, end] (4 kHz unless the
        search samples faster), e.g. ``H1_GWOSC_O2_4KHZ_R1``."""
        import gwdatafind
        site = ifo[0]
        if site not in self._types:
            self._types[site] = sorted(gwdatafind.find_types(site=site, host=self.host))
        rate = self._rate()
        for ft in self._types[site]:
            if not re.fullmatch(rf"{ifo}_GWOSC_\w+_{rate}_R\d+", ft):
                continue
            if _covers([(float(a), float(b)) for a, b in
                        gwdatafind.find_times(site, ft, int(start), int(-(-end // 1)), host=self.host)],
                       start, end):
                return ft
        return None

    def _rate(self) -> str:
        k = self.rate_khz or (16 if self.sample_rate > 4096 else 4)
        return f"{k}KHZ"

    def channel(self, ifo):
        return f"{ifo}:GWOSC-{self._rate()}_R1_STRAIN"

    def frame_urls(self, ifo, start, end):
        import gwdatafind
        ft = self.frame_type(ifo, start, end)
        if ft is None:
            return []
        urls = gwdatafind.find_urls(ifo[0], ft, int(start), int(-(-end // 1)), host=self.host,
                                    urltype=None)
        by_file: dict[str, list[str]] = {}
        for u in urls:
            by_file.setdefault(u.rsplit("/", 1)[-1], []).append(u)
        self.used.setdefault(ifo, {"channel": self.channel(ifo), "frames": []})["frame_type"] = ft
        return [self._best(v) for _, v in sorted(by_file.items())]

    @staticmethod
    def _best(urls: Sequence[str]) -> str:
        """CVMFS copy when mounted, else the https (OSDF director) URL."""
        for u in urls:
            path = re.sub(r"^(osdf|pelican)://[^/]*", "", u) if u.startswith(("osdf", "pelican")) else None
            if path and os.path.exists(GWOSC_CVMFS + path):
                return "file://" + GWOSC_CVMFS + path
        for scheme in ("file://", "https://"):
            for u in urls:
                if u.startswith(scheme):
                    return u
        return urls[0]

    def _read(self, ifo, start, end):
        try:
            return super()._read(ifo, start, end)
        except DataNotAvailable:
            if not self.snippets:
                raise
        # gwosc.org's own file archive (https, port 443): the same released strain, for
        # machines that reach neither CVMFS nor the OSDF caches
        from gwpy.timeseries import TimeSeries
        ts = TimeSeries.fetch_open_data(ifo, start, end, cache=True)
        self.used[ifo] = {"frame_type": self.used.get(ifo, {}).get("frame_type"),
                          "read_from": "gwosc.org archive (gwpy fetch_open_data)"}
        return ts.to_pycbc()

    def provenance(self):
        return {**super().provenance(), "datafind": self.host,
                "segments": "GWOSC timelines (<IFO>_CBC_CAT1; CAT2 reported)"}


# --------------------------------------------------------------------------
# LVK data (proprietary)
# --------------------------------------------------------------------------

# Frame types, channels and science flags that were checked on real data.
LVK_PRESETS = {
    # O2 final calibration, cleaned h(t) (LVK's O2 PyGRB configuration); CIT, 2026-11
    "O2": {
        "frame_types": {"H1": "H1_CLEANED_HOFT_C02", "L1": "L1_CLEANED_HOFT_C02", "V1": "V1O2Repro2A"},
        "channels": {"H1": "H1:DCH-CLEAN_STRAIN_C02", "L1": "L1:DCH-CLEAN_STRAIN_C02",
                     "V1": "V1:Hrec_hoft_V1O2Repro2A_16384Hz"},
        "science_flags": {"H1": "H1:DCH-CLEAN_SCIENCE_C02:1", "L1": "L1:DCH-CLEAN_SCIENCE_C02:1",
                          "V1": "V1:ITF_SCIENCE:2"},
    },
}


@dataclass
class LVKSource(DataSource):
    """Collaboration h(t): DQSegDB segments and IGWN datafind. ``frame_types``,
    ``channels`` and ``science_flags`` are per detector (see
    :data:`LVK_PRESETS`); ``veto_definer`` is the CBC veto-definer XML of the
    run. ``datafind_host``: None uses ``$GWDATAFIND_SERVER`` (on a cluster,
    its local server, which returns files readable in place)."""
    name: str = "lvk"
    frame_types: dict = field(default_factory=dict)
    channels: dict = field(default_factory=dict)
    science_flags: dict = field(default_factory=dict)
    veto_definer: Optional[str] = None
    datafind_host: Optional[str] = None
    segment_url: str = DQSEGDB
    active: Optional[dict] = field(default=None, repr=False)   # {flag: segments} instead of DQSegDB (tests)

    @classmethod
    def preset(cls, run: str, **kw) -> "LVKSource":
        if run not in LVK_PRESETS:
            raise KeyError(f"no LVK preset for {run!r} (have {sorted(LVK_PRESETS)}): "
                           "give frame_types, channels and science_flags")
        return cls(**{**LVK_PRESETS[run], **kw})

    def _science(self, ifo, start, end) -> Segs:
        flag = self.science_flags[ifo]
        if self.active is not None:
            return _clip(self.active.get(flag, []), start, end)
        from gwpy.segments import DataQualityFlag
        q0, q1 = int(start // 1), int(-(-end // 1))          # DQSegDB answers on whole seconds
        return _clip([(float(a), float(b)) for a, b in
                      DataQualityFlag.query(flag, q0, q1, url=self.segment_url).active], start, end)

    def _cat(self, ifo, start, end, cat) -> Segs:
        if not self.veto_definer:
            return []
        from .segments import veto_definer_segments
        return _clip(veto_definer_segments(self.veto_definer, [ifo], int(start // 1), int(-(-end // 1)),
                                           categories=(cat,), url=self.segment_url, active=self.active),
                     start, end)

    def analysable(self, ifo, start, end):
        if ifo not in self.science_flags:
            return []
        return _subtract(self._science(ifo, start, end), self._cat(ifo, start, end, 1))

    def vetoed(self, ifo, start, end):
        return _clip(self._cat(ifo, start, end, 2), start, end)

    def channel(self, ifo):
        return self.channels[ifo]

    def frame_urls(self, ifo, start, end):
        import gwdatafind
        kw = {"host": self.datafind_host} if self.datafind_host else {}
        self.used.setdefault(ifo, {"channel": self.channel(ifo), "frames": []})["frame_type"] = \
            self.frame_types[ifo]
        urls = gwdatafind.find_urls(ifo[0], self.frame_types[ifo], int(start), int(-(-end // 1)),
                                    urltype=None, **kw)
        by_file: dict[str, list[str]] = {}
        for u in urls:
            by_file.setdefault(u.rsplit("/", 1)[-1], []).append(u)
        # a cluster's datafind gives files readable in place; datafind.igwn.org gives OSDF
        # URLs, downloaded with the SciToken (scopes read:/ligo, read:/virgo, ...)
        return [next((u for sch in ("file://", "https://") for u in v if u.startswith(sch)), v[0])
                for _, v in sorted(by_file.items())]

    def _token(self):
        return bearer_token()

    def provenance(self):
        return {**super().provenance(), "datafind": self.datafind_host or os.environ.get("GWDATAFIND_SERVER"),
                "segments": self.segment_url, "science_flags": self.science_flags,
                "veto_definer": os.path.basename(self.veto_definer) if self.veto_definer else None}


def write_manifest(source: DataSource, ifos: Sequence[str], start: float, end: float,
                   path: str | Path, margin: float = 256.0) -> dict:
    """Resolve segments and frame files of ``source`` for [start - margin,
    end + margin] and write them to a JSON file. Done where credentials and
    datafind are (a login node); a job without either then reads the data
    with :class:`ManifestSource`. ``margin`` covers the data the search needs
    before and after the window (template length, padding)."""
    import json
    a, b = start - margin, end + margin
    man = {"source": source.name, "span": [a, b], "ifos": {}, "provenance": None}
    for ifo in ifos:
        try:
            man["ifos"][ifo] = {"analysable": source.analysable(ifo, a, b),
                                "cat2": source.vetoed(ifo, a, b),
                                "channel": source.channel(ifo),
                                "frames": source.frame_urls(ifo, a, b)}
        except Exception as e:
            man["ifos"][ifo] = {"error": f"{type(e).__name__}: {e}"}
    man["provenance"] = source.provenance()
    Path(path).write_text(json.dumps(man, indent=1))
    return man


@dataclass
class ManifestSource(DataSource):
    """A source resolved beforehand by :func:`write_manifest`: no network, no
    credentials. Frames are read from ``frame_dir`` when they are there (files
    shipped with a job), else from their original paths."""
    manifest: dict = field(default_factory=dict)
    frame_dir: Optional[str] = None

    @classmethod
    def load(cls, path: str | Path, sample_rate: int = 2048, frame_dir: Optional[str] = None
             ) -> "ManifestSource":
        import json
        man = json.loads(Path(path).read_text())
        return cls(sample_rate=sample_rate, name=man["source"], manifest=man, frame_dir=frame_dir)

    def _ifo(self, ifo, start, end) -> dict:
        m = self.manifest["ifos"].get(ifo)
        a, b = self.manifest["span"]
        if not m or "error" in m:
            raise DataNotAvailable(f"{ifo} not in the manifest" + (f": {m['error']}" if m else ""))
        if start < a or end > b:
            raise DataNotAvailable(f"[{start:.0f}, {end:.0f}] outside the manifest's span "
                                   f"[{a:.0f}, {b:.0f}]")
        return m

    def analysable(self, ifo, start, end):
        return _clip([tuple(s) for s in self._ifo(ifo, start, end)["analysable"]], start, end)

    def vetoed(self, ifo, start, end):
        return _clip([tuple(s) for s in self._ifo(ifo, start, end)["cat2"]], start, end)

    def channel(self, ifo):
        return self.manifest["ifos"][ifo]["channel"]

    def frame_urls(self, ifo, start, end):
        out = []
        for u in self._ifo(ifo, start, end)["frames"]:
            name = u.rsplit("/", 1)[-1]
            m = re.search(r"-(\d+)-(\d+)\.gwf$", name)
            if m and (int(m[1]) >= end or int(m[1]) + int(m[2]) <= start):
                continue
            local = os.path.join(self.frame_dir or ".", name)
            out.append("file://" + os.path.abspath(local) if os.path.exists(local) else u)
        return out

    def provenance(self):
        return {**(self.manifest.get("provenance") or {}), "sample_rate": self.sample_rate,
                "read": self.used, "resolved_beforehand": True}


def bearer_token() -> Optional[str]:
    """The IGWN SciToken, where the IGWN tools look for one. ``htgettoken -a
    vault.ligo.org -i igwn`` writes one valid for DQSegDB (dqsegdb.read),
    datafind (gwdatafind.read) and the frames (read:/ligo, read:/virgo,
    read:/kagra); a token for one service only can be asked with
    ``--audience https://datafind.igwn.org --scopes gwdatafind.read``."""
    for var in ("BEARER_TOKEN", "SCITOKEN"):
        if os.environ.get(var):
            return os.environ[var].strip()
    uid = os.getuid()
    for c in (os.environ.get("BEARER_TOKEN_FILE"), os.environ.get("SCITOKEN_FILE"),
              os.path.join(os.environ.get("XDG_RUNTIME_DIR", f"/run/user/{uid}"), f"bt_u{uid}"),
              f"/tmp/bt_u{uid}"):
        if c and os.path.isfile(c):
            return Path(c).read_text().strip()
    return None


def has_lvk_credentials() -> bool:
    return bearer_token() is not None


def select_source(mode: str, ifos: Sequence[str], start: float, end: float,
                  public: Optional[PublicSource] = None, lvk: Optional[LVKSource] = None,
                  min_ifos: int = 2) -> DataSource:
    """The source for a search over [start, end].

    ``mode``: "public", "lvk", or "auto": public data if at least ``min_ifos``
    detectors are analysable somewhere in the window; else ``lvk`` if it is
    configured and a SciToken is present; else :class:`DataNotAvailable`
    saying what is missing."""
    if mode == "public":
        return public or PublicSource()
    if mode == "lvk":
        if lvk is None:
            raise DataNotAvailable("LVK data asked for but not configured (frame types, channels, "
                                   "science flags)")
        return lvk
    if mode != "auto":
        raise ValueError(f"unknown data source {mode!r}")
    public = public or PublicSource()
    n_pub, why = 0, ""
    try:
        n_pub = sum(bool(public.analysable(i, start, end)) for i in ifos)
    except Exception as e:                            # GWOSC unreachable
        why = f" (GWOSC query failed: {type(e).__name__})"
    if n_pub >= min_ifos:
        return public
    if lvk is not None and has_lvk_credentials():
        return lvk
    raise DataNotAvailable(
        f"no public data for [{start:.0f}, {end:.0f}] in {min_ifos} detectors{why}: not released yet "
        "or detectors off. " + ("LVK data are configured but no SciToken was found (htgettoken -a "
                                "vault.ligo.org -i igwn)." if lvk is not None else
                                "The LVK back-end needs credentials and its frame types (--source lvk)."))


def segment_summary(source: DataSource, ifos: Sequence[str], start: float, end: float) -> dict:
    """Per detector: analysable and CAT2 time in the window (for the result file)."""
    out = {}
    for ifo in ifos:
        try:
            out[ifo] = {"analysable": source.analysable(ifo, start, end),
                        "cat2": source.vetoed(ifo, start, end)}
        except Exception as e:
            out[ifo] = {"error": f"{type(e).__name__}: {e}"}
    return out
