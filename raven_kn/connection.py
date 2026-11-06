"""Connection module: kilonova candidate events and data connectors.

Analogue of ``ligo.raven.gracedb_events``. :class:`KNCandidate` plays the
role of RAVEN's ``ExtTrig``: one object gathering everything known about an
external (here: kilonova) trigger, whichever service it came from.

Connectors only translate external formats into :class:`KNCandidate` fields;
no physics is done here. Network clients are injectable so that tests and
offline use never touch the network.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Iterable, Optional

import numpy as np
from astropy.time import Time


# --------------------------------------------------------------------------
# Data model
# --------------------------------------------------------------------------

@dataclass
class PhotometryPoint:
    """One photometric measurement or non-detection.

    ``upper_limit=True`` means a non-detection: ``mag`` is then the limiting
    magnitude (5 sigma unless stated). Non-detections *before* the first
    detection are what bracket the merger time (Jamsin et al.,
    arXiv:2610.28775).
    """
    mjd: float
    band: str
    mag: float
    mag_err: Optional[float] = None
    upper_limit: bool = False
    survey: str = ""

    @property
    def gps(self) -> float:
        return float(Time(self.mjd, format="mjd", scale="utc").gps)


@dataclass
class Redshift:
    z: float
    z_err: float = 1e-4
    kind: str = "spectroscopic"      # spectroscopic | photometric
    source: str = ""                 # e.g. "host NGC 4993", "DESI", "kilonova spectrum"


@dataclass
class DistanceMeasurement:
    """Redshift-independent distance (host SBF, TF, Cepheids...) or a
    distance from a kilonova/afterglow light-curve fit."""
    mean_mpc: float
    std_mpc: float
    method: str = ""


@dataclass
class ViewingAngleMeasurement:
    """Constraint on the viewing angle theta_obs = min(iota, pi - iota),
    e.g. from VLBI superluminal motion, afterglow modelling, or a kilonova
    light-curve fit."""
    mean_deg: float
    std_deg: float
    method: str = ""


@dataclass
class GRBAssociation:
    gps: float
    instrument: str = ""
    graceid: Optional[str] = None


@dataclass
class KNCandidate:
    """A kilonova candidate seen as an external trigger (RAVEN ``ExtTrig``).

    Only ``name``, ``ra`` and ``dec`` are mandatory; every other field is a
    piece of information that, when present, tightens the merger-time
    window, the distance or the inclination.
    """
    name: str
    ra: float
    dec: float
    pos_err_arcsec: float = 0.5
    source: str = "manual"              # broker / survey / circular
    kn_score: Optional[float] = None    # classifier probability of being a KN
    graceid: Optional[str] = None       # if uploaded to GraceDB as External event
    photometry: list[PhotometryPoint] = field(default_factory=list)
    redshift: Optional[Redshift] = None
    distance: Optional[DistanceMeasurement] = None
    viewing_angle: Optional[ViewingAngleMeasurement] = None
    grb: Optional[GRBAssociation] = None
    t0_samples_gps: Optional[list[float]] = None   # e.g. NMMA/Redback t0 posterior
    lc_fit_samples: Optional[dict[str, list[float]]] = None  # luminosity_distance, inclination...
    lc_fit_informed: bool = False   # fit used the EM distance/angle as priors: don't reuse its D/iota
    notes: str = ""

    # --- convenience -----------------------------------------------------
    @property
    def detections(self) -> list[PhotometryPoint]:
        return sorted((p for p in self.photometry if not p.upper_limit),
                      key=lambda p: p.mjd)

    @property
    def upper_limits(self) -> list[PhotometryPoint]:
        return sorted((p for p in self.photometry if p.upper_limit),
                      key=lambda p: p.mjd)

    @property
    def first_detection(self) -> Optional[PhotometryPoint]:
        d = self.detections
        return d[0] if d else None

    @property
    def gpstime(self) -> float:
        """Reference time, as in RAVEN ``ExtTrig.gpstime``: the GRB time if
        known, else the first detection."""
        if self.grb is not None:
            return self.grb.gps
        fd = self.first_detection
        if fd is None:
            raise ValueError(f"{self.name}: no detection and no GRB time")
        return fd.gps

    def add_photometry(self, points: Iterable[PhotometryPoint]) -> None:
        """Merge photometry from another connector, dropping exact duplicates."""
        seen = {(round(p.mjd, 6), p.band, p.upper_limit, p.survey)
                for p in self.photometry}
        for p in points:
            key = (round(p.mjd, 6), p.band, p.upper_limit, p.survey)
            if key not in seen:
                self.photometry.append(p)
                seen.add(key)

    # --- (de)serialisation ------------------------------------------------
    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "KNCandidate":
        d = dict(d)
        d["photometry"] = [PhotometryPoint(**p) for p in d.get("photometry", [])]
        for key, klass in (("redshift", Redshift),
                           ("distance", DistanceMeasurement),
                           ("viewing_angle", ViewingAngleMeasurement),
                           ("grb", GRBAssociation)):
            if d.get(key) is not None:
                d[key] = klass(**d[key])
        return cls(**d)

    def to_raven_event_dict(self, pipeline: str = "RAVEN-KN",
                            search: str = "Kilonova") -> dict[str, Any]:
        """Dictionary shaped like a GraceDB External event, so that
        ``ligo.raven.gracedb_events.ExtTrig(graceid, event_dict=...,
        use_radec=True)`` accepts it unchanged. RAVEN reads the position
        from ``extra_attributes['GRB']``, hence that key."""
        return {
            "graceid": self.graceid or self.name,
            "pipeline": pipeline,
            "search": search,
            "gpstime": self.gpstime,
            "extra_attributes": {"GRB": {"ra": self.ra, "dec": self.dec,
                                         "error_radius": self.pos_err_arcsec / 3600.}},
        }


# --------------------------------------------------------------------------
# Connectors
# --------------------------------------------------------------------------

def load_candidate(path: str | Path) -> KNCandidate:
    """Load a candidate from the package's JSON format (see examples/)."""
    with open(path) as f:
        return KNCandidate.from_dict(json.load(f))


def save_candidate(kn: KNCandidate, path: str | Path) -> None:
    with open(path, "w") as f:
        json.dump(kn.to_dict(), f, indent=2)


class FinkConnector:
    """Fetch light curves (with upper limits) from the Fink broker REST API.

    The default endpoint and column names are those of the ZTF Fink API
    (``/api/v1/objects`` with ``withupperlim``). The LSST Fink API uses a
    different base URL and Rubin alert schema: pass ``base_url`` and a
    ``parser`` adapted to it; check the schema against the live service
    before relying on it.
    """

    ZTF_BANDS = {1: "g", 2: "r", 3: "i"}

    def __init__(self, base_url: str = "https://api.fink-portal.org",
                 session: Any = None, parser=None):
        self.base_url = base_url.rstrip("/")
        self._session = session
        self.parser = parser or self.parse_ztf_rows

    @property
    def session(self):
        if self._session is None:
            import requests
            self._session = requests.Session()
        return self._session

    def fetch_rows(self, object_id: str) -> list[dict[str, Any]]:
        r = self.session.post(f"{self.base_url}/api/v1/objects",
                              json={"objectId": object_id,
                                    "withupperlim": "True",
                                    "output-format": "json"},
                              timeout=60)
        r.raise_for_status()
        return r.json()

    @classmethod
    def parse_ztf_rows(cls, rows: list[dict[str, Any]]) -> list[PhotometryPoint]:
        pts = []
        for row in rows:
            tag = row.get("d:tag", "valid")
            if tag == "badquality":
                continue
            mjd = float(row["i:jd"]) - 2400000.5
            band = cls.ZTF_BANDS.get(int(row["i:fid"]), str(row["i:fid"]))
            if tag == "upperlim":
                pts.append(PhotometryPoint(mjd, band, float(row["i:diffmaglim"]),
                                           upper_limit=True, survey="ZTF"))
            else:
                pts.append(PhotometryPoint(mjd, band, float(row["i:magpsf"]),
                                           float(row["i:sigmapsf"]), survey="ZTF"))
        return pts

    def candidate(self, object_id: str, **extra) -> KNCandidate:
        rows = self.fetch_rows(object_id)
        if not rows:
            raise LookupError(f"Fink returned no data for {object_id}")
        valid = [r for r in rows if r.get("d:tag", "valid") == "valid"]
        ref = valid[0] if valid else rows[0]
        return KNCandidate(name=object_id, ra=float(ref["i:ra"]),
                           dec=float(ref["i:dec"]), source="Fink",
                           photometry=self.parser(rows), **extra)

    def add_to(self, kn: KNCandidate, object_id: str) -> KNCandidate:
        """Add this survey's photometry (and pre-discovery limits) to a
        candidate found elsewhere: the multi-survey last-non-detection step."""
        kn.add_photometry(self.parser(self.fetch_rows(object_id)))
        return kn


class ALeRCEConnector:
    """ZTF light curves with non-detections from the ALeRCE broker (public API,
    ``/ztf/v1/objects/<oid>/lightcurve``). Non-detections carry the
    difference-image limiting magnitude ``diffmaglim``: exactly what the
    t0 bracket needs."""

    BANDS = {1: "g", 2: "r", 3: "i"}

    def __init__(self, base_url: str = "https://api.alerce.online/ztf/v1", session: Any = None):
        self.base_url = base_url.rstrip("/")
        self._session = session

    @property
    def session(self):
        if self._session is None:
            import requests
            self._session = requests.Session()
        return self._session

    def fetch(self, oid: str) -> dict[str, Any]:
        r = self.session.get(f"{self.base_url}/objects/{oid}/lightcurve", timeout=60)
        r.raise_for_status()
        return r.json()

    @classmethod
    def parse(cls, lc: dict[str, Any]) -> list[PhotometryPoint]:
        pts = []
        for d in lc.get("detections", []):
            if d.get("dubious") or d.get("isdiffpos") in (-1, "f", "0", 0) or d.get("magpsf") is None:
                continue
            pts.append(PhotometryPoint(float(d["mjd"]), cls.BANDS.get(int(d["fid"]), str(d["fid"])),
                                       float(d["magpsf"]), float(d.get("sigmapsf") or 0.1),
                                       survey="ZTF"))
        for n in lc.get("non_detections", []):
            if n.get("diffmaglim") is None:
                continue
            pts.append(PhotometryPoint(float(n["mjd"]), cls.BANDS.get(int(n["fid"]), str(n["fid"])),
                                       float(n["diffmaglim"]), upper_limit=True, survey="ZTF"))
        return pts

    def candidate(self, oid: str, **extra) -> KNCandidate:
        lc = self.fetch(oid)
        dets = [d for d in lc.get("detections", []) if d.get("ra") is not None]
        if not dets:
            raise LookupError(f"ALeRCE: no detection for {oid}")
        return KNCandidate(name=oid, ra=float(np.median([d["ra"] for d in dets])),
                           dec=float(np.median([d["dec"] for d in dets])), source="ALeRCE",
                           photometry=self.parse(lc), **extra)

    def add_to(self, kn: KNCandidate, oid: str) -> KNCandidate:
        kn.add_photometry(self.parse(self.fetch(oid)))
        return kn


class ATLASForcedPhotometryConnector:
    """ATLAS forced photometry server (https://fallingstar-data.com/forcedphot/),
    queue-based REST API with a personal token. Gives forced flux at the
    candidate position on every ATLAS image: detections and 5-sigma limits
    (``mag5sig``) in cyan 'c' / orange 'o', mapped to NMMA ``atlasc``/``atlaso``.

    Interface written from the server's API guide; needs a token to run
    (not tested against the live service here).
    """

    def __init__(self, token: str, base_url: str = "https://fallingstar-data.com/forcedphot",
                 session: Any = None, snr_detect: float = 5.0, poll_s: float = 10.0,
                 max_wait_s: float = 1800.0):
        self.token, self.base_url = token, base_url.rstrip("/")
        self._session = session
        self.snr_detect, self.poll_s, self.max_wait_s = snr_detect, poll_s, max_wait_s

    @property
    def session(self):
        if self._session is None:
            import requests
            self._session = requests.Session()
        return self._session

    def _headers(self):
        return {"Authorization": f"Token {self.token}", "Accept": "application/json"}

    def request(self, ra: float, dec: float, mjd_min: float, mjd_max: Optional[float] = None) -> str:
        import time
        data = {"ra": ra, "dec": dec, "mjd_min": mjd_min, "send_email": False}
        if mjd_max is not None:
            data["mjd_max"] = mjd_max
        r = self.session.post(f"{self.base_url}/queue/", headers=self._headers(), data=data, timeout=60)
        r.raise_for_status()
        task_url = r.json()["url"]
        waited = 0.0
        while True:
            r = self.session.get(task_url, headers=self._headers(), timeout=60)
            r.raise_for_status()
            j = r.json()
            if j.get("finishtimestamp"):
                if not j.get("result_url"):
                    raise RuntimeError(f"ATLAS task failed: {j.get('error_msg')}")
                return self.session.get(j["result_url"], headers=self._headers(), timeout=60).text
            if waited > self.max_wait_s:
                raise TimeoutError("ATLAS forced photometry still queued")
            time.sleep(self.poll_s)
            waited += self.poll_s

    def parse(self, text: str) -> list[PhotometryPoint]:
        lines = [ln for ln in text.splitlines() if ln.strip()]
        head = lines[0].lstrip("#").split()
        col = {k: i for i, k in enumerate(head)}
        pts = []
        for ln in lines[1:]:
            v = ln.split()
            try:
                mjd, flux, dflux = float(v[col["MJD"]]), float(v[col["uJy"]]), float(v[col["duJy"]])
                band, lim = v[col["F"]], float(v[col["mag5sig"]])
            except (KeyError, ValueError, IndexError):
                continue
            if dflux > 0 and flux / dflux >= self.snr_detect and flux > 0:
                mag = 23.9 - 2.5 * np.log10(flux)
                pts.append(PhotometryPoint(mjd, band, float(mag), float(1.0857 * dflux / flux),
                                           survey="ATLAS"))
            elif lim > 0:
                pts.append(PhotometryPoint(mjd, band, lim, upper_limit=True, survey="ATLAS"))
        return pts

    def add_to(self, kn: KNCandidate, days_before: float = 10.0) -> KNCandidate:
        fd = kn.first_detection
        t_ref = fd.mjd if fd else Time.now().mjd
        kn.add_photometry(self.parse(self.request(kn.ra, kn.dec, t_ref - days_before)))
        return kn


def last_nondetection(kn: KNCandidate, depth_mag: Optional[float] = None) -> Optional[PhotometryPoint]:
    """Latest non-detection before the first detection, over all surveys merged
    into the candidate; with ``depth_mag``, only limits at least that deep."""
    fd = kn.first_detection
    if fd is None:
        return None
    lims = [p for p in kn.upper_limits if p.mjd < fd.mjd and (depth_mag is None or p.mag >= depth_mag)]
    return max(lims, key=lambda p: p.mjd) if lims else None


class GraceDbConnector:
    """Thin wrapper over ``ligo.gracedb.rest.GraceDb`` with the calls RAVEN
    uses: superevent queries, log messages, file downloads.

    Any object with the same methods (e.g. a mock) can be passed as
    ``client``.
    """

    def __init__(self, client: Any = None,
                 service_url: str = "https://gracedb.ligo.org/api/"):
        if client is None:
            from ligo.gracedb.rest import GraceDb
            client = GraceDb(service_url=service_url)
        self.client = client

    def superevents(self, gps_start: float, gps_end: float,
                    far_max: Optional[float] = None) -> list[dict[str, Any]]:
        res = list(self.client.superevents(query=f"{gps_start} .. {gps_end}"))
        if far_max is not None:
            res = [s for s in res if s.get("far") is not None and s["far"] <= far_max]
        return res

    def write_log(self, graceid: str, message: str, filename: Optional[str] = None,
                  filecontents: Optional[str] = None, tags: Iterable[str] = ()):
        return self.client.writeLog(graceid, message, filename=filename,
                                    filecontents=filecontents, tag_name=list(tags))

    def download(self, graceid: str, filename: str, out: str | Path) -> Path:
        out = Path(out)
        out.write_bytes(self.client.files(graceid, filename, raw=True).read())
        return out


def candidate_from_gcn_json(msg: dict[str, Any]) -> GRBAssociation:
    """Extract a GRB time from a GCN Kafka JSON notice (unified schema:
    ``trigger_time`` ISO string, ``instrument``)."""
    t = Time(msg["trigger_time"].replace("Z", ""), format="isot", scale="utc")
    return GRBAssociation(gps=float(t.gps), instrument=msg.get("instrument", ""),
                          graceid=msg.get("id"))


def mjd_to_gps(mjd: float | np.ndarray) -> float | np.ndarray:
    return Time(mjd, format="mjd", scale="utc").gps


def gps_to_mjd(gps: float | np.ndarray) -> float | np.ndarray:
    return Time(gps, format="gps").utc.mjd
