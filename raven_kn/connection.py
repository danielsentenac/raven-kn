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
    classification: Optional[str] = None   # spectroscopic class (e.g. TNS "Kilonova")
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


def helio_to_cmb(z_helio: float, ra: float, dec: float) -> float:
    """Heliocentric -> CMB-frame redshift (solar motion 369.82 km/s towards
    l = 264.021, b = 48.253 deg, Planck 2018): 1 + z_cmb = (1 + z_h)(1 + v.n/c).
    GLADE+ ``zcmb`` is already in this frame; NED and TNS give heliocentric."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    n = SkyCoord(ra * u.deg, dec * u.deg).galactic
    apex = SkyCoord(l=264.021 * u.deg, b=48.253 * u.deg, frame="galactic")
    cos = float(np.cos(n.separation(apex).rad))
    return float((1 + z_helio) * (1 + 369.82 * cos / 299792.458) - 1)


@dataclass
class HostGalaxy:
    """A galaxy near the candidate (GLADE+ row) with its association measures."""
    name: str
    ra: float
    dec: float
    z: Optional[float]               # CMB-frame redshift (GLADE+ zcmb; NED converted)
    z_err: Optional[float]           # measurement error (peculiar velocity added downstream)
    z_kind: str                      # spectroscopic | photometric | none
    d_l_mpc: Optional[float]
    b_mag: Optional[float]
    sep_arcsec: float
    offset_kpc: Optional[float]
    p_chance: Optional[float]        # probability of a chance alignment (Bloom et al. 2002)
    catalog: str = "GLADE+"


class GladePlusConnector:
    """Host-galaxy redshifts from GLADE+ (Dalya et al. 2022, VizieR VII/291),
    the galaxy catalogue LVK uses for GW host searches: 22 million galaxies,
    spectroscopic redshifts where available (``f_dL`` 2-3), photometric
    otherwise (1), CMB-frame redshifts corrected for peculiar velocities.

    Queried over HTTP through a VizieR mirror (no astroquery dependency).
    Hosts are ranked by the probability of chance coincidence
    P_cc = 1 - exp(-pi r^2 sigma(<= m)) (Bloom et al. 2002, number counts of
    Hogg et al. 1997 with the GLADE+ B magnitude, r = offset), the usual
    criterion for short-GRB and kilonova hosts, after a cut on the projected
    offset (``max_offset_kpc``, default 50 kpc: short-GRB hosts, Fong &
    Berger 2013, reach a few tens of kpc).
    """

    COLUMNS = ("GLADE+", "GWGC", "HyperLEDA", "RAJ2000", "DEJ2000", "Bmag", "zcmb",
               "e_zhelio", "dL", "f_dL")

    MIRRORS = ("https://vizier.cfa.harvard.edu/viz-bin/asu-tsv",
               "https://vizier.idia.ac.za/viz-bin/asu-tsv",
               "https://vizier.cds.unistra.fr/viz-bin/asu-tsv")

    def __init__(self, base_url: Optional[str] = None, session: Any = None, timeout: float = 60.0):
        self.mirrors = (base_url,) if base_url else self.MIRRORS
        self.timeout = timeout
        self._session = session

    @property
    def session(self):
        if self._session is None:
            import requests
            self._session = requests.Session()
        return self._session

    def fetch(self, ra: float, dec: float, radius_arcsec: float, max_rows: int = 500) -> str:
        params = {"-source": "VII/291/gladep", "-c": f"{ra:.6f} {dec:+.6f}",
                  "-c.rs": f"{radius_arcsec:.1f}", "-out": ",".join(("_r",) + self.COLUMNS),
                  "-out.max": str(max_rows), "-sort": "_r"}
        errors = []
        for url in self.mirrors:     # VizieR mirrors are sometimes up with their database down
            try:
                r = self.session.get(url, params=params, timeout=self.timeout)
                r.raise_for_status()
            except Exception as e:
                errors.append(f"{url}: {type(e).__name__}")
                continue
            err = [ln for ln in r.text.splitlines()
                   if ln.startswith("#INFO") and "Error=" in ln and ln.split("Error=", 1)[1].strip()]
            if err:
                errors.append(f"{url}: {err[0].split('Error=', 1)[1].strip()}")
                continue
            return r.text
        raise ConnectionError("GLADE+ query failed on all VizieR mirrors: " + "; ".join(errors))

    @staticmethod
    def p_chance(b_mag: Optional[float], r_arcsec: float) -> Optional[float]:
        if b_mag is None or not np.isfinite(b_mag):
            return None
        sigma = 10 ** (0.33 * (b_mag - 24) - 2.44) / (0.33 * np.log(10))   # arcsec^-2
        return float(1 - np.exp(-np.pi * max(r_arcsec, 1.0) ** 2 * sigma))

    @classmethod
    def parse(cls, text: str) -> list[HostGalaxy]:
        from astropy.cosmology import Planck18
        lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
        if len(lines) < 3:
            return []
        head = lines[0].split("\t")
        out = []
        for ln in lines[3:]:                   # header, units, dashes
            row = dict(zip(head, (x.strip() for x in ln.split("\t"))))

            def num(k):
                try:
                    return float(row.get(k, ""))
                except ValueError:
                    return None
            flag = int(num("f_dL") or 0)
            z = num("zcmb")
            kind = {1: "photometric", 2: "spectroscopic", 3: "spectroscopic"}.get(flag, "none")
            name = next((row[k] for k in ("GWGC", "HyperLEDA") if row.get(k, "-") not in ("", "-")),
                        f"GLADE+ {row.get('GLADE+', '?')}")
            sep = num("_r") or 0.0
            off = None
            if z and z > 0:
                d_a = Planck18.angular_diameter_distance(z).to("kpc").value
                off = float(np.radians(sep / 3600) * d_a)
            b = num("Bmag")
            out.append(HostGalaxy(name, num("RAJ2000"), num("DEJ2000"), z, num("e_zhelio"), kind,
                                  num("dL"), b, sep, off, cls.p_chance(b, sep)))
        return out

    def hosts(self, ra: float, dec: float, radius_arcsec: float = 600.0,
              max_offset_kpc: float = 50.0) -> list[HostGalaxy]:
        """Galaxies with a redshift within ``max_offset_kpc`` (projected),
        most probable host (lowest P_cc) first."""
        gals = [g for g in self.parse(self.fetch(ra, dec, radius_arcsec))
                if g.z is not None and g.offset_kpc is not None and g.offset_kpc <= max_offset_kpc]
        return sorted(gals, key=lambda g: (g.p_chance if g.p_chance is not None else 1.0,
                                           g.offset_kpc))

    def add_to(self, kn: KNCandidate, max_p_chance: float = 0.1, overwrite: bool = False,
               **kw) -> Optional[HostGalaxy]:
        """Set ``kn.redshift`` from the most probable host if its P_cc is below
        ``max_p_chance``; a redshift already present (e.g. from a spectrum)
        is kept unless ``overwrite``. Returns the host used, if any."""
        hosts = self.hosts(kn.ra, kn.dec, **kw)
        if not hosts or (hosts[0].p_chance is not None and hosts[0].p_chance > max_p_chance):
            return None
        h = hosts[0]
        if kn.redshift is None or overwrite:
            kn.redshift = Redshift(h.z, h.z_err or (1e-4 if h.z_kind == "spectroscopic" else 0.03),
                                   kind=h.z_kind,
                                   source=f"host {h.name} ({h.catalog}, offset {h.offset_kpc:.1f} kpc, "
                                          f"P_cc {h.p_chance:.2g})")
        return h


class NEDConnector:
    """NASA/IPAC Extragalactic Database: redshifts of galaxies near the
    candidate (TAP, ``NEDTAP.objdir``: SDSS, DESI, 6dF, literature, often
    more complete in spectroscopic redshifts than GLADE+ beyond 100 Mpc) and
    redshift-independent distances of a galaxy (NED-D: SBF, TRGB, Cepheids,
    SNe Ia, Tully-Fisher, FP...).

    NED-D distances obtained *from gravitational waves* are always dropped:
    the EM distance must stay independent of the GW data it is compared to.
    """

    TAP = "https://ned.ipac.caltech.edu/tap/sync"
    NEDD = "https://ned.ipac.caltech.edu/cgi-bin/nDistance"
    GALAXY_TYPES = ("G", "GPair", "GTrpl", "GGroup", "PofG")
    # preferred redshift-independent methods, most reliable first
    METHODS = ("Cepheids", "TRGB", "SBF", "SNIa", "Maser", "FP", "Tully-Fisher", "GCLF")

    def __init__(self, session: Any = None):
        self._session = session

    @property
    def session(self):
        if self._session is None:
            import requests
            self._session = requests.Session()
        return self._session

    def fetch_cone(self, ra: float, dec: float, radius_arcsec: float) -> str:
        q = (f"SELECT TOP 200 prefname, prefphytype, ra, dec, z, zunc, zflag, n_dist, "
             f"DISTANCE(POINT('ICRS',ra,dec), POINT('ICRS',{ra:.7f},{dec:.7f}))*3600 AS sep "
             f"FROM NEDTAP.objdir WHERE CONTAINS(POINT('ICRS',ra,dec), "
             f"CIRCLE('ICRS',{ra:.7f},{dec:.7f},{radius_arcsec / 3600:.6f}))=1 "
             f"AND z IS NOT NULL ORDER BY sep")
        r = self.session.post(self.TAP, data={"QUERY": q, "LANG": "ADQL", "REQUEST": "doQuery",
                                              "FORMAT": "csv"}, timeout=120)
        r.raise_for_status()
        return r.text

    @classmethod
    def parse_cone(cls, text: str) -> list[tuple[HostGalaxy, int]]:
        """(galaxy, number of NED-D distances) for galaxy-type objects."""
        import csv
        import io
        from astropy.cosmology import Planck18
        out = []
        for row in csv.DictReader(io.StringIO(text)):
            if row.get("prefphytype", "") not in cls.GALAXY_TYPES or row["prefphytype"] == "PofG":
                continue
            try:
                z = float(row["z"])
            except (TypeError, ValueError):
                continue
            zunc = float(row["zunc"]) if row.get("zunc") else None
            # NED redshift flag: first letter P = photometric; spectroscopic otherwise
            kind = "photometric" if row.get("zflag", "").startswith("P") else "spectroscopic"
            sep = float(row["sep"])
            off = (float(np.radians(sep / 3600) * Planck18.angular_diameter_distance(z).to("kpc").value)
                   if z > 0 else None)
            z = helio_to_cmb(z, float(row["ra"]), float(row["dec"]))
            out.append((HostGalaxy(row["prefname"], float(row["ra"]), float(row["dec"]), z, zunc, kind,
                                   None, None, sep, off, None, catalog="NED"),
                        int(row.get("n_dist") or 0)))
        return out

    def hosts(self, ra: float, dec: float, radius_arcsec: float = 600.0,
              max_offset_kpc: float = 50.0) -> list[tuple[HostGalaxy, int]]:
        """Galaxies with a redshift within ``max_offset_kpc``, spectroscopic
        first, then by projected offset."""
        gals = [(g, n) for g, n in self.parse_cone(self.fetch_cone(ra, dec, radius_arcsec))
                if g.offset_kpc is not None and g.offset_kpc <= max_offset_kpc]
        return sorted(gals, key=lambda gn: (gn[0].z_kind != "spectroscopic", gn[0].offset_kpc))

    def fetch_distances(self, name: str) -> str:
        r = self.session.get(self.NEDD, params={"name": name}, timeout=120)
        r.raise_for_status()
        return r.text

    @classmethod
    def parse_distances(cls, html: str) -> list[dict]:
        """NED-D rows: distance modulus, its error, D (Mpc), method, reference."""
        import re
        text = re.sub(r"<[^>]*>", " ", html)
        text = re.sub(r"\s+", " ", text)
        i = text.find("REFCODE")
        rows = []
        pat = re.compile(r"(\d{2}\.\d{2}) (?:(\d\.\d{2}) )?(\d+\.\d{3}) ([A-Za-z][A-Za-z .()\-]*?) "
                         r"(\d{4}[A-Za-z&.]{4,6}[^ ]*)")
        for m in pat.finditer(text[i:] if i >= 0 else text):
            mu, err, d, method, ref = m.groups()
            rows.append({"mu": float(mu), "mu_err": float(err) if err else None,
                         "d_mpc": float(d), "method": method.strip(), "ref": ref})
        return rows

    def distance(self, name: str) -> Optional[DistanceMeasurement]:
        """Best redshift-independent distance: most reliable method present
        (``METHODS`` order), its most precise measurement; GW distances excluded."""
        rows = [r for r in self.parse_distances(self.fetch_distances(name))
                if r["mu_err"] and "grav" not in r["method"].lower()]
        for meth in self.METHODS:
            sel = [r for r in rows if r["method"].lower().startswith(meth.lower())]
            if sel:
                r = min(sel, key=lambda r: r["mu_err"])
                return DistanceMeasurement(r["d_mpc"], r["d_mpc"] * np.log(10) / 5 * r["mu_err"],
                                           method=f"{r['method']} ({r['ref']}, NED-D)")
        return None


class TNSConnector:
    """Transient Name Server: spectroscopic classification (e.g. "Kilonova"),
    redshift, host and the list of reported spectra (ASCII/FITS files).

    TNS requires a bot: ``bot_id``, ``bot_name`` and ``api_key`` (created by
    the user on wis-tns.org). Requests follow the TNS API
    (``POST <base>/api/get/object``, form fields ``api_key`` and ``data``,
    ``User-Agent: tns_marker{...}``); replies are parsed whether the object is
    under ``data`` or ``data.reply`` (both API versions).
    """

    def __init__(self, bot_id: int | str, bot_name: str, api_key: str,
                 base_url: str = "https://www.wis-tns.org", session: Any = None):
        self.bot_id, self.bot_name, self.api_key = str(bot_id), bot_name, api_key
        self.base_url = base_url.rstrip("/")
        self._session = session

    @property
    def session(self):
        if self._session is None:
            import requests
            self._session = requests.Session()
        return self._session

    def _headers(self):
        return {"User-Agent": 'tns_marker{"tns_id": ' + self.bot_id
                + ', "type": "bot", "name": "' + self.bot_name + '"}'}

    def fetch(self, objname: str) -> dict[str, Any]:
        name = objname.replace("AT", "").replace("SN", "").strip()
        data = {"objname": name, "objid": "", "photometry": "0", "spectra": "1"}
        r = self.session.post(f"{self.base_url}/api/get/object", headers=self._headers(),
                              data={"api_key": self.api_key, "data": json.dumps(data)}, timeout=60)
        r.raise_for_status()
        return r.json()

    @staticmethod
    def parse(reply: dict[str, Any]) -> dict[str, Any]:
        d = reply.get("data", reply)
        d = d.get("reply", d)
        otype = d.get("object_type") or {}
        tname = otype.get("name") if isinstance(otype, dict) else otype
        z = d.get("redshift")
        hz = d.get("host_redshift")
        spectra = []
        for sp in d.get("spectra") or []:
            inst = sp.get("instrument")
            spectra.append({"jd": sp.get("jd"), "obsdate": sp.get("obsdate"),
                            "instrument": inst.get("name") if isinstance(inst, dict) else inst,
                            "ascii": sp.get("asciifile"), "fits": sp.get("fitsfile"),
                            "group": (sp.get("source_group") or {}).get("name")
                            if isinstance(sp.get("source_group"), dict) else sp.get("source_group")})
        return {"name": f"{d.get('name_prefix', '')}{d.get('objname', '')}".strip(),
                "ra": d.get("radeg"), "dec": d.get("decdeg"), "type": tname,
                "redshift": float(z) if z not in (None, "") else None,
                "host": d.get("hostname"),
                "host_redshift": float(hz) if hz not in (None, "") else None,
                "spectra": spectra}

    def add_to(self, kn: KNCandidate, objname: str, overwrite: bool = False) -> dict[str, Any]:
        """Classification and redshift into the candidate. A redshift reported
        with a spectroscopic classification counts as spectroscopic (host or
        transient lines); the host redshift is preferred when both exist."""
        info = self.parse(self.fetch(objname))
        if info["type"]:
            kn.classification = info["type"]
        z = info["host_redshift"] if info["host_redshift"] is not None else info["redshift"]
        if z is not None:
            z = helio_to_cmb(z, kn.ra, kn.dec)
        if z is not None and (kn.redshift is None or overwrite or kn.redshift.kind != "spectroscopic"):
            src = f"TNS {info['name']}: " + (f"host {info['host']}" if info["host_redshift"] is not None
                                             else "transient redshift")
            kind = "spectroscopic" if (info["type"] or info["spectra"]) else "unknown"
            kn.redshift = Redshift(z, 1e-3 if info["host_redshift"] is None else 1e-4, kind=kind,
                                   source=src)
        return info


def find_host(kn: KNCandidate, ned: Optional[NEDConnector] = None,
              glade: Optional[GladePlusConnector] = None, max_offset_kpc: float = 50.0,
              max_p_chance: float = 0.1, nedd_max_z: float = 0.03,
              overwrite: bool = False) -> dict[str, Any]:
    """Host galaxy, redshift and (if nearby) redshift-independent distance.

    Order: a spectroscopic redshift already in the candidate (from a
    spectrum, TNS or a circular) is kept unless ``overwrite``; otherwise
    NED's nearest galaxy with a spectroscopic redshift within
    ``max_offset_kpc``; otherwise GLADE+ (lowest chance probability, may be
    photometric). For a host closer than ``nedd_max_z`` with NED-D entries,
    the best redshift-independent distance is added (peculiar velocities
    dominate the redshift error there). Each step is logged in the result.
    """
    ned = ned or NEDConnector()
    log: dict[str, Any] = {"steps": []}
    have_spec = kn.redshift is not None and kn.redshift.kind == "spectroscopic" and not overwrite
    host, n_dist = None, 0
    try:
        cands = ned.hosts(kn.ra, kn.dec, max_offset_kpc=max_offset_kpc)
        spec = [(g, n) for g, n in cands if g.z_kind == "spectroscopic"]
        if spec:
            host, n_dist = spec[0]
        log["steps"].append(f"NED: {len(cands)} galaxies with z within {max_offset_kpc} kpc, "
                            f"{len(spec)} spectroscopic")
    except Exception as e:                       # service down: fall back
        log["steps"].append(f"NED failed: {type(e).__name__}")
    if host is None:
        try:
            glade = glade or GladePlusConnector()
            gl = glade.hosts(kn.ra, kn.dec, max_offset_kpc=max_offset_kpc)
            if gl and (gl[0].p_chance is None or gl[0].p_chance <= max_p_chance):
                host = gl[0]
            log["steps"].append(f"GLADE+: {len(gl)} galaxies within {max_offset_kpc} kpc")
        except Exception as e:
            log["steps"].append(f"GLADE+ failed: {type(e).__name__}")
    if host is not None:
        log["host"] = asdict(host)
        if not have_spec and (kn.redshift is None or overwrite or host.z_kind == "spectroscopic"):
            kn.redshift = Redshift(host.z, host.z_err or (1e-4 if host.z_kind == "spectroscopic" else 0.03),
                                   kind=host.z_kind,
                                   source=f"host {host.name} ({host.catalog}, offset "
                                          f"{host.offset_kpc:.1f} kpc)")
            log["steps"].append(f"redshift from {host.catalog} host {host.name}: z={host.z:.5f} "
                                f"({host.z_kind})")
        if host.catalog == "NED" and n_dist > 0 and host.z < nedd_max_z and kn.distance is None:
            try:
                d = ned.distance(host.name)
                if d is not None:
                    kn.distance = d
                    log["steps"].append(f"NED-D distance {d.mean_mpc:.1f} +- {d.std_mpc:.1f} Mpc "
                                        f"({d.method})")
            except Exception as e:
                log["steps"].append(f"NED-D failed: {type(e).__name__}")
    else:
        log["steps"].append("no host with a redshift found")
    return log


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
