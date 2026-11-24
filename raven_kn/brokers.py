"""What the alert brokers know at the position of an optical transient.

:func:`crossmatch` asks each broker for the object nearest to a position and
writes one record per broker: the object name, its separation, the number of
detections and their dates, the classification of the broker, the host
galaxy it proposes and the address of its page. Brokers (interfaces checked
on 2026-11-24):

* **Lasair** (ZTF), ``https://lasair-ztf.lsst.ac.uk/api``: needs a personal
  token (``LASAIR_TOKEN``; 100 calls per hour). Cone search, then the object
  with its Sherlock record (context classification from catalogues: SN, NT,
  AGN, CV, VS, BS, ORPHAN, UNCLEAR; host name, redshift, separation) and its
  TNS name. Without a ZTF object at the position, Sherlock is asked for the
  position alone: the host it proposes is still given.
* **Fink** (ZTF ``https://api.ztf.fink-portal.org``, Rubin
  ``https://api.lsst.fink-portal.org``): cone search, classification.
* **ALeRCE** (ZTF): cone search, classifier and probability.
* **ANTARES**: cone search (locus).

Nothing of this is applied to the transient or enters a decision: a host
redshift of a broker is a proposal from a catalogue match, shown as such.
The records are kept in ``<store>/brokers/<name>.json`` and given by the
service with the transient.

``raven-kn brokers [NAME ...] --store DIR`` runs it by hand; the listener
runs it when a transient is assessed and its record is older than
:data:`REFRESH_H` hours.
"""

from __future__ import annotations

import json
import logging
import os
import re
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Callable, Optional

from .connection import ALeRCEConnector, ANTARESConnector, FinkLSSTConnector, PhotometryPoint
from .journal import _safe, _write_atomic, iso, parse_iso, utcnow

log = logging.getLogger("raven_kn.brokers")

RADIUS_ARCSEC = 3.0
REFRESH_H = 12.0
TIMEOUT_S = 20.0
SHERLOCK = {"SN": "supernova (near a galaxy, off its centre)", "NT": "nuclear transient (at the centre of a galaxy)",
            "AGN": "active galactic nucleus", "CV": "cataclysmic variable", "VS": "variable star",
            "BS": "bright star", "ORPHAN": "no catalogued source nearby", "UNCLEAR": "unclear"}


def _mjd_iso(mjd: Optional[float]) -> Optional[str]:
    if mjd is None:
        return None
    return iso(datetime(1858, 11, 17, tzinfo=timezone.utc) + timedelta(days=float(mjd)))


def _plain(html: Optional[str]) -> Optional[str]:
    """Text of a Sherlock description (it comes with HTML links)."""
    return " ".join(re.sub(r"<[^>]+>", "", html).split()) if html else None


def _pos(x: Any) -> Optional[float]:
    """A positive number, else None (Sherlock writes 0 for unknown)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if v > 0 else None


class LasairConnector:
    """Lasair (ZTF) REST interface. ``token``: the API token of a Lasair
    account (default: the environment variable ``LASAIR_TOKEN``)."""

    BANDS = {1: "g", 2: "r", 3: "i"}
    PAGE = "https://lasair-ztf.lsst.ac.uk/objects/{id}/"

    def __init__(self, token: Optional[str] = None, base_url: str = "https://lasair-ztf.lsst.ac.uk/api",
                 session: Any = None, timeout: float = TIMEOUT_S):
        self.token = token if token is not None else os.environ.get("LASAIR_TOKEN", "")
        self.base_url = base_url.rstrip("/")
        self._session = session
        self.timeout = timeout

    @property
    def session(self):
        if self._session is None:
            import requests
            self._session = requests.Session()
        return self._session

    def _get(self, endpoint: str, **params) -> Any:
        r = self.session.get(f"{self.base_url}/{endpoint}/", params=params, timeout=self.timeout,
                             headers={"Authorization": f"Token {self.token}"})
        r.raise_for_status()
        return r.json()

    def cone(self, ra: float, dec: float, radius_arcsec: float = RADIUS_ARCSEC) -> list[dict[str, Any]]:
        """[{"object", "separation" (arcsec)}], nearest first."""
        rows = self._get("cone", ra=ra, dec=dec, radius=radius_arcsec, requestType="all")
        return sorted(rows, key=lambda r: r["separation"])

    def object(self, object_id: str) -> dict[str, Any]:
        return self._get("object", objectId=object_id, lasair_added="true")

    def sherlock_position(self, ra: float, dec: float) -> dict[str, Any]:
        return self._get("sherlock/position", ra=ra, dec=dec, lite="false")

    @classmethod
    def parse(cls, obj: dict[str, Any]) -> list[PhotometryPoint]:
        """Detections (``candid``, positive difference) and non-detections
        (``diffmaglim``) of an object."""
        pts = []
        for c in obj.get("candidates") or []:
            band = cls.BANDS.get(int(c["fid"]), str(c["fid"]))
            if c.get("candid") is None:
                if c.get("diffmaglim") is not None:
                    pts.append(PhotometryPoint(float(c["mjd"]), band, float(c["diffmaglim"]),
                                               upper_limit=True, survey="ZTF"))
            elif str(c.get("isdiffpos", "t")).lower() in ("t", "1") and c.get("magpsf") is not None:
                pts.append(PhotometryPoint(float(c["mjd"]), band, float(c["magpsf"]),
                                           float(c.get("sigmapsf") or 0.1), survey="ZTF"))
        return sorted(pts, key=lambda p: p.mjd)

    @staticmethod
    def host_of(s: dict[str, Any]) -> Optional[dict[str, Any]]:
        """Host proposed by a Sherlock record or crossmatch row."""
        if not s or not s.get("catalogue_object_id"):
            return None
        return {"name": str(s["catalogue_object_id"]), "catalogue": s.get("catalogue_table_name"),
                "type": s.get("catalogue_object_type"), "z": _pos(s.get("z")), "photo_z": _pos(s.get("photoZ")),
                "photo_z_err": _pos(s.get("photoZErr")), "sep_arcsec": _pos(s.get("separationArcsec")),
                "sep_kpc": _pos(s.get("physical_separation_kpc")),
                "distance_mpc": _pos(s.get("direct_distance")) or _pos(s.get("distance"))}

    def match(self, ra: float, dec: float, radius_arcsec: float = RADIUS_ARCSEC) -> dict[str, Any]:
        if not self.token:
            return {"status": "not_configured", "note": "no Lasair token (LASAIR_TOKEN) on the server"}
        near = self.cone(ra, dec, radius_arcsec)
        if not near:
            # no ZTF object: what the catalogues have at this position
            s = self.sherlock_position(ra, dec)
            cls = next(iter((s.get("classifications") or {}).values()), None)
            best = min(s.get("crossmatches") or [], key=lambda x: x.get("rank") or 99, default=None)
            return {"status": "none", "classification": cls[0] if cls else None,
                    "classification_text": SHERLOCK.get(cls[0]) if cls else None,
                    "description": _plain(cls[1]) if cls and len(cls) > 1 else None,
                    "host": self.host_of(best),
                    "note": "no ZTF object at this position; the classification is Sherlock's, "
                            "from the catalogued sources nearby"}
        oid = near[0]["object"]
        o = self.object(oid)
        data, s, tns = o.get("objectData") or {}, o.get("sherlock") or {}, o.get("TNS") or {}
        c = s.get("classification")
        return {"status": "match", "object": oid, "sep_arcsec": float(near[0]["separation"]),
                "others": [r["object"] for r in near[1:5]], "ra": data.get("ramean"), "dec": data.get("decmean"),
                "n_detections": data.get("ncand"), "first_detection": _mjd_iso(data.get("discMjd")),
                "last_detection": _mjd_iso(data.get("latestMjd")),
                "peak": (f"{data['peakMag']} ({data.get('peakFilter')})" if data.get("peakMag") else None),
                "classification": c, "classification_text": SHERLOCK.get(c),
                "description": _plain(s.get("description")), "host": self.host_of(s),
                "tns": ({"name": f"{tns.get('tns_prefix') or ''} {tns['tns_name']}".strip(),
                         "type": tns.get("type") or None, "z": _pos(tns.get("z"))}
                        if tns.get("tns_name") else None),
                "url": self.PAGE.format(id=oid)}


def _fink_ztf(ra: float, dec: float, radius: float, session: Any = None,
              base_url: str = "https://api.ztf.fink-portal.org") -> dict[str, Any]:
    if session is None:
        import requests
        session = requests
    r = session.post(f"{base_url}/api/v1/conesearch", timeout=TIMEOUT_S,
                     json={"ra": ra, "dec": dec, "radius": radius, "output-format": "json"})
    r.raise_for_status()
    rows = sorted(r.json(), key=lambda x: x.get("v:separation_degree") or 0)
    if not rows:
        return {"status": "none"}
    x = rows[0]
    return {"status": "match", "object": x["i:objectId"],
            "sep_arcsec": float(x.get("v:separation_degree") or 0) * 3600, "ra": x.get("i:ra"), "dec": x.get("i:dec"),
            "n_detections": x.get("d:nalerthist"), "last_detection": _mjd_iso(float(x["i:jd"]) - 2400000.5),
            "classification": x.get("d:classification"), "others": [y["i:objectId"] for y in rows[1:5]],
            "url": f"https://ztf.fink-portal.org/{x['i:objectId']}"}


def _fink_lsst(ra: float, dec: float, radius: float, session: Any = None) -> dict[str, Any]:
    ids = FinkLSSTConnector(session=session).conesearch(ra, dec, radius)
    if not ids:
        return {"status": "none"}
    return {"status": "match", "object": ids[0], "others": ids[1:5],
            "url": f"https://lsst.fink-portal.org/{ids[0]}"}


def _alerce(ra: float, dec: float, radius: float, session: Any = None) -> dict[str, Any]:
    c = ALeRCEConnector(session=session)
    r = c.session.get(f"{c.base_url}/objects/", params={"ra": ra, "dec": dec, "radius": radius}, timeout=TIMEOUT_S)
    r.raise_for_status()
    items = r.json().get("items") or []
    if not items:
        return {"status": "none"}
    from .connection import _sep_arcsec
    items.sort(key=lambda x: _sep_arcsec(ra, dec, x["meanra"], x["meandec"]))
    x = items[0]
    return {"status": "match", "object": x["oid"], "sep_arcsec": _sep_arcsec(ra, dec, x["meanra"], x["meandec"]),
            "ra": x["meanra"], "dec": x["meandec"], "n_detections": x.get("ndet"),
            "first_detection": _mjd_iso(x.get("firstmjd")), "last_detection": _mjd_iso(x.get("lastmjd")),
            "classification": x.get("class"),
            "classification_text": (f"{x['classifier']}, probability {float(x['probability']):.2f}"
                                    if x.get("classifier") and x.get("probability") is not None else None),
            "others": [y["oid"] for y in items[1:5]], "url": f"https://alerce.online/object/{x['oid']}"}


def _antares(ra: float, dec: float, radius: float, session: Any = None) -> dict[str, Any]:
    loci = ANTARESConnector(session=session, timeout=TIMEOUT_S).cone(ra, dec, radius)
    if not loci:
        return {"status": "none"}
    x = loci[0]
    return {"status": "match", "object": x["ztf_object_id"] or x["locus_id"], "sep_arcsec": x["sep_arcsec"],
            "ra": x["ra"], "dec": x["dec"], "n_detections": x.get("num_alerts"),
            "others": [y["ztf_object_id"] or y["locus_id"] for y in loci[1:5]],
            "url": f"https://antares.noirlab.edu/loci/{x['locus_id']}"}


def _lasair(ra: float, dec: float, radius: float, session: Any = None) -> dict[str, Any]:
    return LasairConnector(session=session).match(ra, dec, radius)


BROKERS: dict[str, tuple[str, Callable[..., dict]]] = {
    "lasair": ("Lasair (ZTF)", _lasair),
    "fink_ztf": ("Fink (ZTF)", _fink_ztf),
    "fink_lsst": ("Fink (Rubin)", _fink_lsst),
    "alerce": ("ALeRCE (ZTF)", _alerce),
    "antares": ("ANTARES", _antares),
}


def crossmatch(ra: float, dec: float, radius_arcsec: float = RADIUS_ARCSEC,
               brokers: Optional[dict[str, tuple[str, Callable[..., dict]]]] = None) -> dict[str, Any]:
    """Ask each broker for the object nearest to a position. A broker that
    does not answer gives ``status: "error"`` with the reason; the others are
    not affected. Statuses: ``match``, ``none``, ``error``, ``not_configured``."""
    out = []
    for key, (label, fn) in (brokers or BROKERS).items():
        try:
            rec = fn(ra, dec, radius_arcsec)
        except Exception as e:
            log.warning("%s: no answer at %.5f %+.5f (%s)", label, ra, dec, type(e).__name__)
            # the address asked is left out: long, and it could carry a token
            reason = re.split(r" for url| with url|: /|\(Caused by", f"{type(e).__name__}: {e}")[0][:160]
            rec = {"status": "error", "note": reason}
        out.append({"broker": key, "label": label, **rec})
    return {"time": iso(utcnow()), "ra": ra, "dec": dec, "radius_arcsec": radius_arcsec, "brokers": out}


def _path(root: str | Path, name: str) -> Path:
    return Path(root) / "brokers" / f"{_safe(name)}.json"


def stored(root: str | Path, name: str) -> Optional[dict]:
    try:
        return json.loads(_path(root, name).read_text())
    except (OSError, ValueError):
        return None


def crossmatch_transient(kn: Any, root: str | Path, max_age_h: Optional[float] = None,
                         radius_arcsec: float = RADIUS_ARCSEC, **kw) -> dict[str, Any]:
    """Cross-match an optical transient and store the record. With
    ``max_age_h``, a stored record younger than that is returned as it is."""
    old = stored(root, kn.name)
    if old and max_age_h is not None:
        age = (utcnow() - parse_iso(old["time"])).total_seconds() / 3600
        if age < max_age_h:
            return old
    rec = {"optical": kn.name, **crossmatch(kn.ra, kn.dec, radius_arcsec, **kw)}
    # a broker that did not answer this time: its last answer is kept, marked as old
    for b in rec["brokers"]:
        prev = next((p for p in (old or {}).get("brokers", []) if p["broker"] == b["broker"]), None)
        if b["status"] == "error" and prev and prev.get("status") in ("match", "none"):
            b.update({**prev, "stale_since": prev.get("stale_since") or old["time"], "note": b["note"]})
    p = _path(root, kn.name)
    p.parent.mkdir(parents=True, exist_ok=True)
    _write_atomic(p, json.dumps(rec, indent=1))
    return rec
