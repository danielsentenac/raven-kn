"""GCN Circulars: what observers report about the alerts of the store.

Notices give the time and position of an alert. Whether anybody looked for
its optical counterpart, and with what result, is written by the observers
in GCN Circulars (free text). The listener reads them from the Kafka topic
``gcn.circulars`` (JSON: ``circularId``, ``subject``, ``eventId``, ``body``,
``createdOn`` in ms, checked on the archive on 2026-11-19) and:

* **links** each circular to the alerts of the store it is about
  (:func:`find_links`);
* **labels** it from its subject, and from its text when the subject does not
  say: the bands observed (optical, X-ray, radio, gamma, GW, neutrino), for an
  optical report whether it is a detection or upper limits, and the redshift
  when one is announced (:func:`classify`). These are keyword rules, not a
  reading of the text: the subject is kept so that the reader can judge;
* stores it in ``<store>/circulars/<number>.json`` and writes a ``circular``
  entry in the journal. A circular linked to nothing is kept 30 days in
  ``circulars/unlinked/`` when it names an event: the circular that gives the
  trigger number of a burst may come after the first reports about it.

Links, in this order:

1. GW: a superevent of the store named in the subject;
2. neutrino: a stored IceCube event named in the subject;
3. optical: the name of a stored optical transient in the subject or the text;
4. GRB and X-ray notices: the trigger number of a stored notice in the
   subject or the text (Swift ``trigger=1263718``, Fermi ``trigger
   751960088.32731 / 241030242``, SVOM ``sb26111302``);
5. else a time of the text within :data:`TRIGGER_TIME_TOL_S` of the trigger
   time of a stored notice, on the date of the event name (``GRB 261113A``:
   2026-11-13). Einstein Probe notices carry no trigger number;
6. else the event name, when an earlier circular linked by 4 or 5 carried it:
   most follow-up reports only say ``GRB 261113A``.

Nothing is changed in the optical transients: a redshift announced in a
circular is recorded with the circular, not applied.

``raven-kn circular 46261 46270-46290 --store DIR`` reads circulars from the
archive (https://gcn.nasa.gov/circulars/<number>.json), to fill what was sent
before the listener subscribed, or to link again.
"""

from __future__ import annotations

import json
import logging
import re
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Iterable, Optional

from .journal import _safe, _write_atomic, iso, journal_of, parse_iso, utcnow

log = logging.getLogger("raven_kn.circulars")

CIRCULARS_TOPIC = "gcn.circulars"
ARCHIVE_URL = "https://gcn.nasa.gov/circulars/{id}"
TRIGGER_TIME_TOL_S = 1.5         # the circulars give the trigger time to the second
UNLINKED_DAYS = 30.0
MIN_ID_LEN = 6                   # shorter identifiers are not searched for in the text

# Keywords of each band, searched as whole words in the subject (lower case).
# The optical list holds the facilities met in the circulars of October and
# November 2026; a report from another one is still found by "optical",
# "afterglow", "photometry", "spectroscopic"...
BANDS: dict[str, tuple[str, ...]] = {
    "gw": ("ligo", "virgo", "kagra", "lvk", "gravitational"),
    "nu": ("icecube", "km3net", "antares", "neutrino", "neutrinos"),
    "radio": ("radio", "vla", "jvla", "meerkat", "atca", "alma", "noema", "gmrt", "ugmrt", "askap",
              "e-merlin", "vlbi", "evn", "ami-la", "ami", "sma", "lofar", "chime", "mwa",
              "millimeter", "millimetre", "submillimeter", "mm"),
    "xray": ("x-ray", "xray", "xrt", "swift-xrt", "nicer", "nustar", "chandra", "xmm", "xmm-newton",
             "maxi", "fxt", "ep-fxt", "wxt", "ep-wxt", "mxt", "erosita", "leia", "einstein probe"),
    "gamma": ("gamma-ray", "gbm", "bat", "swift-bat", "bat-guano", "konus", "konus-wind", "gecam",
              "gecam-a", "gecam-b", "gecam-c", "astrosat", "czti", "insight-hxmt", "hxmt", "agile",
              "integral", "lat", "fermi-lat", "eclairs", "grm", "calet", "glowbug", "ipn",
              "triangulation", "grbalpha", "vzlusat", "hebs", "grid", "localization", "localisation",
              "burst"),
    "optical": ("optical", "nir", "near-infrared", "infrared", "afterglow", "photometry",
                "photometric", "spectroscopic", "spectroscopy", "spectrum", "redshift",
                "counterpart", "r-band", "i-band", "g-band", "z-band", "j-band", "kilonova",
                "uvot", "telescope",
                "goto", "ztf", "atlas", "master", "master-net", "not", "vlt", "gtc", "gemini",
                "keck", "lco", "colibri", "colibrí", "grond", "kait", "lulin", "mephisto", "trt",
                "xinglong", "sao", "gmg", "grandma", "knc", "kilonova-catcher", "bootes", "tng",
                "lbt", "hst", "jwst", "liverpool", "svom/vt", "vt", "c-gft", "f-gft", "fm-gft",
                "gwac", "kmtnet", "ddoti", "rem", "nuttela-tao", "abao", "mondy", "kitab", "soar",
                "pan-starrs", "decam", "subaru", "wfst", "7dt", "ohp", "nanshan", "next", "mitsume",
                "azt-20", "dfot", "ondrejov", "d50", "fram-auger", "fram", "vst", "omegacam",
                "itelescope", "panzhihua", "x-shooter", "fors2", "osiris", "git", "growth-india",
                "tshao", "crao", "zeiss-1000", "skynet", "pairitel", "ratir", "saao", "lesedi",
                "meerlicht", "blackgem", "gotoxb", "ls4", "wendelstein", "calar", "caha", "tnot",
                "hmt", "ast3", "tido", "spm", "oaj", "t80", "njust", "sysu", "jinshan", "alt-100c",
                "assy-turgen", "gwac-f60", "zwicky", "euclid", "bootes-5", "winter", "kinder",
                "ntt", "terskol", "mosfire"),
}
_UPPER = re.compile(r"upper[- ]limits?|\blimits?\b|non[- ]detections?|\bno (?:optical |new |credible |"
                    r"significant )?(?:counterpart|afterglow|source|candidate|transient)s?\b"
                    r"|not detect|no detection|undetected", re.I)
_DETECT = re.compile(r"detections?\b|detected|discover|counterpart|afterglow|candidate|redshift|"
                     r"spectroscop|confirmation|brightening|fading|plateau|kilonova", re.I)
_BODY_UPPER = re.compile(r"(?:do|did|was|were|is|are) not (?:\w+ )?detect|\bno (?:new |optical |"
                         r"uncatalog\w* |credible |significant )+(?:source|counterpart|afterglow|"
                         r"transient|candidate|ot)\b|upper limit|non-detection|not detected", re.I)
_BODY_DETECT = re.compile(r"\bwe (?:clearly |marginally |also )?detect|(?:is|was) (?:clearly |well |"
                          r"marginally )?detected|we (?:found|find|identify|identified|discovered)"
                          r"|optical (?:counterpart|afterglow) (?:is|was) ", re.I)
_REDSHIFT = re.compile(r"\bz\s*(?:=|~|≈|of|is)\s*(\d{1,2}\.\d{1,5})\b")
_NOT_REDSHIFT = re.compile(r"photometric redshift|photo-?z|z\s*[<>]", re.I)

_SUPEREVENT = re.compile(r"(?<![A-Za-z0-9])([MT]?S\d{6}[a-z]{1,3})(?![A-Za-z0-9])")
_ICECUBE = re.compile(r"IceCube[- _]?(\d{6}[A-Z])", re.I)
_EVENT_DATE = re.compile(r"(?<!\d)(\d{2})(\d{2})(\d{2})(?:[A-Za-z]|\b)")
_ISO_TIME = re.compile(r"(\d{4})-(\d{2})-(\d{2})[T ](\d{2}):(\d{2}):(\d{2})(\.\d+)?")
_UT_TIME = re.compile(r"(?<![\d:T-])(\d{1,2}):(\d{2}):(\d{2})(\.\d+)?\s*(?:\(\s*)?UTC?\b")


# --------------------------------------------------------------------------
# Parsing and labels
# --------------------------------------------------------------------------

def parse_circular(msg: dict | bytes | str) -> dict[str, Any]:
    """Message of ``gcn.circulars`` (or of the archive) -> record. The
    submitter is kept without the e-mail address."""
    if isinstance(msg, (bytes, str)):
        msg = json.loads(msg)
    cid = int(msg["circularId"])
    created = msg.get("createdOn")
    t = (iso(datetime.fromtimestamp(float(created) / 1000, timezone.utc))
         if created is not None else None)
    return {"id": cid, "subject": " ".join(str(msg.get("subject") or "").split()),
            "event": (msg.get("eventId") or None), "event_types": list(msg.get("eventType") or []),
            "submitter": re.sub(r"\s*<[^>]*>", "", str(msg.get("submitter") or "")).strip(),
            "time": t, "version": msg.get("version"), "url": ARCHIVE_URL.format(id=cid),
            "body": str(msg.get("body") or "")}


def _words(text: str) -> set[str]:
    """Lower-case words of a subject, with the pairs joined by ``/`` or a
    space also as one word (``svom/vt``, ``einstein probe``)."""
    low = text.lower()
    toks = re.findall(r"[a-z0-9í][a-z0-9í.+-]*", low)
    toks = [t.strip(".-+") for t in toks]
    out = set(toks)
    out |= {f"{a} {b}" for a, b in zip(toks, toks[1:])}
    out |= set(re.findall(r"[a-z0-9í-]+/[a-z0-9í-]+", low))
    return out


def classify(subject: str, body: str = "") -> dict[str, Any]:
    """Labels of a circular: ``bands`` (see :data:`BANDS`; empty when no
    keyword is found), ``result`` of an optical report (``detection``,
    ``upper_limit``, or ``observation`` when neither the subject nor the text
    says), ``redshift`` when the subject or the text announces one."""
    words = _words(subject)
    bands = [b for b, keys in BANDS.items() if words & set(keys)]
    # "burst", "localization"... are gamma only when nothing more precise is there
    if "gamma" in bands and len(bands) > 1 and not (words & (set(BANDS["gamma"])
                                                           - {"burst", "localization", "localisation"})):
        bands.remove("gamma")
    order = ["optical", "xray", "radio", "gamma", "gw", "nu"]
    bands.sort(key=order.index)
    result = redshift = None
    if "optical" in bands:
        if _UPPER.search(subject) and not re.search(r"detections? and upper", subject, re.I):
            result = "upper_limit"
        elif _DETECT.search(subject):
            result = "detection"
        else:
            up, det = bool(_BODY_UPPER.search(body)), bool(_BODY_DETECT.search(body))
            result = "detection" if det and not up else "upper_limit" if up and not det \
                else "observation"
        if re.search(r"redshift|spectroscop|\bz\s*[=~≈]", subject, re.I):
            for text in (subject, body):
                m = next((m for m in _REDSHIFT.finditer(text)
                          if not _NOT_REDSHIFT.search(text[max(0, m.start() - 40):m.end()])), None)
                if m:
                    redshift = float(m.group(1))
                    break
    return {"bands": bands, "result": result, "redshift": redshift}


def event_date(event: Optional[str]) -> Optional[datetime]:
    """UTC day of an event name: ``GRB 261113A``, ``EP261115a``,
    ``IceCube-261113A`` -> 2026-11-13."""
    m = _EVENT_DATE.search(event or "")
    if not m:
        return None
    try:
        return datetime(2000 + int(m.group(1)), int(m.group(2)), int(m.group(3)), tzinfo=timezone.utc)
    except ValueError:
        return None


def times_in(text: str, day: Optional[datetime]) -> list[datetime]:
    """UTC times written in a text: full dates (``2026-11-13T05:40:35``), and
    ``05:48:03 UT`` on ``day`` (the date of the event name)."""
    out = []
    for m in _ISO_TIME.finditer(text):
        try:
            out.append(datetime(*(int(m.group(i)) for i in range(1, 7)), tzinfo=timezone.utc)
                       + timedelta(seconds=float(m.group(7) or 0)))
        except ValueError:
            pass
    if day is not None:
        for m in _UT_TIME.finditer(text):
            h, mi, s = int(m.group(1)), int(m.group(2)), int(m.group(3))
            if h < 24 and mi < 60 and s < 61:
                out.append(day + timedelta(hours=h, minutes=mi, seconds=s + float(m.group(4) or 0)))
    return out


def _gps_utc(gps: float) -> datetime:
    from astropy.time import Time
    return Time(float(gps), format="gps").utc.to_datetime(timezone.utc)


def _has_token(text: str, ident: str) -> bool:
    """``ident`` in ``text`` as a whole identifier (``816645341`` is found in
    ``trigger 816645341.32731`` and not in ``8166453410``)."""
    return re.search(r"(?<![A-Za-z0-9])" + re.escape(ident) + r"(?![A-Za-z0-9])", text) is not None


def _norm(name: Optional[str]) -> str:
    return "".join(c for c in str(name or "").upper() if c.isalnum())


# --------------------------------------------------------------------------
# Store
# --------------------------------------------------------------------------

def _dir(root: str | Path, linked: bool = True) -> Path:
    d = Path(root) / "circulars"
    return d if linked else d / "unlinked"


def _read(p: Path) -> Optional[dict]:
    try:
        return json.loads(p.read_text())
    except (OSError, ValueError):
        return None


def stored_circulars(root: str | Path, linked: bool = True) -> list[dict]:
    """Circulars of a store, oldest first (``linked=False``: those waiting in
    ``circulars/unlinked/``)."""
    out = [c for c in (_read(p) for p in _dir(root, linked).glob("*.json")) if c]
    return sorted(out, key=lambda c: c.get("id") or 0)


def event_names(root: str | Path) -> dict[str, list[dict]]:
    """Event name -> the GRB notices it designates, learnt from the circulars
    linked by a trigger number or a trigger time."""
    out: dict[str, list[dict]] = {}
    for c in stored_circulars(root):
        if not c.get("event"):
            continue
        for ln in c.get("links") or []:
            if ln["type"] == "grb" and ln["how"] in ("trigger number", "trigger time"):
                known = out.setdefault(_norm(c["event"]), [])
                if not any(k["id"] == ln["id"] for k in known):
                    known.append({"type": "grb", "id": ln["id"]})
    return out


def find_links(circ: dict, root: str | Path, names: Optional[dict] = None) -> list[dict]:
    """Alerts of the store the circular is about: ``[{"type": gw | nu |
    optical | grb, "id": <id in the store>, "how": ...}]`` (rules in the
    module documentation)."""
    root = Path(root)
    subject, body = circ.get("subject") or "", circ.get("body") or ""
    head = f"{subject} {circ.get('event') or ''}"
    text = f"{subject}\n{body}"
    links: list[dict] = []

    def add(kind, ident, how):
        if not any(ln["type"] == kind and ln["id"] == ident for ln in links):
            links.append({"type": kind, "id": ident, "how": how})

    known_gw = {p.stem for p in (root / "gw").glob("*.json")}
    known_gw |= {p.name for p in (root / "gw_versions").glob("*") if p.is_dir()}
    for sid in _SUPEREVENT.findall(head):
        if sid in known_gw or _safe(sid) in known_gw:
            add("gw", sid, "superevent name")
    nus = {_norm(p.stem): p.stem for p in (root / "nu").glob("*.json")}
    for m in _ICECUBE.finditer(head):
        key = _norm("IceCube" + m.group(1))
        if key in nus:
            add("nu", nus[key], "event name")
    from .listener import optical_files
    for p in optical_files(root):
        if len(p.stem) >= 5 and _has_token(text, p.stem):
            add("optical", p.stem, "name of the optical transient")

    day = event_date(circ.get("event") or subject)
    times = times_in(text, day)
    by_time = []
    for p in sorted((root / "grb").glob("*.json")):
        rec = _read(p) or {}
        g = rec.get("grb") or {}
        ident = str(g.get("graceid") or "")
        if len(ident) >= MIN_ID_LEN and _has_token(text, ident):
            add("grb", p.stem, "trigger number")
        elif g.get("gps") is not None and day is not None and times:
            t0 = _gps_utc(g["gps"])
            if t0.date() == day.date():
                off = min(abs((t - t0).total_seconds()) for t in times)
                if off <= TRIGGER_TIME_TOL_S:
                    by_time.append(p.stem)
    for stem in by_time:
        add("grb", stem, "trigger time")
    if not any(ln["type"] == "grb" for ln in links) and circ.get("event"):
        names = event_names(root) if names is None else names
        for k in names.get(_norm(circ["event"]), []):
            add(k["type"], k["id"], "event name")
    return links


def _prune_unlinked(root: Path, now: datetime) -> None:
    for p in _dir(root, False).glob("*.json"):
        c = _read(p) or {}
        try:
            old = (now - parse_iso(c["time"])).total_seconds() > UNLINKED_DAYS * 86400
        except (KeyError, TypeError, ValueError):
            old = True
        if old:
            p.unlink(missing_ok=True)


def handle_circular(msg: dict | bytes | str, store: Any, _relink: bool = True) -> list[dict]:
    """Store a circular with its labels and its links, and write it in the
    journal when it is linked for the first time or to other alerts than
    before (a corrected circular is sent again with the same number).
    Returns the links."""
    circ = parse_circular(msg) if not (isinstance(msg, dict) and "links" in msg) else dict(msg)
    root = Path(store.root)
    before = event_names(root)
    circ.update(classify(circ["subject"], circ["body"]))
    circ["links"] = find_links(circ, root, before)
    circ["received"] = iso(utcnow())
    name = f"{circ['id']}.json"
    linked_p, wait_p = _dir(root) / name, _dir(root, False) / name
    previous = _read(linked_p)
    if not circ["links"]:
        if previous is not None:                       # a correction never unlinks
            circ["links"] = previous.get("links") or []
        elif circ.get("event"):
            wait_p.parent.mkdir(parents=True, exist_ok=True)
            _write_atomic(wait_p, json.dumps(circ, indent=1))
        if not circ["links"]:
            return []
    linked_p.parent.mkdir(parents=True, exist_ok=True)
    _write_atomic(linked_p, json.dumps(circ, indent=1))
    wait_p.unlink(missing_ok=True)

    def key(ls):
        return sorted((ln["type"], ln["id"]) for ln in ls or [])
    if previous is None or key(previous.get("links")) != key(circ["links"]):
        what = ", ".join(f"{ln['type']} {ln['id']}" for ln in circ["links"])
        label = "/".join(circ["bands"]) or "unlabelled"
        if circ.get("result"):
            label += f" {circ['result'].replace('_', ' ')}"
        journal_of(store).log(
            "circular", f"GCN {circ['id']}: {circ['subject']} [{label}] -> {what}",
            optical=next((ln["id"] for ln in circ["links"] if ln["type"] == "optical"), None),
            superevent=next((ln["id"] for ln in circ["links"] if ln["type"] == "gw"), None),
            circular=circ["id"], subject=circ["subject"], event=circ.get("event"),
            submitter=circ["submitter"], published=circ["time"], bands=circ["bands"],
            result=circ["result"], redshift=circ["redshift"], links=circ["links"], url=circ["url"])
        log.info("circular %s linked to %s", circ["id"], what)
    # a name learnt from this circular: the earlier reports that only gave the name
    if _relink and circ.get("event") and _norm(circ["event"]) not in before \
            and _norm(circ["event"]) in event_names(root):
        for w in stored_circulars(root, linked=False):
            if _norm(w.get("event")) == _norm(circ["event"]):
                handle_circular(w, store, _relink=False)
        _prune_unlinked(root, utcnow())
    return circ["links"]


# --------------------------------------------------------------------------
# What the circulars say about one alert
# --------------------------------------------------------------------------

def summary(circs: Iterable[dict]) -> Optional[dict]:
    """Digest of the circulars linked to one alert: how many, in which bands,
    and for the optical: whether a search was reported, its best outcome
    (``detection`` over ``upper_limit`` over ``observation``), the time of
    the first report, the latest redshift announced."""
    circs = sorted(circs, key=lambda c: c.get("id") or 0)
    if not circs:
        return None
    bands: dict[str, int] = {}
    for c in circs:
        for b in c.get("bands") or []:
            bands[b] = bands.get(b, 0) + 1
    opt = [c for c in circs if "optical" in (c.get("bands") or [])]
    rank = {"detection": 0, "upper_limit": 1, "observation": 2}
    z = [c for c in circs if c.get("redshift") is not None]
    return {"n": len(circs), "ids": [c["id"] for c in circs], "bands": bands,
            "name": next((c["event"] for c in circs if c.get("event")), None),
            "optical": (min((c.get("result") or "observation" for c in opt), key=rank.get)
                        if opt else None),
            "n_optical": len(opt), "first_optical": opt[0].get("time") if opt else None,
            "redshift": z[-1]["redshift"] if z else None,
            "redshift_circular": z[-1]["id"] if z else None,
            "latest": circs[-1].get("time")}


# --------------------------------------------------------------------------
# Archive
# --------------------------------------------------------------------------

def fetch_circular(circular_id: int, timeout: float = 30.0) -> dict:
    """One circular from the archive, as the Kafka message."""
    from urllib.request import Request, urlopen
    req = Request(ARCHIVE_URL.format(id=int(circular_id)) + ".json",
                  headers={"User-Agent": "raven-kn", "Accept": "application/json"})
    with urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())


def expand_ids(specs: Iterable[str]) -> list[int]:
    """``["46261", "46270-46272"]`` -> ``[46261, 46270, 46271, 46272]``."""
    out: list[int] = []
    for s in specs:
        a, _, b = str(s).partition("-")
        out.extend(range(int(a), int(b or a) + 1))
    return out
