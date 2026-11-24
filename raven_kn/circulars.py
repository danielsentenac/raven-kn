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
   most follow-up reports only say ``GRB 261113A``;
7. an optical report about a burst is also linked to the optical transient
   attached to that burst: the afterglow it describes is that transient.

Nothing is changed in the optical transients: a redshift announced in a
circular is recorded with the circular, not applied.

**Distance gate.** A redshift announced for a GRB gives its luminosity
distance (:func:`distance_gate`, Planck 2018), compared with the distance
beyond which a search is not worth it (:data:`D_MAX_USEFUL_MPC`, the default
of :func:`raven_kn.search.assess`). It is written in the journal
(``grb_distance``) and given with the digest of the alert; the service makes
a GRB within reach a candidate and marks the others as too far. GRB 261114A
at z = 3.716 is at 33.6 Gpc: no GW search can be about it.

**Optical transient of a reported counterpart.** A circular labelled as an
optical detection, linked to one GRB or X-ray notice and giving a position
(:func:`position_in`) creates the optical transient of that burst when it has
none: source ``GCN circular <number>``, the burst attached, so that its
merger-time window is the one of the GRB. It has no photometry. A redshift
that a circular gives as spectroscopic is applied to the transient attached to
the burst, unless it already has a spectroscopic redshift from another source.
Reports about stellar flares and other Galactic sources create nothing.
Each time such a transient is created, gets another redshift or gets accepted
photometry, the search decision is run again (:func:`assess_transient`), so
that a revised redshift changes the decision without anybody asking.

**Reported photometry.** The magnitudes written as ``r = 19.95 +/- 0.06`` or
``r > 21.56`` are extracted (:func:`photometry_in`) with the time found near
them (a delay after the trigger, else a date), and kept with the circular as
*reported, unverified*: they are shown, with the sentence they come from, and
are not part of the photometry of the transient, so that they do not enter the
kilonova test. ``raven-kn accept-photometry`` copies the chosen ones into the
transient, once a person has read them.

``raven-kn circular 46261 46270-46290 --store DIR`` reads circulars from the
archive (https://gcn.nasa.gov/circulars/<number>.json), to fill what was sent
before the listener subscribed, or to link again.
"""

from __future__ import annotations

import json
import logging
import os
import re
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Iterable, Optional

from .circular_facts import facts_in, links_in
from .journal import _safe, _write_atomic, iso, journal_of, parse_iso, utcnow

log = logging.getLogger("raven_kn.circulars")

CIRCULARS_TOPIC = "gcn.circulars"
ARCHIVE_URL = "https://gcn.nasa.gov/circulars/{id}"
TRIGGER_TIME_TOL_S = 1.5         # the circulars give the trigger time to the second
UNLINKED_DAYS = 30.0
MIN_ID_LEN = 6                   # shorter identifiers are not searched for in the text
D_MAX_USEFUL_MPC = 500.0         # as d_max_useful_mpc of raven_kn.search.assess

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


def distance_gate(z: float, circular: Optional[int] = None,
                  d_max_mpc: float = D_MAX_USEFUL_MPC) -> dict[str, Any]:
    """Luminosity distance of a redshift (Planck 2018) and whether a GW
    search can reach it: ``{"z", "circular", "d_mpc", "limit_mpc",
    "within_reach", "text"}``."""
    from astropy.cosmology import Planck18
    d = float(Planck18.luminosity_distance(float(z)).to_value("Mpc"))
    within = d <= d_max_mpc
    dist = f"{d / 1000:.3g} Gpc" if d >= 1000 else f"{d:.3g} Mpc"
    text = (f"redshift {float(z):g}" + (f" (GCN {circular})" if circular else "")
            + f", luminosity distance {dist}: "
            + (f"within the range of a GW search ({d_max_mpc:.0f} Mpc)" if within
               else f"too far for a GW search (limit {d_max_mpc:.0f} Mpc)"))
    return {"z": float(z), "circular": circular, "d_mpc": d, "limit_mpc": float(d_max_mpc),
            "within_reach": within, "text": text}


# --------------------------------------------------------------------------
# Position and photometry written in the text
# --------------------------------------------------------------------------

_RA_SEX = re.compile(r"R\.?\s?A\.?\s*(?:\(?J?2000(?:\.0)?\)?)?\s*[=:]\s*"
                     r"(\d{1,2})\s*[h: ]\s*(\d{1,2})\s*[m: ]\s*(\d{1,2}(?:\.\d+)?)", re.I)
_DEC_SEX = re.compile(r"Dec(?:l)?\.?\s*(?:\(?J?2000(?:\.0)?\)?)?\s*[=:]\s*"
                      r"([+\-−–]?)\s*(\d{1,2})\s*[d°: ]\s*(\d{1,2})\s*[m'′: ]\s*(\d{1,2}(?:\.\d+)?)",
                      re.I)
_RADEC_DEG = re.compile(r"R\.?A\.?\s*,\s*Dec(?:l)?\.?\s*(?:\(?J?2000(?:\.0)?\)?)?\s*[=:]\s*"
                        r"(\d{1,3}\.\d+)\s*,\s*([+\-−–]?\d{1,2}\.\d+)", re.I)
_RA_DEG = re.compile(r"R\.?A\.?\s*(?:\(?J?2000(?:\.0)?\)?)?\s*[=:]\s*(\d{1,3}\.\d+)\s*(?:deg|d\b|°)?"
                     r"[\s,;]*Dec(?:l)?\.?\s*(?:\(?J?2000(?:\.0)?\)?)?\s*[=:]\s*([+\-−–]?\d{1,2}\.\d+)",
                     re.I)
_POS_ERR = re.compile(r"(?:uncertaint\w+|error|accuracy|radius)\s*(?:of|is|:)?\s*(?:about|~|±|\+/-)?\s*"
                      r"(\d+(?:\.\d+)?)\s*(arcsec|\"|''|arcmin|')", re.I)
_GALACTIC = re.compile(r"flar(?:e|ing) star|stellar flare|flare of|known (?:AGN|star|source)|"
                       r"cataclysmic|dwarf nova|\bCV\b|Galactic (?:source|transient|origin)|"
                       r"M[- ]dwarf|not (?:a|an) (?:GRB|extragalactic)", re.I)
_AT_NAME = re.compile(r"\b(AT|SN)\s?(20\d{2}[a-z]{2,5})\b")
_SPECTROSCOPIC = re.compile(r"spectroscop|spectrum|spectra\b|absorption|emission line", re.I)

_BANDS = {**{b: f"ps1::{b}" for b in "grizy"}, "u": "sdssu", "J": "2massj", "H": "2massh",
          "K": "2massks", "Ks": "2massks"}
_MAG = re.compile(r"(?<![A-Za-z0-9_.])(Ks|[ugrizyUBVRIJHKL])['’′]?(?:[- ]band)?"
                  r"(?:\s*\((?:AB|Vega)\))?\s*(=|~|≈|>=|≥|>|<)\s*(\d{1,2}\.\d{1,3})(?!\d)"
                  r"(?:\s*(?:\+/-|±|\+-)\s*(\d?\.\d{1,3}))?")
_UNITS = {"ks": 1000.0, "ksec": 1000.0, "s": 1.0, "sec": 1.0, "secs": 1.0, "second": 1.0, "seconds": 1.0,
          "min": 60.0, "mins": 60.0, "minute": 60.0, "minutes": 60.0,
          "h": 3600.0, "hr": 3600.0, "hrs": 3600.0, "hour": 3600.0, "hours": 3600.0,
          "d": 86400.0, "day": 86400.0, "days": 86400.0}
_U = r"(ksec|ks|seconds?|secs?|s|minutes?|mins?|hours?|hrs?|h|days?|d)"
_AFTER = r"\s*(?:after|post|since|from)[- ]?\s*(?:the\s+)?(?:[\w/-]+\s+){0,3}?(?:trigger|burst|T0|GRB|event|alert|notice)"
_DELAY_RANGE = re.compile(r"(?:from\s+)?~?(\d+(?:\.\d+)?)\s*(?:to|-|–)\s*~?(\d+(?:\.\d+)?)\s*" + _U + r"\b"
                          + _AFTER, re.I)
_DELAY = re.compile(r"(?<![\d.])~?\+?(\d+(?:\.\d+)?)\s*" + _U + r"\b" + _AFTER, re.I)
# "3 ks after the GWAC discovery time": counted from the time the text gave for that discovery
_DELAY_FOUND = re.compile(r"(?<![\d.])~?\+?(\d+(?:\.\d+)?)\s*" + _U + r"\b\s*(?:after|since|from)\s+(?:the\s+)?"
                          r"(?:[\w/-]+\s+){0,3}?(?:discovery|detection|first detection)", re.I)
_DELAY_PLUS = re.compile(r"(?:\(|T0?\s*)\+\s*(\d+(?:\.\d+)?)\s*" + _U + r"\b")     # "(+0.11h)", "T0+3.2 h"


def position_in(text: str) -> Optional[dict[str, float]]:
    """First position written in a text: ``{"ra", "dec"}`` in degrees
    (sexagesimal ``RA (J2000) = 22h29m24.76s``, ``Dec = +27d01m34.4s``, or
    decimal degrees), with ``err_arcsec`` when an uncertainty follows it
    (else 1 arcsec)."""
    ra = dec = None
    end = 0
    m, d = _RA_SEX.search(text), None
    if m:
        d = _DEC_SEX.search(text, m.end())
    if m and d and d.start() - m.end() < 200:
        ra = 15 * (int(m.group(1)) + int(m.group(2)) / 60 + float(m.group(3)) / 3600)
        dec = int(d.group(2)) + int(d.group(3)) / 60 + float(d.group(4)) / 3600
        if d.group(1) in ("-", "−", "–"):
            dec = -dec
        end = d.end()
    else:
        m = _RADEC_DEG.search(text) or _RA_DEG.search(text)
        if m:
            ra = float(m.group(1))
            dec = float(m.group(2).replace("−", "-").replace("–", "-"))
            end = m.end()
    if ra is None or not (0 <= ra < 360 and -90 <= dec <= 90):
        return None
    e = _POS_ERR.search(text, end, end + 250)
    err = 1.0
    if e:
        err = float(e.group(1)) * (60.0 if e.group(2) in ("arcmin", "'") else 1.0)
    return {"ra": round(ra, 6), "dec": round(dec, 6), "err_arcsec": err}


def _segments(body: str) -> list[tuple[int, str]]:
    """Sentences of a text with their offset (paragraphs are often one line)."""
    out, pos = [], 0
    for part in re.split(r"(\n\s*\n|(?<=[.;:])\s+(?=[A-Z(]))", body):
        if part.strip() and not re.fullmatch(r"\s+", part):
            out.append((pos, part))
        pos += len(part)
    return out


def photometry_in(body: str, trigger: Optional[datetime] = None) -> list[dict[str, Any]]:
    """Magnitudes written in a text as ``r = 19.95 +/- 0.06`` or ``r > 21.56``.

    Each one: ``band`` (as written), ``mag``, ``mag_err``, ``upper_limit``
    (``>``, or a sentence about a limit), ``system`` (``AB``, ``Vega`` or
    None), ``time`` (ISO, or None) with ``time_from``: ``delay after the
    trigger`` (the last ``N hours after the trigger`` before the magnitude,
    in its sentence or the three before; needs ``trigger``), ``delay after
    the discovery time of the text`` (``3 ks after the GWAC discovery time``,
    from the last time given before the sentence) or ``date in the text``; and ``text``, the sentence it comes from. Magnitudes given in a
    table are not read. ``z = 2.4`` is a
    redshift, not a z-band magnitude: magnitudes are between 8 and 30."""
    segs = _segments(body)
    out, seen = [], set()
    for k, (off, seg) in enumerate(segs):
        if seg.count("|") >= 2:                        # a table: columns cannot be told apart
            continue
        for m in _MAG.finditer(seg):
            band, op, mag = m.group(1), m.group(2), float(m.group(3))
            if not 8.0 <= mag <= 30.0 or (band == "z" and mag < 12.0):
                continue
            err = float(m.group(4)) if m.group(4) else None
            limit = op in (">", ">=", "≥") or (err is None and bool(
                re.search(r"upper limit|limiting mag|not detect|no \w+ ?(?:source|counterpart)", seg, re.I)))
            if op == "<":                              # "brighter than": not a measurement
                continue
            # the time: before the magnitude in its sentence, else after it, else the sentence before
            before = seg[:m.start()]
            when, how = None, None
            # a delay after the discovery: from the last time the text gave before this sentence
            found = [h for h in _DELAY_FOUND.finditer(seg)]
            given = times_in(body[:off], None) if found else []
            if found and given:
                h = min(found, key=lambda h: abs(h.start() - m.start()))
                when = given[-1] + timedelta(seconds=float(h.group(1)) * _UNITS[h.group(2).lower()])
                how = "delay after the discovery time of the text"
            for text, last in (() if when else ((before, True), (seg[m.end():], False),
                               *((segs[j][1], True) for j in range(k - 1, max(k - 4, -1), -1)))):
                if trigger is not None:
                    hits = [(h.start(), (float(h.group(1)) + float(h.group(2))) / 2 * _UNITS[h.group(3).lower()])
                            for h in _DELAY_RANGE.finditer(text)]
                    hits += [(h.start(), float(h.group(1)) * _UNITS[h.group(2).lower()])
                             for rx in (_DELAY, _DELAY_PLUS) for h in rx.finditer(text)]
                    if hits:
                        hits.sort()
                        when = trigger + timedelta(seconds=(hits[-1] if last else hits[0])[1])
                        how = "delay after the trigger"
                        break
                dates = times_in(text, None)
                if dates:
                    when, how = (dates[-1] if last else dates[0]), "date in the text"
                    break
            key = (band, mag, err, limit)
            if key in seen:
                continue
            seen.add(key)
            ctx = seg[max(0, m.start() - 200):m.end() + 120]
            out.append({"band": band, "mag": mag, "mag_err": err, "upper_limit": limit,
                        "system": "Vega" if re.search(r"\bVega\b", ctx) else
                                  "AB" if re.search(r"\bAB\b", ctx) else None,
                        "time": iso(when) if when else None, "time_from": how,
                        "text": " ".join(seg.split())[:300]})
    return out


def pipeline_band(band: str) -> str:
    """Band name of the pipeline for a letter of a circular (``r`` ->
    ``ps1::r``, ``J`` -> ``2massj``); other letters are kept as written."""
    return _BANDS.get(band, band)


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
    circ["distance_logged"] = bool((_read(_dir(root) / f"{circ['id']}.json") or {}).get("distance_logged"))
    # what an optical report says: where, how bright, and whether its redshift is from a spectrum
    text = f"{circ['subject']}\n{circ['body']}"
    optical = "optical" in circ["bands"]
    circ["position"] = position_in(circ["body"]) if optical and circ["result"] == "detection" else None
    circ["photometry"] = photometry_in(circ["body"], _trigger_time(root, circ["links"])) if optical else []
    circ["redshift_kind"] = (None if circ["redshift"] is None else
                             "spectroscopic" if _SPECTROSCOPIC.search(text) else "unknown")
    circ["galactic"] = bool(_GALACTIC.search(text))
    # the quantities and the web addresses of the text (reported, not verified)
    circ["facts"] = facts_in(circ["subject"], circ["body"], circ["bands"])
    circ["web_links"] = links_in(circ["body"])
    _link_transient(circ, store)
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
    # the distance of the GRB, once per circular (also for one stored before this rule)
    grbs = [ln["id"] for ln in circ["links"] if ln["type"] == "grb"]
    if circ["redshift"] is not None and grbs and not (previous or {}).get("distance_logged"):
        gate = distance_gate(circ["redshift"], circ["id"])
        journal_of(store).log(
            "grb_distance", f"{circ.get('event') or 'GRB'} ({', '.join(grbs)}): {gate['text']}",
            level="info" if not gate["within_reach"] else "warning",
            event=circ.get("event"), triggers=grbs, **gate)
        circ["distance_logged"] = True
        _write_atomic(linked_p, json.dumps(circ, indent=1))
    try:
        _counterpart(circ, store)
    except Exception:                                  # the circular itself is stored and written
        log.exception("circular %s: optical transient not created or updated", circ["id"])
    # a name learnt from this circular: the earlier reports that only gave the name
    if _relink and circ.get("event") and _norm(circ["event"]) not in before \
            and _norm(circ["event"]) in event_names(root):
        for w in stored_circulars(root, linked=False):
            if _norm(w.get("event")) == _norm(circ["event"]):
                handle_circular(w, store, _relink=False)
        _prune_unlinked(root, utcnow())
    return circ["links"]


# --------------------------------------------------------------------------
# Optical transient of a reported counterpart, reported photometry
# --------------------------------------------------------------------------

def _grb_record(root: Path, links: Iterable[dict]) -> Optional[tuple[str, dict]]:
    """(id in the store, stored record) of the burst a circular is about,
    when it is about one burst only and its notice was received."""
    grbs = [ln["id"] for ln in links if ln["type"] == "grb"]
    if len(grbs) != 1:
        return None
    rec = _read(root / "grb" / f"{grbs[0]}.json") or {}
    return (grbs[0], rec) if rec.get("grb") else None


def _trigger_time(root: Path, links: Iterable[dict]) -> Optional[datetime]:
    got = _grb_record(root, links)
    return _gps_utc(got[1]["grb"]["gps"]) if got else None


def _attached_to(store: Any, g: dict):
    """The optical transient the burst ``g`` (stored ``grb`` record) is attached to."""
    for kn in store.transients():
        a = kn.grb
        if a is None or a.instrument != g.get("instrument"):
            continue
        if (a.graceid is not None and a.graceid == g.get("graceid")) or abs(a.gps - g["gps"]) < 1.0:
            return kn
    return None


def _spectroscopic_redshift(root: Path, grb_id: str) -> Optional[dict]:
    """Latest circular linked to the burst that gives a spectroscopic redshift."""
    z = [c for c in stored_circulars(root)
         if c.get("redshift") is not None and c.get("redshift_kind") == "spectroscopic"
         and any(ln["type"] == "grb" and ln["id"] == grb_id for ln in c.get("links") or [])]
    return z[-1] if z else None


def _link_transient(circ: dict, store: Any, kn: Any = None) -> bool:
    """Link an optical report about one burst to the optical transient
    attached to that burst (``kn``, else looked up). True if a link was added."""
    if "optical" not in (circ.get("bands") or []) or not hasattr(store, "transients"):
        return False
    got = _grb_record(Path(store.root), circ.get("links") or [])
    if got is None:
        return False
    if kn is None:
        kn = _attached_to(store, got[1]["grb"])
    elif kn.grb is None or _attached_to_same(kn, got[1]["grb"]) is False:
        return False
    if kn is None or any(ln["type"] == "optical" and ln["id"] == kn.name for ln in circ["links"]):
        return False
    circ["links"].append({"type": "optical", "id": kn.name, "how": "optical report about its burst"})
    return True


def _attached_to_same(kn: Any, g: dict) -> bool:
    a = kn.grb
    return a.instrument == g.get("instrument") and (
        (a.graceid is not None and a.graceid == g.get("graceid")) or abs(a.gps - g["gps"]) < 1.0)


def assess_transient(kn: Any, store: Any, **assess_kw) -> Optional[dict]:
    """Run the search decision of an optical transient
    (:func:`raven_kn.search.assess`) into ``<store>/assess/`` and write it in
    the journal when it differs from the previous one. The 3D sky map (100 MB)
    is kept only when the search is triggered: the request refers to it. A
    failure is written in the journal and returns None."""
    from .journal import log_decision
    adir = Path(store.root) / "assess"
    stem = kn.name.replace(" ", "_")
    previous = _read(adir / f"{stem}_decision.json")
    try:
        from .search import assess
        d = assess(kn, out_dir=adir, **assess_kw)
    except Exception as e:
        log.exception("%s: search decision failed", kn.name)
        journal_of(store).log("decision_failed", f"{kn.name}: search decision failed "
                              f"({type(e).__name__}: {e})"[:300], level="error", optical=kn.name,
                              error=repr(e)[:500])
        return None
    if not d["trigger_targeted_search"]:
        (adir / f"{stem}_skymap.fits").unlink(missing_ok=True)
    log_decision(journal_of(store), d, previous)
    try:                                               # what the brokers have at its position
        from .brokers import REFRESH_H, crossmatch_transient
        if os.environ.get("RAVEN_KN_BROKERS", "1") != "0":   # "0": no request to the brokers (tests)
            crossmatch_transient(kn, store.root, max_age_h=REFRESH_H)
    except Exception:
        log.exception("%s: brokers not asked", kn.name)
    return d


def _counterpart(circ: dict, store: Any) -> Optional[str]:
    """Create the optical transient of the burst from a reported counterpart
    with a position, and apply a spectroscopic redshift to the transient
    attached to the burst (rules in the module documentation). Returns the
    name of the transient created or updated."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from .connection import GRBAssociation, OpticalTransient, Redshift
    from .listener import add_optical
    root = Path(store.root)
    got = _grb_record(root, circ["links"])
    if got is None or not hasattr(store, "transients"):
        return None
    grb_id, rec = got
    g = rec["grb"]
    kn = _attached_to(store, g)
    changed = False
    pos = circ.get("position")
    if kn is None and pos and circ.get("result") == "detection" and not circ.get("galactic"):
        x = rec.get("xray") or {}
        ra0, dec0, err0 = ((x.get("ra"), x.get("dec"), x.get("err_deg")) if x.get("ra") is not None
                           else (rec.get("ra"), rec.get("dec"), rec.get("err_deg")))
        if ra0 is not None:
            sep = SkyCoord(ra0 * u.deg, dec0 * u.deg).separation(
                SkyCoord(pos["ra"] * u.deg, pos["dec"] * u.deg)).deg
            if sep > 2.0 * (err0 or 0.0) + 1.0 / 60:
                log.info("circular %s: position %.3f deg from the burst, outside its error circle: "
                         "no optical transient", circ["id"], sep)
                return None
        at = _AT_NAME.search(f"{circ['subject']}\n{circ['body']}")
        base = circ.get("event") or f"GCN{circ['id']}"
        name = f"{at.group(1)}{at.group(2)}" if at else _safe("".join(base.split())) + "-OT"
        if any(t.name == name for t in store.transients()):
            return None                                # a transient of that name, attached elsewhere
        kn = OpticalTransient(
            name=name, ra=pos["ra"], dec=pos["dec"], pos_err_arcsec=pos["err_arcsec"],
            source=f"GCN circular {circ['id']}", grb=GRBAssociation(**g),
            notes=f"Optical counterpart of {circ.get('event') or 'the burst'} reported in GCN circular "
                  f"{circ['id']} ({circ['subject']}); position read from its text. "
                  f"GRB {g.get('graceid') or ''} ({g.get('instrument')}) attached from the circular.")
        changed = True
    if kn is None:
        return None
    zc = _spectroscopic_redshift(root, grb_id)
    if zc is not None:
        src = f"GCN circular {zc['id']}"
        old = kn.redshift
        mine = old is not None and old.source.startswith("GCN circular")
        if (old is None or old.kind != "spectroscopic" or mine) \
                and not (mine and abs(old.z - float(zc["redshift"])) < 1e-9):   # the same value again
            kn.redshift = Redshift(z=float(zc["redshift"]), z_err=1e-3, kind="spectroscopic", source=src)
            kn.notes = (kn.notes + f" Redshift {zc['redshift']:g} from {src}.").strip()
            changed = True
    if changed:
        add_optical(kn, store)
        # the optical reports about the burst, this one included, are about this transient
        for c in [circ, *stored_circulars(root)]:
            if _link_transient(c, store, kn):
                _write_atomic(_dir(root) / f"{c['id']}.json", json.dumps(c, indent=1))
        assess_transient(kn, store)
        log.info("circular %s: optical transient %s %s", circ["id"], kn.name,
                 "created" if kn.source == f"GCN circular {circ['id']}" else "updated")
        return kn.name
    return None


def photometry_of(root: str | Path, c: dict) -> list[dict]:
    """Magnitudes of a stored circular, read again from its text (the rules
    change, the stored text does not); those stored with it when it has no text."""
    if not c.get("body"):
        return c.get("photometry") or []
    if "optical" not in (c.get("bands") or []):
        return []
    return photometry_in(c["body"], _trigger_time(Path(root), c.get("links") or []))


def alert_keys(root: str | Path, kn: Any) -> set[tuple[str, str]]:
    """``(type, id)`` of the alerts whose circulars are about an optical
    transient: itself and the stored notices of its burst."""
    root = Path(root)
    keys = {("optical", kn.name)}
    if kn.grb is not None:
        for p in (root / "grb").glob("*.json"):
            g = (_read(p) or {}).get("grb") or {}
            if g.get("instrument") == kn.grb.instrument and (
                    (kn.grb.graceid is not None and g.get("graceid") == kn.grb.graceid)
                    or abs((g.get("gps") or 0) - kn.grb.gps) < 1.0):
                keys.add(("grb", p.stem))
    return keys


def reported(root: str | Path, keys: Iterable[tuple[str, str]]) -> dict[str, Any]:
    """What the circulars linked to the alerts ``keys`` report, oldest
    circular first, by category: ``prompt``, ``afterglow``, ``spectroscopy``
    (the facts of :func:`raven_kn.circular_facts.facts_in`, and the redshift
    announced), ``photometry`` (:func:`photometry_in`) and ``links``
    (:func:`raven_kn.circular_facts.links_in`). Each item carries
    ``circular``, ``subject``, ``submitter``, ``published`` and ``circular_url``.
    Nothing here is verified. The text is read again at each call."""
    keys = set(keys)
    out: dict[str, Any] = {"prompt": [], "afterglow": [], "spectroscopy": [], "photometry": [],
                           "links": [], "circulars": []}
    for c in stored_circulars(root):
        if not any((ln["type"], ln["id"]) in keys for ln in c.get("links") or []):
            continue
        src = {"circular": c["id"], "subject": c.get("subject"), "submitter": c.get("submitter"),
               "published": c.get("time"), "circular_url": c.get("url")}
        out["circulars"].append(c["id"])
        # read again from the text: the rules change, the stored text does not
        facts = facts_in(c.get("subject") or "", c.get("body") or "", c.get("bands"))
        if c.get("redshift") is not None:
            out["spectroscopy"].append({
                "category": "spectroscopy", "kind": "redshift", "value": c["redshift"],
                "err_minus": None, "err_plus": None, "approx": False, "unit": "", "band": None,
                "instrument": None, "limit": False, "what": c.get("redshift_kind"),
                "text": c.get("subject") or "", **src})
        for f in facts:
            if f.get("category") in ("prompt", "afterglow", "spectroscopy"):
                out[f["category"]].append({**f, **src})
        for i, pt in enumerate(photometry_of(root, c)):
            out["photometry"].append({**pt, "index": i, **src})
        for ln in links_in(c.get("body") or ""):
            out["links"].append({**ln, "bands": c.get("bands") or [], **src})
    return out


def reported_photometry(root: str | Path, kn: Any) -> list[dict]:
    """Magnitudes reported in the circulars about an optical transient or
    about its burst, oldest circular first: the points of
    :func:`photometry_in`, each with ``circular``, ``index`` (its place in
    that circular) and ``accepted`` (already copied into the transient)."""
    root = Path(root)
    keys = alert_keys(root, kn)
    have = {(p.survey, p.band, round(p.mag, 3)) for p in kn.photometry}
    out = []
    for c in stored_circulars(root):
        if not any((ln["type"], ln["id"]) in keys for ln in c.get("links") or []):
            continue
        for i, pt in enumerate(photometry_of(root, c)):
            out.append({**pt, "circular": c["id"], "index": i, "subject": c["subject"],
                        "accepted": (f"GCN {c['id']}", pipeline_band(pt["band"]), round(pt["mag"], 3)) in have})
    return out


def accept_photometry(store: Any, name: str, circulars: Optional[Iterable[int]] = None) -> list[dict]:
    """Copy reported magnitudes into the photometry of the optical transient
    ``name`` (survey ``GCN <number>``): those of the given circulars, or all.
    Points without a time cannot be copied. Returns the points added."""
    from astropy.time import Time
    from .connection import PhotometryPoint
    from .listener import add_optical
    kn = next((t for t in store.transients() if t.name == name), None)
    if kn is None:
        raise ValueError(f"no optical transient {name} in the store")
    want = set(int(c) for c in circulars) if circulars else None
    added = []
    for pt in reported_photometry(store.root, kn):
        if pt["accepted"] or pt["time"] is None or (want is not None and pt["circular"] not in want):
            continue
        kn.photometry.append(PhotometryPoint(
            mjd=float(Time(parse_iso(pt["time"])).mjd), band=pipeline_band(pt["band"]), mag=pt["mag"],
            mag_err=pt["mag_err"], upper_limit=bool(pt["upper_limit"]), survey=f"GCN {pt['circular']}"))
        added.append(pt)
    if added:
        add_optical(kn, store)
        assess_transient(kn, store)
    return added


# --------------------------------------------------------------------------
# What the circulars say about one alert
# --------------------------------------------------------------------------

def summary(circs: Iterable[dict]) -> Optional[dict]:
    """Digest of the circulars linked to one alert: how many, in which bands,
    and for the optical: whether a search was reported, its best outcome
    (``detection`` over ``upper_limit`` over ``observation``), the time of
    the first report, the latest redshift announced and its distance gate
    (:func:`distance_gate`)."""
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
            "distance": distance_gate(z[-1]["redshift"], z[-1]["id"]) if z else None,
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
