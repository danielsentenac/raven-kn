"""Quantities and links that observers write in GCN Circulars.

:func:`raven_kn.circulars.classify` labels a circular; this module reads the
numbers of its text, by keyword rules written from the circulars of November
2026. Each one is a *fact*: what was *reported*, with the sentence it comes
from, never verified and never used by the pipeline.

Categories and kinds:

``prompt``        ``t90`` (s), ``fluence`` (erg/cm^2), ``peak_flux``
                  (erg/cm^2/s or ph/cm^2/s), ``epeak`` (keV), ``epeak_rest``,
                  ``photon_index``, ``eiso`` (erg), ``liso`` (erg/s),
                  ``burst_class`` (short or long, as written)
``afterglow``     ``decay_index`` (flux ~ t^index), ``xray_flux``
                  (erg/cm^2/s), ``nh`` (cm^-2), ``radio_flux`` (mJy or uJy)
``spectroscopy``  ``lines`` (the species named), ``host`` (a sentence about
                  the host galaxy)

A fact: ``category``, ``kind``, ``value`` (None for ``lines``, ``host``,
``burst_class``), ``err_minus``, ``err_plus`` (positive numbers, None when not
given), ``unit``, ``band`` (the energy range or the frequency given in the
sentence), ``instrument`` (named in the sentence, else in the subject),
``approx`` (written with ``~`` or ``about``), ``limit`` (the sentence says an
upper limit), ``what`` (text of a fact without a number) and ``text``, the
sentence.

Uncertainties are read as written: ``7.2 +/-0.2``, ``41 +3/-1``,
``46 -9 / +6``, ``7.447 (-0.084 +0.121)``, ``106(-9,+10)``,
``(2.05 +0.15/-0.20) E-05``, ``(3.51 ± 0.61)x10^-6``. A value in a table is
not read.

:func:`links_in` lists the web addresses of the text (references to other
circulars excepted) with the words that introduce them and a kind: ``image``
(shown by the application), ``data`` or ``page``.
"""

from __future__ import annotations

import re
from typing import Any, Optional

_N = r"\d+(?:\.\d+)?"
_SIGN = r"[-+−–]"
_EXP = (r"(?:\s*\)?\s*(?:[x×*]\s*10\s*\^?\s*\(?\s*(?P<e1>" + _SIGN + r"?\s*\d+)\s*\)?"
        r"|\s?[eE]\s?(?P<e2>" + _SIGN + r"?\d+)(?![\w.])))?")
# a value with its uncertainty: "(2.05 +0.15/-0.20) E-05", "106(-9,+10)", "7.2 +/-0.2", "~0.2"
_QTY = re.compile(
    r"(?P<approx>~|≈|about\s+|around\s+|approximately\s+)?\(?\s*(?P<v>" + _SIGN + r"?\s?" + _N + r")"
    r"(?:\s*(?:\+\s*/\s*-|±|\+-)\s*(?P<sym>" + _N + r")"
    r"|\s*\(?\s*(?P<s1>" + _SIGN + r")\s*(?P<a1>" + _N + r")\s*[/,]?\s*(?P<s2>" + _SIGN + r")\s*(?P<a2>" + _N + r")\s*\)?"
    r")?\s*\)?" + _EXP)
_RANGE = re.compile(r"(?:>\s*)?" + _N + r"\s*(?:keV|MeV|GeV)?\s*(?:-|–|to)\s*" + _N + r"\s*(?:keV|MeV|GeV)"
                    r"|>\s*" + _N + r"\s*(?:keV|MeV|GeV)", re.I)
_FREQ = re.compile(_N + r"\s*(?:GHz|MHz)")
_GLUE = (r"(?:\s*(?:,|\)|=|:|-value\b|(?:duration|is|was|of|to|with|a|an|the|value|measured|estimated)\b))*\s*")

_INSTRUMENTS = re.compile(
    r"SVOM/(?:ECLAIRs|GRM|VT|MXT)|ECLAIRs|Fermi[- /]?GBM|Fermi[- /]?LAT|\bGBM\b|\bLAT\b|Swift[-/ ]?(?:BAT|XRT|UVOT)"
    r"|\bBAT\b|\bXRT\b|\bUVOT\b|Konus-Wind|GECAM-[ABCD]|Insight-HXMT(?:/HE)?|HXMT(?:/HE)?|AstroSat[ -]CZTI"
    r"|CALET|AGILE|INTEGRAL|EP[-/ ]?(?:WXT|FXT)|\bWXT\b|\bFXT\b|GRBAlpha|Glowbug|MAXI|NICER|NuSTAR|Chandra"
    r"|XMM-Newton|ATCA|\bJ?VLA\b|MeerKAT|ALMA|NOEMA|u?GMRT|e-MERLIN|AMI-LA|VLBA|\bEVN\b"
    r"|VLT/[A-Za-z0-9-]+|X-shooter|FORS2|\bGTC\b|\bNOT\b|Keck|Gemini|GWAC(?:-F60)?|GROND|\bLBT\b|\bTNG\b")

_LABELS = {  # kind -> (category, anchor, unit that must follow or None, unit when none is written)
    "t90": ("prompt", r"T[_ ]?90", None, "s"),
    "fluence": ("prompt", r"fluence", r"erg\s*/?\s*cm", "erg/cm^2"),
    "peak_flux": ("prompt", r"peak (?:energy |photon )?flux|photon flux",
                  r"(?:erg|ph(?:otons?)?)\s*[ /]?\s*(?:cm\^?-?2\s*[ /]?\s*s|s\s*[ /]?\s*cm)", None),
    "epeak_rest": ("prompt", r"E_?p(?:eak)?,\s?z|rest[- ]frame peak energy(?! flux)[^.;]{0,60}?E_?p(?:eak)?,\s?z"
                   r"|rest[- ]frame peak energy(?! flux)", None, "keV"),
    "epeak": ("prompt", r"E_?peak\b|\bEp\b(?!,\s?z)|(?<!rest-frame )(?<!rest frame )peak energy(?! flux)", None, "keV"),
    "photon_index": ("prompt", r"(?:power[- ]?law|photon) index(?: alpha| beta)?|(?<!\^)\balpha\b(?!\s*\))",
                     None, ""),
    "eiso": ("prompt", r"E_?iso\b", r"erg(?!\s*/)", "erg"),
    "liso": ("prompt", r"L_?iso\b", r"erg\s*/?\s*s", "erg/s"),
    "decay_index": ("afterglow", r"(?:decay|temporal) index|decay with an? index|slope(?: index)?|\bt\s?\^", None, ""),
    "nh": ("afterglow", r"column density(?: \(N_?H\))?|\bN_?H\b", r"cm\s?\^?\s?-2", "cm^-2"),
}
_FLUX_UNIT = re.compile(r"erg\s*[ /]?\s*(?:cm\^?-?2\s*[ /]?\s*s|s\s*[ /]?\s*cm)", re.I)
_RADIO_UNIT = re.compile(r"(mJy|[uμµ]Jy)")
_CLASS = re.compile(r"\b(?:(likely|possible|probable)\s+)?(short|long)[- ](?:duration[- ])?(?:hard[- ])?"
                    r"(?:GRB|burst|gamma[- ]ray burst)", re.I)
_SPECIES = re.compile(
    r"Ly(?:man)?[- ]?(?:alpha|beta|gamma|α|β)|Lyman[- ](?:limit|forest|break)"
    r"|H[- ]?(?:alpha|beta|gamma|α|β)\b|\[(?:O|N|S|Ne) ?I{1,3}\]"
    r"|\b(?:C|Si|Fe|Mg|Al|O|N|Ca|Zn|S|Ni|Cr|Mn|Na|He) ?(?:IV|VI|V|I{1,3})\*?(?![A-Za-z])")
_LINE_WORDS = re.compile(r"absorption|emission|\blines?\b|features?|trough|doublet", re.I)
_URL = re.compile(r"https?://[^\s)>\]\"'<]+")
_IMAGE = re.compile(r"\.(?:png|gif|jpe?g|webp)$", re.I)
_DATA = re.compile(r"\.(?:fits?|fits\.gz|fit\.gz|txt|dat|csv|json|pdf|ecsv|amon)$|/(?:healpix|json)$", re.I)


def _num(s: str) -> float:
    return float(re.sub(r"\s", "", s).replace("−", "-").replace("–", "-"))


def _quantity(m: re.Match) -> dict[str, Any]:
    exp = m.group("e1") or m.group("e2")
    k = 10.0 ** _num(exp) if exp else 1.0
    lo = hi = None
    if m.group("sym"):
        lo = hi = _num(m.group("sym"))
    elif m.group("a1"):
        for sign, a in ((m.group("s1"), m.group("a1")), (m.group("s2"), m.group("a2"))):
            if sign == "+":
                hi = _num(a)
            else:
                lo = _num(a)

    def scale(x):
        return None if x is None else float(f"{x * k:.6g}")
    return {"value": scale(_num(m.group("v"))), "err_minus": scale(lo), "err_plus": scale(hi),
            "approx": bool(m.group("approx"))}


def _sentences(body: str) -> list[str]:
    out = []
    for part in re.split(r"\n\s*\n|(?<=[.;])\s+(?=[A-Z(])", body):
        s = " ".join(part.split())
        if s and s.count("|") < 2:                     # a table: columns cannot be told apart
            out.append(s)
    return out


def _mask_ranges(s: str) -> str:
    """Energy ranges replaced by blanks: their numbers are not the value."""
    return _RANGE.sub(lambda m: " " * len(m.group(0)), s)


def _band(s: str, at: int) -> Optional[str]:
    """Energy range written in the sentence, the nearest to ``at``."""
    hits = [(abs(m.start() - at), " ".join(m.group(0).split())) for m in _RANGE.finditer(s)]
    return min(hits)[1] if hits else None


def _soft(s: str) -> bool:
    """The sentence gives an energy range in keV that ends below 150 keV."""
    for m in _RANGE.finditer(s):
        top = re.search(r"(" + _N + r")\s*(keV|MeV|GeV)\s*$", m.group(0), re.I)
        if top and top.group(2).lower() == "kev" and float(top.group(1)) <= 150 and ">" not in m.group(0):
            return True
    return False


def _instrument(sentence: str, subject: str) -> Optional[str]:
    for text in (sentence, subject):
        m = _INSTRUMENTS.search(text)
        if m:
            return m.group(0)
    return None


def facts_in(subject: str, body: str, bands: Optional[list[str]] = None) -> list[dict[str, Any]]:
    """Facts of a circular (see the module). ``bands`` are its labels
    (:func:`raven_kn.circulars.classify`): a flux in erg/cm^2/s that is not a
    peak flux is an X-ray flux in an X-ray circular, or when its sentence
    gives a range below 150 keV."""
    bands = bands or []
    out: list[dict[str, Any]] = []
    seen: set[tuple] = set()

    def add(category, kind, s, q=None, unit="", band=None, what=None, span=None):
        key = (kind, None if q is None else q["value"], what)
        if key in seen:
            return
        seen.add(key)
        out.append({"category": category, "kind": kind,
                    **(q or {"value": None, "err_minus": None, "err_plus": None, "approx": False}),
                    "unit": unit, "band": band, "instrument": _instrument(s, subject),
                    # "an upper limit on ... of -3.1", "< 2.6e-10": said just before the value
                    "limit": span is not None and bool(re.search(
                        r"(?:upper|lower) limit|<\s*$", s[max(0, span[0] - 60):span[1]], re.I)),
                    "what": what, "text": s[:400]})

    for s in _sentences(body):
        masked = _mask_ranges(s)
        used: list[tuple[int, int]] = []
        for kind, (category, anchor, unit_rx, unit) in _LABELS.items():
            if kind == "decay_index" and not re.search(r"decay|fading|fade|\bt\s?\^", s, re.I):
                continue
            if kind == "photon_index" and re.search(r"decay|fading", s, re.I) and "index" not in s.lower():
                continue
            for a in re.finditer(anchor, masked, re.I if kind not in ("epeak",) else 0):
                if any(lo <= a.start() < hi for lo, hi in used):
                    continue
                if unit_rx is None:                    # the value follows the words at once
                    g = re.compile(_GLUE).match(masked, a.end())
                    m = _QTY.match(masked, g.end())
                    if m is None:
                        continue
                    tail = masked[m.end():m.end() + 14]
                    u = unit
                    if kind in ("epeak", "epeak_rest"):
                        um = re.match(r"\s*(keV|MeV)", tail)
                        u = um.group(1) if um else unit
                    elif kind == "t90" and re.match(r"\s*(?:ms|msec)\b", tail):
                        u = "ms"
                else:                                  # the first value followed by the unit
                    m = u = None
                    for c in _QTY.finditer(masked, a.end(), min(len(masked), a.end() + 160)):
                        um = re.compile(r"\s*" + unit_rx, re.I).match(masked, c.end())
                        if um:
                            m = c
                            u = unit or ("ph/cm^2/s" if um.group(0).strip().lower().startswith("ph")
                                         else "erg/cm^2/s")
                            break
                    if m is None:
                        continue
                if any(lo <= m.start() < hi for lo, hi in used):
                    continue
                used.append((a.start(), m.end()))
                add(category, kind, s, _quantity(m), u, _band(s, a.start()), span=(a.start(), m.start()))
                break
        # a flux that is not the peak flux of the burst: X-ray, in an X-ray report
        if "xray" in bands or _soft(s):
            for c in _QTY.finditer(masked):
                um = _FLUX_UNIT.match(masked, c.end() + len(masked[c.end():]) - len(masked[c.end():].lstrip()))
                if um and not any(lo <= c.start() < hi for lo, hi in used) \
                        and re.search(r"flux", masked[:c.start()], re.I):
                    add("afterglow", "xray_flux", s, _quantity(c), "erg/cm^2/s", _band(s, c.start()))
                    break
        for c in _QTY.finditer(masked):
            um = _RADIO_UNIT.match(masked, c.end() + len(masked[c.end():]) - len(masked[c.end():].lstrip()))
            if um:
                f = _FREQ.search(s)
                add("afterglow", "radio_flux", s, _quantity(c), um.group(1).replace("μ", "u").replace("µ", "u"),
                    " ".join(f.group(0).split()) if f else None)
        m = _CLASS.search(s)
        if m:
            add("prompt", "burst_class", s, what=" ".join(x for x in (m.group(1), m.group(2)) if x).lower())
        if _LINE_WORDS.search(s):
            species = []
            for x in _SPECIES.finditer(s):
                name = " ".join(x.group(0).split())
                if name not in species:
                    species.append(name)
            if species:
                add("spectroscopy", "lines", s, what=", ".join(species))
        if re.search(r"\bhost(?: galaxy)?\b", s, re.I) and not re.search(r"hosted|on behalf", s, re.I):
            add("spectroscopy", "host", s, what=s[:400])
    return out


def links_in(body: str) -> list[dict[str, Any]]:
    """Web addresses of a text: ``url``, ``kind`` (``image``, ``data`` or
    ``page``), ``about`` (``light curve``, ``sky map``, ``diagram`` or None,
    from the words before it and its name), ``text``, the words that
    introduce it, and ``general``: True for an address that is not about the
    event (the page of a telescope or of a collaboration: no number of five
    digits in it). References to other circulars are left out."""
    out, seen = [], set()
    for m in _URL.finditer(body):
        url = m.group(0).rstrip(".,;:")
        if url in seen or re.search(r"gcn\.nasa\.gov/circulars", url):
            continue
        seen.add(url)
        before = body[max(0, m.start() - 300):m.start()]
        intro = re.split(r"(?<=[.;])\s+(?=[A-Z(])|\n\s*\n|https?://\S+", before)[-1]
        intro = " ".join(intro.split()).rstrip(":( ").strip()
        hint = f"{intro} {url.rsplit('/', 1)[-1]}".lower()
        about = ("diagram" if re.search(r"diagram|relation|amati|yonetoku", hint) else
                 "sky map" if re.search(r"sky ?map|locali[sz]ation|healpix|cover map", hint) else
                 "light curve" if re.search(r"light ?curve|_lc_|lc\.", hint) else None)
        kind = "image" if _IMAGE.search(url) else "data" if _DATA.search(url) else "page"
        # an address about the event carries its number or its date (GRB261113A, bn261111477,
        # id=3532560); the home page of a telescope, of a collaboration or of a program does not
        path = re.sub(r"^https?://[^/]+", "", url)
        out.append({"url": url, "kind": kind, "about": about, "text": intro[:200],
                    "general": kind != "image" and not re.search(r"\d{5}", path)})
    return out
