"""Time to exclude from targeted-search on-source windows.

PyGRB aborts if its on-source window contains vetoed time or lies outside
science data. Before placing sub-segments
(:func:`raven_kn.targeted.pygrb.veto_free_subsegments`), the window is
cleaned of, for every detector analysed:

* time outside the science flag (e.g. ``H1:DMT-ANALYSIS_READY:1``, or
  ``H1:DATA`` for open data);
* the veto-definer flags of the requested categories (PyGRB uses CAT1 and
  CAT2: ``segments-vetoes = +CBC_CAT2`` in LVK's configs), with each flag's
  start/end padding and validity interval.

Queries use ``gwpy`` and DQSegDB (``https://segments.ligo.org``, which needs
IGWN credentials; on LDG clusters they are picked up automatically).
"""

from __future__ import annotations

from typing import Iterable, Optional, Sequence

DQSEGDB = "https://segments.ligo.org"


def _as_list(segs) -> list[tuple[float, float]]:
    return [(float(a), float(b)) for a, b in segs]


def veto_definer_segments(veto_definer: str, ifos: Iterable[str], start: float, end: float,
                          categories: Sequence[int] = (1, 2), url: str = DQSEGDB,
                          active: Optional[dict] = None) -> list[tuple[float, float]]:
    """Union of the active, padded veto segments of ``categories`` for
    ``ifos`` in [start, end]. ``active``: {flag name: segment list} to use
    instead of querying (tests, offline use)."""
    from gwpy.segments import DataQualityDict, Segment, SegmentList
    vdf = DataQualityDict.from_veto_definer_file(veto_definer, start=start, end=end,
                                                 ifo=None)
    ifos = set(ifos)
    flags = DataQualityDict({k: f for k, f in vdf.items()
                             if f.category in categories and k.split(":")[0] in ifos})
    if not flags:
        return []
    if active is None:
        flags.populate(source=url, segments=[Segment(start, end)], on_error="raise")
    else:
        for k in list(flags):
            f = flags[k]
            f.known = SegmentList([Segment(start, end)])
            f.active = SegmentList(Segment(*s) for s in active.get(k, []))
            flags[k] = f.pad()          # the flag's own (start, end) padding
    union = SegmentList()
    for f in flags.values():
        union |= f.active & SegmentList([Segment(start, end)])
    return _as_list(union.coalesce())


def science_gaps(science_flags: dict[str, str], start: float, end: float,
                 url: str = DQSEGDB, active: Optional[dict] = None) -> list[tuple[float, float]]:
    """Times in [start, end] where at least one detector is not in science
    (complement of each science flag, unioned). ``science_flags``: {ifo:
    flag}, e.g. {"H1": "H1:DMT-ANALYSIS_READY:1"}."""
    from gwpy.segments import DataQualityFlag, Segment, SegmentList
    span = SegmentList([Segment(start, end)])
    gaps = SegmentList()
    for ifo, name in science_flags.items():
        if active is not None:
            sci = SegmentList(Segment(*s) for s in active.get(name, []))
        else:
            sci = DataQualityFlag.query(name, start, end, url=url).active
        gaps |= span - sci
    return _as_list(gaps.coalesce())


def excluded_segments(ifos: Iterable[str], start: float, end: float,
                      veto_definer: Optional[str] = None, categories: Sequence[int] = (1, 2),
                      science_flags: Optional[dict[str, str]] = None, url: str = DQSEGDB,
                      active: Optional[dict] = None) -> list[tuple[float, float]]:
    """Science gaps and vetoes, merged: pass to ``write_pygrb_overrides(...,
    vetoes=...)`` or ``raven-kn pygrb --vetoes``."""
    from gwpy.segments import Segment, SegmentList
    out = SegmentList()
    if science_flags:
        out |= SegmentList(Segment(*s) for s in science_gaps(science_flags, start, end, url, active))
    if veto_definer:
        out |= SegmentList(Segment(*s) for s in veto_definer_segments(
            veto_definer, ifos, start, end, categories, url, active))
    return _as_list(out.coalesce())


def write_segments(segs: Sequence[tuple[float, float]], path: str) -> None:
    with open(path, "w") as f:
        f.write("# excluded GPS intervals (start end), from raven_kn.targeted.segments\n")
        for a, b in segs:
            f.write(f"{a:.3f} {b:.3f}\n")
