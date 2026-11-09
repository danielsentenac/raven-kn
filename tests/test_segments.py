"""Excluded time for PyGRB on-source windows: veto definer + science gaps (offline flags)."""
import pytest

pytest.importorskip("gwpy")
lsctables = pytest.importorskip("igwn_ligolw.lsctables")

from raven_kn.targeted.segments import excluded_segments, veto_definer_segments  # noqa: E402


def write_vdf(path):
    from igwn_ligolw import ligolw, utils
    doc = ligolw.Document()
    doc.appendChild(ligolw.LIGO_LW())
    tab = lsctables.New(lsctables.VetoDefTable)
    doc.childNodes[0].appendChild(tab)
    rows = [("H1", "DCH-GLITCH", 1, 2, 0, 0, -1, 1),     # CAT2, padded 1 s each side (start_pad < 0)
            ("L1", "DCH-BAD", 1, 1, 0, 0, 0, 0),         # CAT1
            ("H1", "DCH-MILD", 1, 3, 0, 0, 0, 0),        # CAT3: not used
            ("V1", "DCH-OLD", 1, 2, 0, 1000, 0, 0)]      # valid before the window only
    for ifo, name, ver, cat, s, e, sp, ep in rows:
        r = tab.RowType()
        r.process_id, r.ifo, r.name, r.version, r.category = 0, ifo, name, ver, cat
        r.start_time, r.end_time, r.start_pad, r.end_pad, r.comment = s, e, sp, ep, ""
        tab.append(r)
    utils.write_filename(doc, str(path))
    return str(path)


ACTIVE = {"H1:DCH-GLITCH:1": [(1000100, 1000102)], "L1:DCH-BAD:1": [(1000500, 1000510)],
          "H1:DCH-MILD:1": [(1000700, 1000800)], "V1:DCH-OLD:1": [(1000300, 1000400)],
          "H1:DMT-ANALYSIS_READY:1": [(1000000, 1000900)],
          "L1:DMT-ANALYSIS_READY:1": [(1000000, 1001000)]}


def test_veto_definer_categories_and_padding(tmp_path):
    vdf = write_vdf(tmp_path / "vdf.xml")
    v = veto_definer_segments(vdf, ["H1", "L1", "V1"], 1000000, 1001000, active=ACTIVE)
    assert v == [(1000099.0, 1000103.0), (1000500.0, 1000510.0)]
    assert veto_definer_segments(vdf, ["L1"], 1000000, 1001000, active=ACTIVE) == [(1000500.0, 1000510.0)]


def test_excluded_with_science_and_subsegments(tmp_path):
    from raven_kn.targeted.pygrb import veto_free_subsegments
    vdf = write_vdf(tmp_path / "vdf.xml")
    ex = excluded_segments(["H1", "L1"], 1000000, 1001000, veto_definer=vdf,
                           science_flags={"H1": "H1:DMT-ANALYSIS_READY:1",
                                          "L1": "L1:DMT-ANALYSIS_READY:1"}, active=ACTIVE)
    assert ex == [(1000099.0, 1000103.0), (1000500.0, 1000510.0), (1000900.0, 1001000.0)]
    segs = veto_free_subsegments(1000000, 1001000, 256, ex)
    assert all(b <= a2 - 1 or a >= b2 + 1 for a, b in segs for a2, b2 in ex)
