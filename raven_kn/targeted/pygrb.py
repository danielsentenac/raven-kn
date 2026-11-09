"""Inputs for the production targeted search, PyGRB (``pycbc_make_offline_grb_workflow``).

PyGRB is LVK's coherent targeted search for GRBs: on-source windows of
seconds to minutes around a trigger time, a sky point (or grid), coherent
SNR with chi-square and null-stream vetoes, off-source background. A
kilonova t0 window lasts hours to days, so it is split into sub-segments
and one PyGRB analysis is configured per sub-segment (as in Jamsin et al.).

This module does not produce a complete PyGRB configuration: executables,
data discovery, segment/veto definitions, injections and site settings come
from the LVK PyGRB base configuration. It writes, per sub-segment, an
*override* file with the event-specific options, to be layered on the base:

    pycbc_make_offline_grb_workflow --config-files base.ini <override>.ini ...

Option names checked against PyCBC 2.8 (``pycbc.workflow.segment``,
``matched_filter``, ``tmpltbank``) and LVK's O4 PyGRB test configuration
(pycbc-config ``O4/pygrb/dev/o2_tests/GRB170817A.ini``): ``[workflow]
trigger-name, trigger-time, ra, dec ("<value> rad"), sky-error, start-time,
end-time``; ``[workflow-exttrig_segments] on-before,
on-after, min-duration, max-duration``; ``[workflow-tmpltbank]
tmpltbank-method = PREGENERATED_BANK, tmpltbank-pregenerated-bank``.
``pad-data`` and ``quanta`` stay in the base file: the off-source length must
be compatible with them (PyGRB checks it).
"""

from __future__ import annotations

import configparser
import math
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from .search import TargetedSearchConfig, _mass_equal, build_bank


def write_bank_hdf(mchirps: Sequence[float], path: str | Path, f_lower: float) -> Path:
    """PyCBC HDF template bank (equal-mass, non-spinning) from chirp masses."""
    import h5py
    path = Path(path)
    m = np.array([_mass_equal(mc) for mc in mchirps])
    with h5py.File(path, "w") as f:
        f["mass1"] = m
        f["mass2"] = m
        f["spin1z"] = np.zeros_like(m)
        f["spin2z"] = np.zeros_like(m)
        f["f_lower"] = np.full_like(m, f_lower)
        # precomputed hashes, as in PyCBC banks (PyCBC otherwise computes them
        # with np.fromstring, which fails with NumPy >= 2)
        import hashlib
        f["template_hash"] = np.array([np.frombuffer(hashlib.blake2b(
            np.array((a, a, 0.0, 0.0)).tobytes(), digest_size=8).digest(), dtype=np.int64)[0]
            for a in m])
        f.attrs["parameters"] = ["mass1", "mass2", "spin1z", "spin2z", "f_lower", "template_hash"]
    return path


def veto_free_subsegments(gps_start: float, gps_end: float, max_len: float,
                          vetoes: Sequence[tuple[float, float]] = (), min_len: float = 16.0,
                          guard_s: float = 1.0) -> list[tuple[float, float]]:
    """Split [gps_start, gps_end] into on-source sub-segments of at most
    ``max_len`` seconds that contain no vetoed time (PyGRB aborts otherwise).

    ``vetoes``: (start, end) GPS intervals, all detectors merged (e.g. CAT1/
    CAT2 from the veto definer via DQSegDB). Each is widened by ``guard_s``;
    the veto-free stretches are cut into equal pieces of <= ``max_len``;
    pieces shorter than ``min_len`` are dropped (their time is lost to the
    search and should be reported).
    """
    vs = sorted((max(a - guard_s, gps_start), min(b + guard_s, gps_end))
                for a, b in vetoes if b > gps_start and a < gps_end)
    free, t = [], gps_start
    for a, b in vs:
        if a > t:
            free.append((t, a))
        t = max(t, b)
    if t < gps_end:
        free.append((t, gps_end))
    out = []
    for a, b in free:
        if b - a < min_len:
            continue
        n = int(math.ceil((b - a) / max_len))
        e = np.linspace(a, b, n + 1)
        out += [(float(x), float(y)) for x, y in zip(e[:-1], e[1:])]
    return out


def write_pygrb_overrides(cfg: TargetedSearchConfig, out_dir: str | Path,
                          bank: Optional[Sequence[float]] = None,
                          offsource_factor: float = 6.0,
                          base_config: str | Sequence[str] = "pygrb_base.ini",
                          event_template: Optional[str | Path] = None,
                          bank_path: Optional[str] = None,
                          extra_overrides: Sequence[str] = (),
                          vetoes: Sequence[tuple[float, float]] = ()) -> list[Path]:
    """Per sub-segment of the t0 window: an event ini file and the
    ``pycbc_make_offline_grb_workflow`` command, in ``run_pygrb.sh``; one
    shared HDF bank.

    PyCBC's config parser rejects an option defined in two files, so the
    event file *replaces* the per-GRB file of the base configuration
    (``event_template``, e.g. LVK's ``GRB170817A.ini``: its ``[workflow]``
    event options are replaced, its other sections, e.g. injection
    settings, kept, and references to the template's own path pointed to
    the new file), and the options that live in the base analysis file are
    passed as ``--config-overrides``:
    ``workflow-exttrig_segments:on-before / on-after / min-duration /
    max-duration`` and ``workflow-tmpltbank:tmpltbank-pregenerated-bank``.

    For a sub-segment of length L centred on T: trigger-time = T,
    on-before = on-after = L/2, off-source between 4 x L/2 and
    ``offsource_factor`` x L. ``bank_path``: path of the bank as seen where
    the workflow runs (default: local path). ``vetoes``: vetoed GPS
    intervals; sub-segment edges are placed so that no on-source contains
    them (:func:`veto_free_subsegments`). Integer GPS edges are used, so
    each on-source is shrunk inward to whole seconds.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    bank = list(bank) if bank is not None else build_bank(
        cfg.mchirp_min, cfg.mchirp_max, cfg.f_low, cfg.min_match, cfg.sample_rate, cfg.max_templates)
    bank_file = write_bank_hdf(bank, out / f"{cfg.label}_bank.hdf", cfg.f_low)
    bank_ref = bank_path or str(bank_file.resolve())
    base = [base_config] if isinstance(base_config, str) else list(base_config)
    if vetoes:
        pieces = veto_free_subsegments(cfg.gps_start, cfg.gps_end, cfg.subsegment_s, vetoes)
    else:
        edges = np.arange(cfg.gps_start, cfg.gps_end, cfg.subsegment_s).tolist() + [cfg.gps_end]
        pieces = [(a, b) for a, b in zip(edges[:-1], edges[1:]) if b - a >= 16]
    files, cmds = [], []
    for k, (a, b) in enumerate(pieces):
        if vetoes:      # stay strictly inside the veto-free piece
            a, b = math.ceil(a), math.floor(b)
            half = int((b - a) // 2)
            t = a + half
        else:
            t = int(round(0.5 * (a + b)))
            half = int(math.ceil(0.5 * (b - a)))
        span = int(math.ceil(offsource_factor * 2 * half))
        name = f"{cfg.label}_seg{k:03d}"
        f = out / f"{name}.ini"
        cp = configparser.ConfigParser(interpolation=None)
        cp.optionxform = str
        if event_template is not None:
            cp.read(event_template)
            tpl = str(event_template)
            for sec in cp.sections():
                for opt, val in cp[sec].items():
                    if tpl in val or Path(tpl).name in val:
                        cp[sec][opt] = val.replace(tpl, str(f.name)).replace(Path(tpl).name, f.name)
        if not cp.has_section("workflow"):
            cp.add_section("workflow")
        # "<value> <unit>" with a space, as in LVK's PyGRB configs: a negative
        # "<value><unit>" is misparsed by argparse (pycbc.types.angle_as_radians)
        cp["workflow"].update({"trigger-name": name, "trigger-time": str(t),
                               "ra": f"{math.radians(cfg.ra_deg):.12f} rad",
                               "dec": f"{math.radians(cfg.dec_deg):.12f} rad",
                               "sky-error": "0.0 rad",
                               "start-time": str(t - span // 2 - 64),
                               "end-time": str(t + span // 2 + 64)})
        with open(f, "w") as fh:
            fh.write(f"; RAVEN-KN event file for PyGRB, sub-segment {k} of {cfg.label}\n"
                     f"; on-source [{a:.1f}, {b:.1f}] GPS ({b - a:.0f} s)"
                     + (f"; from template {Path(event_template).name}" if event_template else "")
                     + "\n")
            cp.write(fh)
        files.append(f)
        overrides = [f"workflow-exttrig_segments:on-before:{half}",
                     f"workflow-exttrig_segments:on-after:{half}",
                     f"workflow-exttrig_segments:min-duration:{4 * half}",
                     f"workflow-exttrig_segments:max-duration:{span}",
                     f"workflow-tmpltbank:tmpltbank-pregenerated-bank:{bank_ref}",
                     *extra_overrides]
        cmds.append(f"mkdir -p {name}\n"   # PyGRB saves its config into --output
                    "pycbc_make_offline_grb_workflow \\\n  --config-files " + " ".join(base)
                    + f" {f.name} \\\n  --workflow-name {name} --output {name} \\\n"
                    "  --config-overrides " + " ".join(overrides))
    script = out / "run_pygrb.sh"
    script.write_text("#!/bin/sh\n# Generated by RAVEN-KN: one PyGRB workflow per sub-segment.\n"
                      f"# Base configuration files: {' '.join(base)} (LVK, not provided here).\n"
                      "set -e\n" + "\n\n".join(cmds) + "\n")
    script.chmod(0o755)
    return files
