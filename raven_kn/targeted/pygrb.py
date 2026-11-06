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
``matched_filter``, ``tmpltbank``): ``[workflow] trigger-name, trigger-time,
ra, dec, start-time, end-time``; ``[workflow-exttrig_segments] on-before,
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


def write_pygrb_overrides(cfg: TargetedSearchConfig, out_dir: str | Path,
                          bank: Optional[Sequence[float]] = None,
                          offsource_factor: float = 6.0,
                          base_config: str = "pygrb_base.ini") -> list[Path]:
    """One override ini per sub-segment of the t0 window, a shared HDF bank,
    and ``run_pygrb.sh`` with the workflow-generation commands.

    For a sub-segment of length L centred on T: trigger-time = T,
    on-before = on-after = L/2, off-source between ``2 L`` and
    ``offsource_factor`` x L around it (start/end-time bound the data used).
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    bank = list(bank) if bank is not None else build_bank(
        cfg.mchirp_min, cfg.mchirp_max, cfg.f_low, cfg.min_match, cfg.sample_rate, cfg.max_templates)
    bank_file = write_bank_hdf(bank, out / f"{cfg.label}_bank.hdf", cfg.f_low)
    edges = np.arange(cfg.gps_start, cfg.gps_end, cfg.subsegment_s).tolist() + [cfg.gps_end]
    files, cmds = [], []
    for k, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        if b - a < 16:
            continue
        t = int(round(0.5 * (a + b)))
        half = int(math.ceil(0.5 * (b - a)))
        span = int(math.ceil(offsource_factor * 2 * half))
        name = f"{cfg.label}_seg{k:03d}"
        cp = configparser.ConfigParser()
        cp.optionxform = str
        cp["workflow"] = {"trigger-name": name, "trigger-time": str(t),
                          "ra": f"{cfg.ra_deg:.6f}deg", "dec": f"{cfg.dec_deg:.6f}deg",
                          "start-time": str(t - span // 2 - 64), "end-time": str(t + span // 2 + 64)}
        cp["workflow-exttrig_segments"] = {"on-before": str(half), "on-after": str(half),
                                           "min-duration": str(4 * half),
                                           "max-duration": str(span)}
        cp["workflow-tmpltbank"] = {"tmpltbank-method": "PREGENERATED_BANK",
                                    "tmpltbank-pregenerated-bank": str(bank_file.resolve())}
        f = out / f"{name}.ini"
        with open(f, "w") as fh:
            fh.write(f"; RAVEN-KN override for PyGRB, sub-segment {k} of {cfg.label}\n"
                     f"; on-source [{a:.1f}, {b:.1f}] GPS ({b - a:.0f} s); layer on the LVK "
                     f"PyGRB base configuration\n")
            cp.write(fh)
        files.append(f)
        cmds.append(f"pycbc_make_offline_grb_workflow --config-files {base_config} {f.name} "
                    f"--output-dir {name}")
    script = out / "run_pygrb.sh"
    script.write_text("#!/bin/sh\n# Generated by RAVEN-KN: one PyGRB workflow per sub-segment.\n"
                      f"# {base_config}: LVK PyGRB base configuration (not provided here).\n"
                      "set -e\n" + "\n".join(cmds) + "\n")
    script.chmod(0o755)
    return files
