"""NMMA wrapper: kilonova light-curve fit -> merger-time (t0) posterior.

Turns a :class:`KNCandidate` into NMMA (>= 1.0, ``lightcurve-analysis``)
inputs, runs it (or prints the command for a cluster), and reads the
posterior back into the candidate:

* ``kn.t0_samples_gps``  <- trigger_time + ``timeshift``
* ``kn.lc_fit_samples``  <- ``luminosity_distance``, ``inclination`` (from
  NMMA's ``inclination_EM``, a viewing angle in [0, pi/2])

Interface checked against the NMMA 1.0.1 source:

* data file: ``<isot> <filter> <mag> <mag_err>`` per line, non-detections
  with ``mag_err = inf`` (``nmma.em.io.strict_read_csv``); the model returns
  no flux before the merger, so a non-detection before t0 is automatically
  consistent and one after t0 penalises it: that is how a pre-merger
  upper limit narrows t0 (Jamsin et al., arXiv:2610.28775);
* data earlier than ``--trigger-time`` are cut (``--data-tmin`` = 0), so the
  trigger time is set *before* the earliest non-detection kept;
* merger time = trigger time + ``timeshift`` (days), detector frame
  (``nmma.em.model``: observed = source x (1+z) + timeshift);
* the first detection must be later than ``(1+z_max) * em_tmin +
  timeshift_max`` (``check_model_time_consistency``), which bounds the
  timeshift prior;
* results: bilby ``<label>_posterior_samples.dat`` or ``<label>_result.json``.

Priors: by default the EM distance/inclination posterior (spectroscopic
redshift, VLBI, ...) is passed to NMMA as priors ("informed"), which helps
t0 a little (Jamsin et al.: ~5 % with a known distance). The fitted distance
and inclination are then NOT reused as an independent constraint (that
would count the redshift twice): ``kn.lc_fit_informed`` is set.
"""

from __future__ import annotations

import json
import math
import shlex
import shutil
import subprocess
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
from astropy.time import Time

from ..connection import KNCandidate
from .distance_inclination import combine, constraints_from_candidate

DAY = 86400.0

# survey -> {band: NMMA filter}
FILTER_MAP = {
    "ZTF": {"g": "ztfg", "r": "ztfr", "i": "ztfi"},
    "LSST": {b: f"lsst{b}" for b in "ugrizy"},
    "RUBIN": {b: f"lsst{b}" for b in "ugrizy"},
    "PS1": {b: f"ps1::{b}" for b in "grizy"},
    "PAN-STARRS": {b: f"ps1::{b}" for b in "grizy"},
    "SDSS": {"u": "sdssu"},
}

# Model-parameter priors (everything except distance, inclination, timeshift),
# copied from the NMMA 1.0.1 prior files (GPL-3, nmma/priors/<model>.prior).
# For other models (e.g. Bu2026, used by Jamsin et al.) pass ``model_prior``.
MODEL_PRIORS = {
    "Bu2019lm": [
        "KNphi = Uniform(name='KNphi', minimum=15., maximum=75.)",
        "log10_mej_dyn = Uniform(name='log10_mej_dyn', minimum=-3., maximum=-1.0)",
        "log10_mej_wind = Uniform(name='log10_mej_wind', minimum=-3., maximum=-0.5)",
    ],
    "Bu2023Ye": [
        "log10_mej_dyn = Uniform(name='log10_mej_dyn', minimum=-3., maximum=-1.7)",
        "vej_dyn = Uniform(name='vej_dyn', minimum=0.12, maximum=0.25)",
        "Yedyn = Uniform(name='Yedyn', minimum=0.15, maximum=0.30)",
        "log10_mej_wind = Uniform(name='log10_mej_wind', minimum=-2., maximum=-0.89)",
        "vej_wind = Uniform(name='vej_wind', minimum=0.03, maximum=0.15)",
        "Yewind = Uniform(name='Yewind', minimum=0.20, maximum=0.40)",
    ],
}
_OURS = ("luminosity_distance", "inclination_EM", "timeshift")


def nmma_filter(survey: str, band: str) -> str:
    m = FILTER_MAP.get(survey.upper(), {})
    return m.get(band, band)        # unknown: assume band is already an NMMA name


@dataclass
class NMMASetup:
    workdir: str
    label: str
    model: str
    trigger_mjd: float
    timeshift_bounds: tuple[float, float]
    informed: bool
    data_file: str
    prior_file: str
    command: list[str]
    n_detections: int
    n_limits: int
    notes: list[str] = field(default_factory=list)

    @property
    def command_line(self) -> str:
        return " ".join(shlex.quote(c) for c in self.command)

    def save(self) -> Path:
        p = Path(self.workdir) / "raven_kn_nmma_setup.json"
        p.write_text(json.dumps(asdict(self), indent=2))
        return p

    @classmethod
    def load(cls, workdir: str | Path) -> "NMMASetup":
        d = json.loads((Path(workdir) / "raven_kn_nmma_setup.json").read_text())
        d["timeshift_bounds"] = tuple(d["timeshift_bounds"])
        return cls(**d)


def _model_prior_lines(model: str, model_prior: Optional[str | Path]) -> list[str]:
    if model_prior is not None:
        lines = Path(model_prior).read_text().splitlines()
    elif model in MODEL_PRIORS:
        lines = MODEL_PRIORS[model]
    else:
        raise ValueError(f"No built-in prior for {model}: pass model_prior=<NMMA prior file>")
    keep = []
    for ln in lines:
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        if s.split("=")[0].strip() in _OURS:
            continue
        keep.append(s)
    return keep


def build_prior(kn: KNCandidate, model: str, timeshift_bounds: tuple[float, float],
                informed: bool = True, model_prior: Optional[str | Path] = None,
                d_max: float = 1500.0) -> tuple[list[str], list[str]]:
    """Prior file lines and notes."""
    lines = _model_prior_lines(model, model_prior)
    notes = []
    a, b = timeshift_bounds
    lines.append(f"timeshift = Uniform(minimum={a:.5f}, maximum={b:.5f}, name='timeshift')")
    cons = [c for c in constraints_from_candidate(kn) if not c.name.startswith("lc_fit")]
    if informed and cons:
        post = combine(cons, d_max=d_max)
        m, s = float(np.mean(post.d_mpc)), float(np.std(post.d_mpc))
        lines.append(f"luminosity_distance = TruncatedGaussian(mu={m:.3f}, sigma={max(s, 0.5):.3f}, "
                     f"minimum={max(m - 6 * s, 0.1):.3f}, maximum={m + 6 * s:.3f}, "
                     f"name='luminosity_distance')")
        notes.append(f"distance prior from EM data: {m:.1f} +- {s:.1f} Mpc")
    else:
        lines.append(f"luminosity_distance = Uniform(minimum=1.0, maximum={d_max:.1f}, "
                     f"name='luminosity_distance')")
    if informed and kn.viewing_angle is not None:
        mu, sg = math.radians(kn.viewing_angle.mean_deg), math.radians(kn.viewing_angle.std_deg)
        lines.append(f"inclination_EM = TruncatedGaussian(mu={mu:.5f}, sigma={sg:.5f}, minimum=0., "
                     f"maximum={math.pi / 2:.6f}, name='inclination_EM')")
        notes.append(f"viewing-angle prior: {kn.viewing_angle.mean_deg} +- {kn.viewing_angle.std_deg} deg")
    else:
        lines.append(f"inclination_EM = Sine(minimum=0., maximum={math.pi / 2:.6f}, name='inclination_EM')")
    return lines, notes


def write_lightcurve(kn: KNCandidate, path: str | Path, t_start_mjd: float,
                     t_end_mjd: float = np.inf) -> tuple[int, int]:
    """NMMA observations file; returns (n_detections, n_limits) written."""
    rows, nd, nl = [], 0, 0
    for p in sorted(kn.photometry, key=lambda p: p.mjd):
        if p.mjd < t_start_mjd or p.mjd > t_end_mjd:
            continue
        isot = Time(p.mjd, format="mjd", scale="utc").isot
        filt = nmma_filter(p.survey, p.band)
        if p.upper_limit:
            rows.append(f"{isot} {filt} {p.mag:.3f} inf")
            nl += 1
        else:
            rows.append(f"{isot} {filt} {p.mag:.3f} {p.mag_err if p.mag_err else 0.1:.3f}")
            nd += 1
    Path(path).write_text("\n".join(rows) + "\n")
    return nd, nl


def setup(kn: KNCandidate, workdir: str | Path, model: str = "Bu2019lm",
          informed: bool = True, model_prior: Optional[str | Path] = None,
          max_age_d: float = 7.0, lookback_d: float = 3.0,
          em_tmin: float = 0.1, em_tmax: float = 14.0, z_max: float = 0.12,
          interpolation_type: str = "tensorflow", sampler: str = "pymultinest",
          svd_path: Optional[str] = None, em_error_budget: float = 1.0,
          extra_args: Optional[list[str]] = None, executable: str = "lightcurve-analysis") -> NMMASetup:
    """Write the NMMA data and prior files and build the command.

    The merger is searched in [first detection - ``max_age_d``, first
    detection - (1 + z_max) * em_tmin]. Non-detections up to ``lookback_d``
    before that are kept (they constrain t0 through the model). The trigger
    time is put at the earliest kept point.
    """
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    fd = kn.first_detection
    if fd is None:
        raise ValueError(f"{kn.name}: NMMA needs at least one detection")
    earliest_merger = fd.mjd - max_age_d
    t_start = earliest_merger - lookback_d
    kept_ul = [p.mjd for p in kn.upper_limits if t_start <= p.mjd < fd.mjd]
    trigger = min([earliest_merger] + kept_ul) - 1e-3
    ts_lo = earliest_merger - trigger
    ts_hi = (fd.mjd - trigger) - (1 + z_max) * em_tmin
    if ts_hi <= ts_lo:
        raise ValueError("timeshift range empty: lower em_tmin or raise max_age_d")
    label = kn.name.replace(" ", "_").replace("+", "_")
    data_file, prior_file = wd / f"{label}.dat", wd / f"{label}.prior"
    # NMMA requires every detection inside the model span for the earliest
    # allowed merger; later points carry no t0 information: drop them.
    t_end = earliest_merger + em_tmax - 0.01
    n_late = sum(p.mjd > t_end for p in kn.photometry)
    nd, nl = write_lightcurve(kn, data_file, trigger, t_end)
    lines, notes = build_prior(kn, model, (ts_lo, ts_hi), informed, model_prior)
    prior_file.write_text("\n".join(lines) + "\n")
    if n_late:
        notes.append(f"dropped {n_late} points later than earliest merger + em_tmax ({em_tmax} d)")
    cmd = [executable, "--em-model", model, "--light-curve-data", str(data_file),
           "--trigger-time", f"{trigger:.6f}", "--prior-file", str(prior_file),
           "--outdir", str(wd / "nmma_out"), "--label", label,
           "--em-tmin", str(em_tmin), "--em-tmax", str(em_tmax),
           "--interpolation-type", interpolation_type, "--sampler", sampler,
           "--em-error-budget", str(em_error_budget), "--plot"]
    if svd_path:
        cmd += ["--svd-path", str(svd_path)]
    cmd += list(extra_args or [])
    s = NMMASetup(str(wd), label, model, trigger, (ts_lo, ts_hi), informed and bool(
        constraints_from_candidate(kn)), str(data_file), str(prior_file), cmd, nd, nl, notes)
    s.save()
    return s


def run(s: NMMASetup, timeout: Optional[float] = None) -> subprocess.CompletedProcess:
    if shutil.which(s.command[0]) is None:
        raise FileNotFoundError(f"{s.command[0]} not found: install nmma (pip install nmma) or "
                                f"run on a cluster: {s.command_line}")
    log = Path(s.workdir) / "nmma.log"
    with open(log, "w") as f:
        return subprocess.run(s.command, stdout=f, stderr=subprocess.STDOUT, timeout=timeout,
                              check=True)


def read_posterior(outdir: str | Path, label: str) -> dict[str, np.ndarray]:
    """bilby posterior: ``<label>_posterior_samples.dat`` or ``<label>_result.json``."""
    outdir = Path(outdir)
    dat = outdir / f"{label}_posterior_samples.dat"
    if dat.exists():
        arr = np.genfromtxt(dat, names=True)
        return {k: np.asarray(arr[k], float) for k in arr.dtype.names}
    js = outdir / f"{label}_result.json"
    if js.exists():
        post = json.loads(js.read_text())["posterior"]
        post = post.get("content", post)
        return {k: np.asarray(v, float) for k, v in post.items()
                if isinstance(v, list) and v and isinstance(v[0], (int, float))}
    raise FileNotFoundError(f"no NMMA posterior for {label} in {outdir}")


def apply_posterior(kn: KNCandidate, s: NMMASetup, post: dict[str, np.ndarray]) -> KNCandidate:
    """Put the NMMA posterior into the candidate (t0 in GPS, distance, inclination)."""
    if "timeshift" not in post:
        raise KeyError("posterior has no 'timeshift'")
    t0_mjd = s.trigger_mjd + post["timeshift"]
    kn.t0_samples_gps = [float(x) for x in Time(t0_mjd, format="mjd", scale="utc").gps]
    lc = {}
    if "luminosity_distance" in post:
        lc["luminosity_distance"] = post["luminosity_distance"].tolist()
    if "inclination_EM" in post:
        lc["inclination"] = post["inclination_EM"].tolist()
    kn.lc_fit_samples = lc or None
    kn.lc_fit_informed = s.informed
    kn.notes = (kn.notes + f"; NMMA {s.model} fit").lstrip("; ")
    return kn


def fit(kn: KNCandidate, workdir: str | Path, dry_run: bool = False, **kw) -> tuple[KNCandidate, NMMASetup]:
    """Setup, run NMMA, read back. With ``dry_run`` only the files and the
    command are produced (run it elsewhere, then :func:`import_result`)."""
    s = setup(kn, workdir, **kw)
    if not dry_run:
        run(s)
        apply_posterior(kn, s, read_posterior(Path(s.workdir) / "nmma_out", s.label))
    return kn, s


def import_result(kn: KNCandidate, workdir: str | Path) -> KNCandidate:
    s = NMMASetup.load(workdir)
    return apply_posterior(kn, s, read_posterior(Path(s.workdir) / "nmma_out", s.label))
