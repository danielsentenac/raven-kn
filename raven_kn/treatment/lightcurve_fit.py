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
  timeshift_max`` (z_max from the distance prior's upper bound) (``check_model_time_consistency``), which bounds the
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
    # NMMA's distributed models have no LSST filters (none of the 9 in its model
    # repository, 2026-11): Rubin bands use the nearest PS1 band (g,r,i,z,y)
    # and SDSS u, as in the Rubin injection calibration.
    "LSST": {**{b: f"ps1::{b}" for b in "grizy"}, "u": "sdssu"},
    "RUBIN": {**{b: f"ps1::{b}" for b in "grizy"}, "u": "sdssu"},
    "PS1": {b: f"ps1::{b}" for b in "grizy"},
    "PAN-STARRS": {b: f"ps1::{b}" for b in "grizy"},
    "SDSS": {"u": "sdssu"},
    "ATLAS": {"c": "atlasc", "o": "atlaso"},
}

# Model-parameter priors (everything except distance, inclination, timeshift),
# copied from the NMMA 1.0.1 prior files (GPL-3, nmma/priors/<model>.prior).
# For other models (e.g. Bu2026, used by Jamsin et al.) pass ``model_prior``.
MODEL_PRIORS = {
    # NMMA's Bu2019lm.prior (dyn up to -1.0, wind -3..-0.5) exceeds the
    # training grid of the distributed Bu2019lm SVD model (dyn -3..-1.699,
    # wind -2..-0.886): outside it the GP extrapolates. Training range used.
    "Bu2019lm": [
        "KNphi = Uniform(name='KNphi', minimum=15., maximum=75.)",
        "log10_mej_dyn = Uniform(name='log10_mej_dyn', minimum=-3., maximum=-1.699)",
        "log10_mej_wind = Uniform(name='log10_mej_wind', minimum=-2., maximum=-0.886)",
    ],
    # Ka2017 (Kasen et al. 2017, one component, no viewing angle):
    # NMMA 1.0.2 priors/Ka2017.prior; clipped to the training range if known.
    "Ka2017": [
        "log10_mej = Uniform(name='log10_mej', minimum=-3., maximum=-1.)",
        "log10_vej = Uniform(name='log10_vej', minimum=-1.52, maximum=-0.53)",
        "log10_Xlan = Uniform(name='log10_Xlan', minimum=-9, maximum=-1)",
    ],
    # LANL SuperNu grid (NMMA "LANLTP2", tensorflow interpolation only; arXiv:2105.11543):
    # toroidal dynamical + peanut-shaped wind ejecta. Ranges = the grid's training
    # range (param_mins/maxs of LANLTP2.joblib); NMMA's LANL2022 prior starts the
    # velocities at 0, outside the grid.
    "LANLTP2": [
        "log10_mej_dyn = Uniform(name='log10_mej_dyn', minimum=-3., maximum=-1.)",
        "vej_dyn = Uniform(name='vej_dyn', minimum=0.05, maximum=0.3)",
        "log10_mej_wind = Uniform(name='log10_mej_wind', minimum=-3., maximum=-1.)",
        "vej_wind = Uniform(name='vej_wind', minimum=0.05, maximum=0.3)",
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

# Order of the training parameters stored in NMMA SVD model files
# (``param_mins`` / ``param_maxs``), from nmma.em.model (v1.0.2).
MODEL_PARAM_NAMES = {
    "Bu2019lm": ["log10_mej_dyn", "log10_mej_wind", "KNphi", "KNtheta"],
    "Bu2019nsbh": ["log10_mej_dyn", "log10_mej_wind", "KNtheta"],
    "Bu2019lm_sparse": ["log10_mej_dyn", "log10_mej_wind"],
    "Ka2017": ["log10_mej", "log10_vej", "log10_Xlan"],
    "LANLTP2": ["log10_mej_dyn", "vej_dyn", "log10_mej_wind", "vej_wind", "KNtheta"],
}


def model_training_ranges(model: str, svd_path: Optional[str | Path],
                          param_names: Optional[list[str]] = None) -> dict[str, tuple[float, float]]:
    """Training range of each model parameter, read from ``<svd_path>/<model>.joblib``.

    ``KNtheta`` (deg) is returned as ``inclination_EM`` (rad). Empty dict if
    the file or the parameter names are unavailable.
    """
    names = param_names or MODEL_PARAM_NAMES.get(model)
    if svd_path is None or names is None:
        return {}
    f = Path(svd_path).expanduser() / f"{model}.joblib"
    if not f.exists():
        return {}
    import joblib
    core = joblib.load(f)
    entry = next(iter(core.values()))
    lo, hi = np.asarray(entry["param_mins"], float), np.asarray(entry["param_maxs"], float)
    if len(lo) != len(names):
        return {}
    out = {}
    for n, a, b in zip(names, lo, hi):
        if n == "KNtheta":
            out["inclination_EM"] = (math.radians(a), math.radians(b))
        else:
            out[n] = (float(a), float(b))
    return out


def clip_prior_lines(lines: list[str], ranges: dict[str, tuple[float, float]]) -> tuple[list[str], list[str]]:
    """Narrow ``minimum``/``maximum`` of each prior to the training range."""
    import re
    out, notes = [], []
    for ln in lines:
        name = ln.split("=")[0].strip()
        if name in ranges:
            a, b = ranges[name]
            mn = re.search(r"minimum\s*=\s*([-+0-9.eE]+)", ln)
            mx = re.search(r"maximum\s*=\s*([-+0-9.eE]+)", ln)
            if mn and mx:
                lo, hi = float(mn.group(1)), float(mx.group(1))
                nlo, nhi = max(lo, a), min(hi, b)
                if nlo >= nhi:
                    raise ValueError(f"prior on {name} [{lo}, {hi}] outside training range [{a}, {b}]")
                if (nlo, nhi) != (lo, hi):
                    ln = ln[:mn.start(1)] + f"{nlo:.6g}" + ln[mn.end(1):]
                    mx = re.search(r"maximum\s*=\s*([-+0-9.eE]+)", ln)
                    ln = ln[:mx.start(1)] + f"{nhi:.6g}" + ln[mx.end(1):]
                    notes.append(f"{name} prior [{lo:g}, {hi:g}] clipped to training range "
                                 f"[{nlo:.4g}, {nhi:.4g}]")
        out.append(ln)
    return out, notes


# NMMA interpolation type each distributed model exists in (sklearn_gp ``.joblib``
# per filter, or a tensorflow network per filter in ``<model>_tf/``)
MODEL_INTERPOLATION = {"Bu2019lm": "sklearn_gp", "Ka2017": "sklearn_gp", "LANLTP2": "tensorflow"}


def nmma_filter(survey: str, band: str) -> str:
    m = FILTER_MAP.get(survey.upper(), {})
    if band in m:
        return m[band]
    if len(band) <= 1:              # a bare band letter NMMA has no model file for
        import warnings
        warnings.warn(f"survey {survey!r}, band {band!r}: no NMMA filter mapping; add it to "
                      "FILTER_MAP or give the NMMA filter name as the band")
    return band                     # assume band is already an NMMA name


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


def distance_prior(kn: KNCandidate, informed: bool = True,
                   d_max: float = 1500.0) -> tuple[str, Optional[str], float]:
    """(prior line, note, upper bound in Mpc) for ``luminosity_distance``."""
    cons = [c for c in constraints_from_candidate(kn) if not c.name.startswith("lc_fit")]
    if informed and cons:
        post = combine(cons, d_max=d_max)
        m, s = float(np.mean(post.d_mpc)), float(np.std(post.d_mpc))
        hi = m + 6 * s
        return (f"luminosity_distance = TruncatedGaussian(mu={m:.3f}, sigma={max(s, 0.5):.3f}, "
                f"minimum={max(m - 6 * s, 0.1):.3f}, maximum={hi:.3f}, "
                f"name='luminosity_distance')",
                f"distance prior from EM data: {m:.1f} +- {s:.1f} Mpc", hi)
    return (f"luminosity_distance = Uniform(minimum=1.0, maximum={d_max:.1f}, "
            f"name='luminosity_distance')", None, d_max)


def build_prior(kn: KNCandidate, model: str, timeshift_bounds: tuple[float, float],
                informed: bool = True, model_prior: Optional[str | Path] = None,
                d_max: float = 1500.0) -> tuple[list[str], list[str]]:
    """Prior file lines and notes."""
    lines = _model_prior_lines(model, model_prior)
    notes = []
    a, b = timeshift_bounds
    lines.append(f"timeshift = Uniform(minimum={a:.5f}, maximum={b:.5f}, name='timeshift')")
    line, note, _ = distance_prior(kn, informed, d_max)
    lines.append(line)
    if note:
        notes.append(note)
    if informed and kn.viewing_angle is not None:
        mu, sg = math.radians(kn.viewing_angle.mean_deg), math.radians(kn.viewing_angle.std_deg)
        lines.append(f"inclination_EM = TruncatedGaussian(mu={mu:.5f}, sigma={sg:.5f}, minimum=0., "
                     f"maximum={math.pi / 2:.6f}, name='inclination_EM')")
        notes.append(f"viewing-angle prior: {kn.viewing_angle.mean_deg} +- {kn.viewing_angle.std_deg} deg")
    else:
        lines.append(f"inclination_EM = Sine(minimum=0., maximum={math.pi / 2:.6f}, name='inclination_EM')")
    return lines, notes


def write_lightcurve(kn: KNCandidate, path: str | Path, t_start_mjd: float,
                     t_end_mjd: float = np.inf,
                     drop_limit_only_filters: bool = False) -> tuple[int, int, list[str]]:
    """NMMA observations file; returns (n_detections, n_limits, dropped filters).

    ``drop_limit_only_filters``: remove filters with only non-detections,
    needed when NMMA is run without :mod:`nmma_launcher` (NMMA 1.0.x
    ``check_model_time_consistency`` crashes on them).
    """
    pts = [p for p in kn.photometry if t_start_mjd <= p.mjd <= t_end_mjd]
    all_f = {nmma_filter(p.survey, p.band) for p in pts}
    with_det = {nmma_filter(p.survey, p.band) for p in pts if not p.upper_limit}
    keep = with_det if drop_limit_only_filters else all_f
    dropped = sorted(all_f - keep)
    rows, nd, nl = [], 0, 0
    for p in sorted(pts, key=lambda p: p.mjd):
        filt = nmma_filter(p.survey, p.band)
        if filt not in keep:
            continue
        isot = Time(p.mjd, format="mjd", scale="utc").isot
        if p.upper_limit:
            rows.append(f"{isot} {filt} {p.mag:.3f} inf")
            nl += 1
        else:
            rows.append(f"{isot} {filt} {p.mag:.3f} {p.mag_err if p.mag_err else 0.1:.3f}")
            nd += 1
    Path(path).write_text("\n".join(rows) + "\n")
    return nd, nl, dropped


def setup(kn: KNCandidate, workdir: str | Path, model: str = "Bu2019lm",
          informed: bool = True, model_prior: Optional[str | Path] = None,
          max_age_d: float = 7.0, lookback_d: float = 3.0,
          em_tmin: float = 0.1, em_tmax: float = 14.0, d_max: float = 1500.0,
          fit_days: Optional[float] = 3.0,
          interpolation_type: str = "tensorflow", sampler: str = "pymultinest",
          svd_path: Optional[str] = None, em_error_budget: float = 1.0,
          extra_args: Optional[list[str]] = None, executable: str = "lightcurve-analysis",
          use_launcher: bool = True) -> NMMASetup:
    """Write the NMMA data and prior files and build the command.

    The merger is searched in [first detection - ``max_age_d``, first
    detection - (1 + z_max) * em_tmin], z_max being the redshift of the
    largest distance allowed by the distance prior. Non-detections up to ``lookback_d``
    before that are kept (they constrain t0 through the model). The trigger
    time is put at the earliest kept point.

    ``fit_days``: only data up to ``fit_days`` after the first detection are
    fitted (None: all). Later points carry little t0 information, and when
    the model cannot follow the late evolution they bias t0: with Bu2019lm
    on AT2017gfo, all data give a merger ~48 h early, the first ~3 days
    ~16 h early (chi2 profile, README).

    ``use_launcher``: if ``executable`` is an NMMA ``lightcurve-analysis``
    with a ``python`` next to it (conda/venv layout), run NMMA through
    :mod:`raven_kn.treatment.nmma_launcher` with that python, which works
    around NMMA's crash on non-detection-only filters; otherwise such
    filters are dropped (noted).
    """
    wd = Path(workdir)
    prefix = [executable]
    exe_path = shutil.which(executable)
    if use_launcher and exe_path and Path(exe_path).name == "lightcurve-analysis":
        py = Path(exe_path).parent / "python"
        if py.exists():
            prefix = [str(py), str(Path(__file__).with_name("nmma_launcher.py"))]
    wd.mkdir(parents=True, exist_ok=True)
    fd = kn.first_detection
    if fd is None:
        raise ValueError(f"{kn.name}: NMMA needs at least one detection")
    earliest_merger = fd.mjd - max_age_d
    t_start = earliest_merger - lookback_d
    kept_ul = [p.mjd for p in kn.upper_limits if t_start <= p.mjd < fd.mjd]
    trigger = min([earliest_merger] + kept_ul) - 1e-3
    ts_lo = earliest_merger - trigger
    # NMMA requires the model start, stretched by (1 + z_max) and shifted by
    # the largest timeshift, to precede the first detection
    from astropy import units as u
    from astropy.cosmology import Planck18, z_at_value
    d_hi = distance_prior(kn, informed, d_max)[2]
    z_max = 1.02 * float(z_at_value(Planck18.luminosity_distance, d_hi * u.Mpc))
    ts_hi = (fd.mjd - trigger) - (1 + z_max) * em_tmin - 1e-3
    if ts_hi <= ts_lo:
        raise ValueError("timeshift range empty: lower em_tmin or raise max_age_d")
    label = kn.name.replace(" ", "_").replace("+", "_")
    data_file, prior_file = wd / f"{label}.dat", wd / f"{label}.prior"
    # NMMA requires every detection inside the model span for the earliest
    # allowed merger; later points carry no t0 information: drop them.
    t_end = earliest_merger + em_tmax - 0.01
    if fit_days is not None:
        t_end = min(t_end, fd.mjd + fit_days)
    n_late = sum(p.mjd > t_end for p in kn.photometry)
    nd, nl, dropped = write_lightcurve(kn, data_file, trigger, t_end,
                                       drop_limit_only_filters=len(prefix) == 1)
    lines, notes = build_prior(kn, model, (ts_lo, ts_hi), informed, model_prior, d_max)
    ranges = model_training_ranges(model, svd_path)
    if ranges:
        lines, clip_notes = clip_prior_lines(lines, ranges)
        notes += clip_notes
    elif svd_path is None:
        notes.append("model training ranges not checked (no svd_path): priors may extrapolate the model")
    prior_file.write_text("\n".join(lines) + "\n")
    if dropped:
        notes.append(f"dropped non-detection-only filters {dropped}: NMMA run without the "
                     f"launcher crashes on them")
    if n_late:
        notes.append(f"dropped {n_late} points later than "
                     + (f"first detection + {fit_days} d" if fit_days is not None and
                        t_end == fd.mjd + fit_days else f"earliest merger + em_tmax ({em_tmax} d)"))
    cmd = prefix + ["--em-model", model, "--light-curve-data", str(data_file),
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
    # read from the given directory (the run may have been moved since setup)
    return apply_posterior(kn, s, read_posterior(Path(workdir) / "nmma_out", s.label))
