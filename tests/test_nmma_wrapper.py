import os
import stat
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.time import Time

from raven_kn.connection import gps_to_mjd
from raven_kn.treatment import lightcurve_fit as lcf
from raven_kn.treatment.distance_inclination import constraints_from_candidate
from raven_kn.treatment.t0 import estimate_window, window_from_bracket

from test_pipeline import T0, make_kn

FAKE = Path(__file__).with_name("fake_lightcurve_analysis.py")


def _fake_nmma_bin(tmp_path):
    """conda-like bin/ with lightcurve-analysis and python."""
    b = tmp_path / "envbin"
    b.mkdir()
    exe = b / "lightcurve-analysis"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    (b / "python").symlink_to(sys.executable)
    return str(exe)


def test_setup_files(tmp_path):
    kn = make_kn()
    s = lcf.setup(kn, tmp_path / "w", executable=_fake_nmma_bin(tmp_path))
    data = Path(s.data_file).read_text().split("\n")
    assert any(l.endswith(" inf") and " ps1::g " in l for l in data)      # pre-merger LSST g UL kept (PS1 g model)
    assert s.n_limits == 1 and s.n_detections == 2
    # trigger before the UL, so NMMA's data-tmin=0 cut keeps it
    ul_mjd = kn.upper_limits[0].mjd
    assert s.trigger_mjd < ul_mjd
    # merger range: from first detection - 7 d to just before the first detection
    lo, hi = s.timeshift_bounds
    fd = kn.first_detection.mjd
    assert s.trigger_mjd + lo == pytest.approx(fd - 7, abs=1e-6)
    assert s.trigger_mjd + hi < fd
    assert s.trigger_mjd + lo < gps_to_mjd(T0) < s.trigger_mjd + hi
    prior = Path(s.prior_file).read_text()
    assert "TruncatedGaussian(mu=4" in prior            # D_L ~ 43 Mpc from the redshift
    assert "inclination_EM = TruncatedGaussian" in prior
    assert s.informed and s.command[1].endswith("nmma_launcher.py")
    assert lcf.NMMASetup.load(tmp_path / "w").trigger_mjd == s.trigger_mjd


def test_limit_only_filter_dropped_without_launcher(tmp_path):
    s = lcf.setup(make_kn(), tmp_path, executable="no-such-nmma-binary")
    data = Path(s.data_file).read_text()
    assert " ps1::g " not in data and s.command[0] == "no-such-nmma-binary"
    assert any("non-detection-only" in n for n in s.notes)


def test_launcher_patch_handles_limit_only_filter(tmp_path):
    """Fake nmma package reproducing the upstream crash; the launcher must avoid it."""
    import subprocess
    pkg = tmp_path / "fake" / "nmma" / "em"
    pkg.mkdir(parents=True)
    (pkg.parent / "__init__.py").write_text("")
    (pkg / "__init__.py").write_text("")
    (pkg / "utils.py").write_text(
        "import numpy as np\n"
        "def check_model_time_consistency(d, m, p, injection=None):\n"
        "    t, mag, e, tt = d\n"
        "    for k in t:\n"
        "        det = np.isfinite(mag[k]) & np.isfinite(e[k])\n"
        "        t[k][det].min()\n"
        "    return d\n")
    (pkg / "analysis.py").write_text(
        "import numpy as np, sys\n"
        "from . import utils\n"
        "def main():\n"
        "    inf = np.inf\n"
        "    d = ({'r': np.array([1., 2.]), 'g': np.array([0.5])},\n"
        "         {'r': np.array([20., 21.]), 'g': np.array([24.])},\n"
        "         {'r': np.array([.1, .1]), 'g': np.array([inf])}, 0.)\n"
        "    out = utils.check_model_time_consistency(d, None, {})\n"
        "    assert 'g' in out[0]\n"
        "    print('OK', sys.argv[1:])\n"
        "    return 0\n")
    launcher = Path(lcf.__file__).with_name("nmma_launcher.py")
    env = dict(os.environ, PYTHONPATH=str(tmp_path / "fake"))
    r = subprocess.run([sys.executable, str(launcher), "--x", "1"], env=env,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "OK ['--x', '1']" in r.stdout


def test_prior_file_parses_with_bilby(tmp_path):
    bilby = pytest.importorskip("bilby")
    s = lcf.setup(make_kn(), tmp_path, model="Bu2023Ye")
    pd = bilby.core.prior.PriorDict(filename=s.prior_file)
    assert {"timeshift", "luminosity_distance", "inclination_EM", "Yewind"} <= set(pd)
    smp = pd.sample(1000)
    assert np.all(smp["inclination_EM"] <= np.pi / 2)


def test_uninformed_priors(tmp_path):
    s = lcf.setup(make_kn(), tmp_path, informed=True)
    assert s.informed
    s2 = lcf.setup(make_kn(z=None, theta=None), tmp_path / "u", informed=True)
    assert not s2.informed                    # nothing to inform with
    assert "Uniform(minimum=1.0" in Path(s2.prior_file).read_text()


def test_unknown_model_needs_prior(tmp_path):
    with pytest.raises(ValueError):
        lcf.setup(make_kn(), tmp_path, model="Bu2026")
    pf = tmp_path / "Bu2026.prior"
    pf.write_text("log10_mej = Uniform(minimum=-3, maximum=-1, name='log10_mej')\n"
                  "timeshift = Uniform(minimum=-2, maximum=1, name='timeshift')\n")
    s = lcf.setup(make_kn(), tmp_path / "w", model="Bu2026", model_prior=pf)
    txt = Path(s.prior_file).read_text()
    assert "log10_mej" in txt and txt.count("timeshift =") == 1


def test_run_with_fake_nmma_narrows_window(tmp_path, monkeypatch):
    exe = tmp_path / "lightcurve-analysis"
    exe.write_text(f"#!/bin/sh\nexec {sys.executable} {FAKE} \"$@\"\n")
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{tmp_path}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_T0_MJD", str(gps_to_mjd(T0)))
    kn = make_kn()
    bracket = window_from_bracket(kn, d_mpc_upper=50)
    kn, s = lcf.fit(kn, tmp_path / "run")
    assert len(kn.t0_samples_gps) == 4000
    assert np.median(kn.t0_samples_gps) == pytest.approx(T0, abs=0.05 * 86400)
    w = estimate_window(kn, d_mpc_upper=50)
    assert w.gps_start < T0 < w.gps_end
    assert w.duration < bracket.duration
    # informed fit: its D/iota must not be reused as a second constraint
    assert kn.lc_fit_informed
    assert not any(c.name.startswith("lc_fit") for c in constraints_from_candidate(kn))
    # re-import from disk gives the same samples
    kn2 = lcf.import_result(make_kn(), tmp_path / "run")
    assert kn2.t0_samples_gps == kn.t0_samples_gps


def test_missing_executable(tmp_path):
    s = lcf.setup(make_kn(), tmp_path, executable="no-such-nmma-binary")
    with pytest.raises(FileNotFoundError, match="install nmma"):
        lcf.run(s)


def test_filter_map():
    assert lcf.nmma_filter("ZTF", "r") == "ztfr"
    assert lcf.nmma_filter("LSST", "y") == "ps1::y"     # NMMA has no LSST filter models
    assert lcf.nmma_filter("RUBIN", "u") == "sdssu"
    assert lcf.nmma_filter("PS1", "z") == "ps1::z"
    assert lcf.nmma_filter("whatever", "besselli") == "besselli"


def test_timeshift_bound_respects_distance_prior(tmp_path):
    """NMMA check: first detection >= (1 + z_max) * em_tmin + timeshift_max,
    with z_max from the distance prior (uninformed: up to 1500 Mpc, z ~ 0.28)."""
    from astropy import units as u
    from astropy.cosmology import Planck18, z_at_value
    kn = make_kn(z=None, theta=None)
    s = lcf.setup(kn, tmp_path, executable="no-such-nmma-binary")
    z_hi = float(z_at_value(Planck18.luminosity_distance, 1500 * u.Mpc))
    first = kn.first_detection.mjd - s.trigger_mjd
    assert (1 + z_hi) * 0.1 + s.timeshift_bounds[1] < first


def test_priors_clipped_to_model_training_range(tmp_path):
    import joblib
    core = {"ps1__g": {"param_mins": np.array([-3., -2., 0., 0.]),
                       "param_maxs": np.array([-1.699, -0.886, 90., 90.])}}
    joblib.dump(core, tmp_path / "Bu2019lm.joblib")
    r = lcf.model_training_ranges("Bu2019lm", tmp_path)
    assert r["log10_mej_wind"] == (-2.0, -0.886)
    assert r["inclination_EM"][1] == pytest.approx(np.pi / 2)
    wide = tmp_path / "wide.prior"          # NMMA 1.0.1 Bu2019lm.prior ranges
    wide.write_text("log10_mej_dyn = Uniform(minimum=-3., maximum=-1.0, name='log10_mej_dyn')\n"
                    "log10_mej_wind = Uniform(minimum=-3., maximum=-0.5, name='log10_mej_wind')\n")
    s = lcf.setup(make_kn(), tmp_path / "w", model_prior=wide, svd_path=str(tmp_path),
                  executable="no-such-nmma-binary")
    txt = Path(s.prior_file).read_text()
    assert "maximum=-1.699" in txt and "minimum=-2" in txt and "maximum=-0.886" in txt
    assert sum("clipped to training range" in n for n in s.notes) == 2


def test_import_result_after_moving_run(tmp_path, monkeypatch):
    exe = tmp_path / "lightcurve-analysis"
    exe.write_text(f"#!/bin/sh\nexec {sys.executable} {FAKE} \"$@\"\n")
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("PATH", f"{tmp_path}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_T0_MJD", str(gps_to_mjd(T0)))
    lcf.fit(make_kn(), tmp_path / "run", use_launcher=False)
    (tmp_path / "run").rename(tmp_path / "moved")
    kn = lcf.import_result(make_kn(), tmp_path / "moved")
    assert len(kn.t0_samples_gps) == 4000


def test_fit_days_cuts_late_points(tmp_path):
    kn = make_kn()                       # detections at +10.9 h and +1.4 d
    s = lcf.setup(kn, tmp_path / "a", fit_days=0.5, executable="x")
    assert s.n_detections == 1 and any("first detection + 0.5 d" in n for n in s.notes)
    s = lcf.setup(kn, tmp_path / "b", fit_days=None, executable="x")
    assert s.n_detections == 2
