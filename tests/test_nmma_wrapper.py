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


def test_setup_files(tmp_path):
    kn = make_kn()
    s = lcf.setup(kn, tmp_path)
    data = Path(s.data_file).read_text().split("\n")
    assert any(l.endswith(" inf") and " lsstg " in l for l in data)       # pre-merger UL kept
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
    assert s.informed and s.command[0] == "lightcurve-analysis"
    assert lcf.NMMASetup.load(tmp_path).trigger_mjd == s.trigger_mjd


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
    assert lcf.nmma_filter("LSST", "y") == "lssty"
    assert lcf.nmma_filter("PS1", "z") == "ps1::z"
    assert lcf.nmma_filter("whatever", "besselli") == "besselli"
