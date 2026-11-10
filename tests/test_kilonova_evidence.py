"""Kilonova confirmation before a GW search is triggered."""

from pathlib import Path

import numpy as np
from astropy.time import Time

from raven_kn.connection import GRBAssociation, KNCandidate, PhotometryPoint, Redshift
from raven_kn.search import assess
from raven_kn.treatment.kilonova_evidence import band_wavelength, kilonova_evidence

ROOT = Path(__file__).resolve().parents[1]
GFO = ROOT / "examples" / "validation" / "AT2017gfo_reduced.dat"
GRB_GPS = 1187008884.47
RA, DEC = 197.450371, -23.381486


def gfo_points(max_days):
    """Real AT2017gfo photometry (NMMA's reduced file) up to ``max_days``
    after the first point."""
    rows = [l.split() for l in GFO.read_text().splitlines() if l.strip()]
    t = Time([r[0] for r in rows], format="isot", scale="utc").mjd
    pts = [PhotometryPoint(float(ti), r[1], float(r[2]), float(r[3]), survey="lit")
           for ti, r in zip(t, rows) if np.isfinite(float(r[3]))]
    t1 = min(p.mjd for p in pts)
    return [p for p in pts if p.mjd - t1 <= max_days]


def kn_with(points, **kw):
    return KNCandidate("SSS17a", RA, DEC, photometry=points,
                       redshift=Redshift(0.009783, 1e-5, "spectroscopic", "NGC 4993"), **kw)


def test_band_wavelengths():
    assert band_wavelength("ps1::g") < band_wavelength("ps1::i")
    assert band_wavelength("ztfr") == band_wavelength("r")
    assert band_wavelength("sdssu") < band_wavelength("lsstg")
    assert band_wavelength("Ks") > band_wavelength("J")
    assert band_wavelength("o") > band_wavelength("c")


def test_at2017gfo_light_curve_confirms_within_two_days():
    early = kilonova_evidence(kn_with(gfo_points(0.3)), 40.0)
    assert not early.confirmed                       # one night: no evolution measurable yet
    ev = kilonova_evidence(kn_with(gfo_points(2.0)), 40.0)
    assert ev.confirmed and ev.by == ["light curve"]
    assert -17 < ev.peak_abs_mag < -15
    assert max(ev.max_decline_mag_per_day, ev.max_reddening_mag_per_day) >= 0.3


def test_single_detection_needs_a_spectrum():
    one = [PhotometryPoint(Time("2017-08-17T23:33:17").mjd, "i", 16.0, 0.1, survey="Swope")]
    ev = kilonova_evidence(kn_with(one), 40.0)
    assert not ev.confirmed and "1 detection epoch" in ev.reasons[-1]
    ev = kilonova_evidence(kn_with(one, classification="Kilonova"), 40.0)
    assert ev.confirmed and ev.by == ["spectrum"]


def test_supernova_like_not_confirmed():
    t0 = 60000.0
    sn = [PhotometryPoint(t0 + d, b, 17.0 + 0.03 * d + (0.1 if b == "r" else 0.0), 0.05)
          for d in (0, 1, 2, 3) for b in ("g", "r")]
    ev = kilonova_evidence(kn_with(sn), 100.0)       # M = -18, slow
    assert not ev.confirmed
    assert "outside the kilonova range" in ev.reasons[-1]


def test_grb_plus_one_point_does_not_trigger_until_confirmed(tmp_path):
    one = [PhotometryPoint(Time("2017-08-17T23:33:17").mjd, "i", 16.0, 0.1, survey="Swope")]
    kn = kn_with(one, grb=GRBAssociation(GRB_GPS, "Fermi-GBM"))
    d = assess(kn, tmp_path)
    assert d["detectability"]["worth_search"] and not d["trigger_targeted_search"]
    assert not (tmp_path / "SSS17a_targeted_search_request.json").exists()
    kn.photometry = gfo_points(2.0) + one
    d = assess(kn, tmp_path)
    assert d["trigger_targeted_search"] and d["t0_window"]["duration_h"] < 0.01
