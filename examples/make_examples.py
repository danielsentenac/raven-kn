"""Write example candidate files.

gw170817_like.json: real position (AT2017gfo), host redshift (NGC 4993,
z = 0.009783) and VLBI viewing angle (~20 deg, Mooley+2018); the
photometry is SYNTHETIC (Rubin-like, Jamsin et al. scenario: g/r
non-detections 3 h before the merger, first detection at +10.9 h).
far_candidate.json: synthetic candidate at z = 0.15 (~ 700 Mpc): beyond a
useful GW search.
"""
import json
from pathlib import Path

from astropy.time import Time

HERE = Path(__file__).parent
T0_GPS = 1187008882.43                       # GW170817 merger
t0_mjd = Time(T0_GPS, format="gps").utc.mjd
h = 1 / 24

phot = [
    {"mjd": t0_mjd - 3 * h, "band": "g", "mag": 24.5, "upper_limit": True, "survey": "LSST"},
    {"mjd": t0_mjd - 3 * h + 0.02, "band": "r", "mag": 24.3, "upper_limit": True, "survey": "LSST"},
    {"mjd": t0_mjd + 10.9 * h, "band": "i", "mag": 17.1, "mag_err": 0.03, "survey": "LSST"},
    {"mjd": t0_mjd + 11.4 * h, "band": "g", "mag": 17.3, "mag_err": 0.03, "survey": "LSST"},
    {"mjd": t0_mjd + 1.4, "band": "r", "mag": 17.6, "mag_err": 0.03, "survey": "LSST"},
    {"mjd": t0_mjd + 2.4, "band": "r", "mag": 18.3, "mag_err": 0.04, "survey": "LSST"},
]

gw170817 = {
    "name": "AT2017gfo-like", "ra": 197.450374, "dec": -23.381495, "pos_err_arcsec": 0.1,
    "source": "example", "kn_score": 0.9, "photometry": phot,
    "redshift": {"z": 0.009783, "z_err": 1e-5, "kind": "spectroscopic", "source": "host NGC 4993"},
    "viewing_angle": {"mean_deg": 20.0, "std_deg": 5.0, "method": "VLBI superluminal motion"},
    "notes": "position/redshift/viewing angle real; photometry synthetic",
}
gw170817_grb = dict(gw170817, name="AT2017gfo-like+GRB",
                    grb={"gps": T0_GPS + 1.74, "instrument": "Fermi-GBM"})

far = {
    "name": "far-candidate", "ra": 150.1, "dec": 2.2, "source": "example", "kn_score": 0.7,
    "photometry": [
        {"mjd": t0_mjd - 2.0, "band": "g", "mag": 24.0, "upper_limit": True, "survey": "LSST"},
        {"mjd": t0_mjd + 0.6, "band": "g", "mag": 23.0, "mag_err": 0.1, "survey": "LSST"},
        {"mjd": t0_mjd + 0.62, "band": "i", "mag": 22.8, "mag_err": 0.1, "survey": "LSST"},
    ],
    "redshift": {"z": 0.15, "z_err": 5e-4, "kind": "spectroscopic", "source": "host"},
    "notes": "synthetic",
}

for name, d in (("gw170817_like", gw170817), ("gw170817_like_grb", gw170817_grb),
                ("far_candidate", far)):
    (HERE / f"{name}.json").write_text(json.dumps(d, indent=2))
    print("wrote", HERE / f"{name}.json")
