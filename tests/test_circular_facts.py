"""Quantities and links read in GCN Circulars: sentences of real circulars of
November 2026 (the number of each is given)."""
import pytest

from raven_kn.circular_facts import facts_in, links_in
from raven_kn.circulars import handle_circular, reported
from raven_kn.serve import Store

from test_circulars import circ, store  # noqa: F401  (fixture)


def one(kind, text, subject="GRB 261113A: observation", bands=None):
    got = [f for f in facts_in(subject, text, bands) if f["kind"] == kind]
    assert len(got) == 1, got
    return got[0]


@pytest.mark.parametrize("text, value, minus, plus", [
    ("The SVOM/ECLAIRs light curve showed a multiple narrow peak structure with a T90 duration "
     "of about 7.447 (-0.084 +0.121) s.", 7.447, 0.084, 0.121),                                   # 46244
    ("this burst mainly consists of multiple pulses with a duration (T90) of 41 +3/-1 s.", 41, 1, 3),   # 46245
    ("a single broad peak structure with a T90 duration of about 46 -9 / +6 s.", 46, 9, 6),       # 46267
    ("The burst duration is T90 = 7.8 -0.1 / +0.1 in the 4-120 keV energy band.", 7.8, 0.1, 0.1),  # 46270
    ("multi-pulses with a T90-value of 10.5 +1.5/-2.5 s in the 15-5000 keV band.", 10.5, 2.5, 1.5),  # 46311
    ("multiple pulses with a T90 duration of 7.2 +/-0.2 s in the 15-5000 keV band.", 7.2, 0.2, 0.2),  # 46287
    ("a single emission episode with a duration (T90) of about 17.7 s (50-300 keV).", 17.7, None, None),  # 46320
])
def test_t90_as_written(text, value, minus, plus):
    f = one("t90", text)
    assert (f["value"], f["err_minus"], f["err_plus"], f["unit"]) == (value, minus, plus, "s")
    assert f["category"] == "prompt" and f["text"] == text


def test_energy_ranges_are_the_band_not_the_value():
    f = one("fluence", "The event fluence (10-1000 keV) in this time interval is "
                       "(2.05 +0.15/-0.20) E-05 erg/cm^2.", "GRB 261111A: GECAM-B observation of a burst")
    assert f["value"] == pytest.approx(2.05e-5) and f["err_minus"] == pytest.approx(2.0e-6)
    assert f["err_plus"] == pytest.approx(1.5e-6) and f["band"] == "10-1000 keV"
    assert f["unit"] == "erg/cm^2" and f["instrument"] == "GECAM-B"
    t = one("t90", "The burst duration is T90 = 7.8 -0.1 / +0.1 in the 4-120 keV energy band.")
    assert t["band"] == "4-120 keV"


def test_konus_wind_sentences():                                                                 # 46268
    body = ("As observed by Konus-Wind, the burst had the total fluence of (3.51 ± 0.61)x10^-6 erg/cm^2 "
            "and a 16-ms peak energy flux, measured from T0, of (1.36 ± 0.24)x10^-6 erg/cm^2/s (both in "
            "the 20 keV - 10 MeV energy range).\n\nThe spectrum of the burst (measured from T0 to "
            "T0+8.448 s) is best fitted in the 20 keV - 15 MeV range by a power law with exponential "
            "cutoff (CPL) model dN/dE ~ (E^alpha)*exp(-E*(2+alpha)/Ep) with alpha = -0.81 (-0.30, + 0.33) "
            "and Ep = 106(-9,+10) keV, PGstat = 49/60 dof.\n\nFitting this spectrum by a Band GRB "
            "function yields the same values of alpha and Ep, and an upper limit on the high energy "
            "photon index beta of -3.1\n\nwe estimate the burst isotropic energy release E_iso to "
            "(9.49 ± 1.65)x10^52 erg, the isotropic peak luminosity L_iso to (1.68 ± 0.29)x10^53 erg/s.")
    got = {(f["kind"], f["value"]): f for f in facts_in("Konus-Wind detection of GRB 261113A", body)}
    assert got[("fluence", 3.51e-6)]["err_plus"] == pytest.approx(6.1e-7)
    assert got[("peak_flux", 1.36e-6)]["unit"] == "erg/cm^2/s"
    assert (got[("epeak", 106.0)]["err_minus"], got[("epeak", 106.0)]["err_plus"]) == (9.0, 10.0)
    assert got[("photon_index", -0.81)]["err_plus"] == 0.33 and not got[("photon_index", -0.81)]["limit"]
    assert got[("photon_index", -3.1)]["limit"]
    assert got[("eiso", 9.49e52)]["unit"] == "erg" and got[("liso", 1.68e53)]["unit"] == "erg/s"
    assert all(f["instrument"] == "Konus-Wind" for f in got.values())


def test_photon_flux_and_peak_energy_without_unit():                                              # 46270
    body = ("The time-averaged spectrum from T0-0.68 s to T0+7.16 s in the 5-120 keV energy range is best "
            "fitted by a cutoff power law model with a power-law index alpha = -0.61 ± 0.05 and peak "
            "energy of 108 ± 10. With this model, the 4-120 keV fluence is (2.81 ± 0.09)e-6 erg/cm^2 and "
            "the 4-120 keV photon flux is 6.89 ± 0.15 ph/cm^2/s.")
    got = {f["kind"]: f for f in facts_in("GRB 261113A: SVOM/ECLAIRs refined analysis", body)}
    assert got["photon_index"]["value"] == -0.61 and got["epeak"]["value"] == 108.0
    assert got["fluence"]["value"] == pytest.approx(2.81e-6)
    assert (got["peak_flux"]["value"], got["peak_flux"]["unit"]) == (6.89, "ph/cm^2/s")


def test_afterglow_and_spectroscopy():
    d = one("decay_index", "It shows a decay with a slope index of -0.8 over our 1.5-hour observation "
                           "span.", "GRB 261113A: SVOM/VT optical decay slope")                 # 46254
    assert (d["category"], d["value"], d["instrument"]) == ("afterglow", -0.8, "SVOM/VT")
    r = one("radio_flux", "we detect a radio counterpart at 5 GHz with a flux density of ~0.2 mJy at a "
                          "position consistent with the optical counterpart position.",
            "GRB 261122B: ATCA Detection of Radio Counterpart")                                   # 46351
    assert (r["value"], r["unit"], r["band"], r["approx"], r["instrument"]) == (0.2, "mJy", "5 GHz", True, "ATCA")
    x = one("xray_flux", "The estimated flux of the flare is around 1.4 x 10^(-10) erg/s/cm^2 in "
                         "0.5-4.0 keV, corresponding to an X-ray luminosity of around 5.6 x 10^(31) erg/s.",
            "EP261114a is likely a flaring star")                                                 # 46291
    assert x["value"] == pytest.approx(1.4e-10) and x["band"] == "0.5-4.0 keV"
    n = one("nh", "an absorbed power law with a Galactic hydrogen column density of 1.20 (-1.03/+1.14)×10^21 "
                  "cm^-2 and a photon index of 2.01 (-0.35/+0.39).", "EP261115a: refined analysis", ["xray"])
    assert n["value"] == pytest.approx(1.2e21) and n["err_plus"] == pytest.approx(1.14e21)        # 46316
    ln = one("lines", "Furthermore we detect several narrow metal absorption lines consistent with C II, "
                      "Si II, C IV, and Fe II, which allow us to accurately measure the redshift z = 3.580.")
    assert ln["what"] == "C II, Si II, C IV, Fe II" and ln["category"] == "spectroscopy"          # 46261
    c = one("burst_class", "The Fermi GBM team reports the detection of a likely LONG GRB")       # 46224
    assert c["what"] == "likely long" and c["value"] is None


def test_nothing_is_read_from_a_table_or_from_plain_text():
    assert facts_in("GRB 261113A: optical observations",
                    "We observed the field of GRB 261113A with the 1 m telescope.\n\n"
                    "T-T0 | filter | mag\n0.5 h | r | 19.2\n") == []


def test_links():                                                                                # 46224, 46287
    body = ("The skymap can be found here:\nhttps://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/"
            "2026/bn261111477/quicklook/glg_skymap_all_bn261111477.png\n\nThe HEALPix FITS file can be "
            "found here:\nhttps://heasarc.gsfc.nasa.gov/x/glg_healpix_all_bn261111477.fit\n\n"
            "The SVOM/GRM light curve can be found here: https://www.bursthub.cn/gcn/2026/11/svgrb261113A.png.\n"
            "The localization of GRB 261113A in the 'Amati' relation diagram is shown at:\n"
            "https://www.bursthub.cn/gcn/2026/11/svgrb261113A_amati.png\n"
            "(Kan et al., GCN [46284](https://gcn.nasa.gov/circulars/46284)); details at "
            "https://grb.mpe.mpg.de/grb/GRB261111477/")
    got = links_in(body)
    assert [(x["kind"], x["about"]) for x in got] == [
        ("image", "sky map"), ("data", "sky map"), ("image", "light curve"), ("image", "diagram"),
        ("page", None)]
    assert got[2]["url"].endswith("svgrb261113A.png") and got[0]["text"] == "The skymap can be found here"
    assert [x["general"] for x in got] == [False] * 5
    # 46283: the plot of the burst is about the event, the page of the telescope is not
    hx = links_in("The HXMT/HE light curve can be found here:\nhttps://www.bursthub.cn/gcn/2026/hxmtgrb261113B.png\n\n"
                  "More information about it could be found at: http://hxmten.ihep.ac.cn/ and "
                  "https://fermi.gsfc.nasa.gov/ssc/data/access/gbm/, https://ror.org/05qajvd42. Real time updated "
                  "cover map available here: https://master.sai.msu.ru/site/master2/observ.php?id=3532560")
    assert [(x["kind"], x["about"], x["general"]) for x in hx] == [
        ("image", "light curve", False), ("page", None, True), ("page", None, True), ("page", None, True),
        ("page", "sky map", False)]


def test_reported_by_alert_and_service(store):  # noqa: F811
    for n in (46244, 46261, 46284):
        handle_circular(circ(n), store)
    got = reported(store.root, [("grb", "sb26111302")])
    assert got["circulars"] == [46244, 46261]
    assert [f["value"] for f in got["prompt"] if f["kind"] == "t90"][0] == 7.447
    assert got["prompt"][0]["circular"] == 46244 and got["prompt"][0]["subject"].startswith("GRB 261113A")
    z = next(f for f in got["spectroscopy"] if f["kind"] == "redshift")
    assert (z["value"], z["what"], z["circular"]) == (3.58, "spectroscopic", 46261)
    assert any(f["kind"] == "lines" for f in got["spectroscopy"])
    assert reported(store.root, [("grb", "nothing")])["circulars"] == []
    # a circular stored before the facts were read is read when asked for
    import json
    p = store.root / "circulars" / "46244.json"
    old = json.loads(p.read_text())
    del old["facts"], old["web_links"]
    p.write_text(json.dumps(old))
    assert reported(store.root, [("grb", "sb26111302")])["prompt"] == got["prompt"]
    assert Store(store.root).reported("grb:sb26111302,optical:none")["circulars"] == [46244, 46261]
