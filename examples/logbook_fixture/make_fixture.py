"""A store with its logbook, for developing and testing the logbook application.

Two cases, written with the pipeline's own handlers and a replay clock:

1. **GW170817** (real, as it happened): the Fermi-GBM notice, the LVK
   alert (reconstructed from the public record, with the two public sky
   maps), the Swope discovery of SSS17a and the AT2017gfo light curve, in
   the order they became public (``gw170817``). Then, dated November 2026,
   what RAVEN-KN did with it: the decision from the optical data and the
   targeted search of the public strain whose results are in
   ``examples/gw170817_reverse/`` (``gw170817_reprocessing``). The host
   redshift is NGC 4993 from NED, set by hand (no network).
2. **FIX-KN-1** (synthetic, to show what the real case does not): a real
   IceCube notice and real Swift-BAT / XRT notices from ``tests/data``, an
   invented candidate, invented superevents ``S261020fx`` (preliminary and
   initial maps, coincidence) and ``S261021fy`` (retracted), a failed run, a
   run on LVK data with PE summaries that fail the checks, and a listener
   outage. Every invented number is marked "synthetic" in its notes.

  python examples/logbook_fixture/make_fixture.py --out /path/to/store
  raven-kn serve --store /path/to/store --token-file token.txt
"""

from __future__ import annotations

import argparse
import base64
import io
import json
import shutil
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import healpy as hp
import numpy as np
from astropy.table import Table
from astropy.time import Time, TimeDelta

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "examples")]

from raven_kn.connection import OpticalTransient, PhotometryPoint, Redshift          # noqa: E402
from raven_kn.journal import Heartbeat, RunLog, log_decision, pe_warnings       # noqa: E402
from raven_kn.listener import (AlertStore, add_optical, handle_grb_notice,  # noqa: E402
                               handle_gw_alert, handle_neutrino_notice)
from raven_kn.search import assess                                               # noqa: E402
import gw170817_reverse as R                                                      # noqa: E402

LIVE = ROOT / "tests" / "data" / "gcn_live"
REV = ROOT / "examples" / "gw170817_reverse"


class Clock:
    def __init__(self):
        self.t = datetime(2017, 8, 17, tzinfo=timezone.utc)

    def __call__(self) -> datetime:
        self.t += timedelta(milliseconds=1)          # keeps the order of entries of one step
        return self.t

    def set(self, when) -> None:
        if isinstance(when, (int, float)):
            when = Time(when, format="mjd", scale="utc").to_datetime(timezone.utc)
        elif isinstance(when, str):
            when = datetime.fromisoformat(when.replace("Z", "+00:00"))
        self.t = max(self.t, when)                   # never backwards: entries keep their order


def gaussian_moc(ra, dec, sigma_deg, order=6, mu=120.0, sig=30.0) -> Table:
    nside = 1 << order
    npix = hp.nside2npix(nside)
    lon, lat = hp.pix2ang(nside, np.arange(npix), nest=True, lonlat=True)
    v0 = hp.ang2vec(ra, dec, lonlat=True)
    sep = np.degrees(np.arccos(np.clip(hp.ang2vec(lon, lat, lonlat=True) @ v0, -1, 1)))
    p = np.exp(-0.5 * (sep / sigma_deg) ** 2)
    p /= p.sum()
    return Table({"UNIQ": (4 * 4 ** order + np.arange(npix)).astype(np.int64),
                  "PROBDENSITY": p / hp.nside2pixarea(nside),
                  "DISTMU": np.full(npix, mu), "DISTSIGMA": np.full(npix, sig),
                  "DISTNORM": np.full(npix, 1 / (mu ** 2 + sig ** 2))})


SYNTHETIC_CLASS = ({"BNS": 0.9, "NSBH": 0.05, "BBH": 0.0, "Terrestrial": 0.05}, {"HasNS": 1.0, "HasRemnant": 0.9})


def gw_alert(sid, alert_type, iso_time, created, far=None, moc=None, pipeline="synthetic",
             instruments=("H1", "L1"), external_coinc=None, classification=SYNTHETIC_CLASS) -> dict:
    ev = None
    if alert_type != "RETRACTION":
        ev = {"time": iso_time, "far": far, "significant": True, "instruments": list(instruments),
              "group": "CBC", "pipeline": pipeline, "search": "AllSky"}
        if classification:
            ev["classification"], ev["properties"] = classification
        if moc is not None:
            buf = io.BytesIO()
            moc.write(buf, format="fits")
            ev["skymap"] = base64.b64encode(buf.getvalue()).decode()
    out = {"superevent_id": sid, "alert_type": alert_type, "time_created": created, "event": ev}
    if external_coinc:
        out["external_coinc"] = external_coinc
    return out


GW170817_ISO = "2017-08-17T12:41:04.43Z"
GW170817_FAR = 1 / (8.0e4 * 365.25 * 86400)     # less than one per 8.0e4 years (PRL 119, 161101)


def gw170817_alert(alert_type: str, name: str | None, created: str) -> dict:
    """GW170817 as an igwn.gwalert message. In 2017 the alert was a GCN
    circular, not a Kafka message, so the message is reconstructed from what
    was public: the merger time, the published FAR bound, the public sky
    maps (fetch_gw170817_maps.py; ``name`` None for the first circular, which
    is given here without a map), and the association with GRB 170817A
    (Fermi-GBM trigger 524666471) with no joint FAR. No source
    classification: none was issued in 2017."""
    moc = Table.read(Path(__file__).resolve().parent / f"GW170817_{name}.multiorder.fits") if name else None
    return gw_alert("GW170817", alert_type, GW170817_ISO, created, GW170817_FAR, moc,
                    pipeline="reconstructed from the public record" + (f" (map {name})" if name else ""),
                    instruments=("H1", "L1", "V1"), classification=None,
                    external_coinc={"gcn_notice_id": 524666471, "observatory": "Fermi", "search": "GRB",
                                    "time_difference": -2.04})


def gw170817(store: AlertStore, clock: Clock) -> None:
    """The event as it happened, in the order the messages became public
    (times of the GCN circulars, from memory of the record: to be checked
    against the GCN archive before being quoted):

    12:41:04 UTC  merger (H1, L1, V1)
    12:41:06      Fermi-GBM trigger 524666471; notice 14 s later
    13:21:42      LVK circular 21505: a BNS candidate coincident with the GRB
    17:54:51      LVK circular 21513: three-detector sky map, 31 deg2
    23:33:17      Swope image of SSS17a (AT 2017gfo) in NGC 4993;
                  circular 21529 on 18 Aug at 01:05:23
    18-27 Aug     the light curve (NMMA's reduced AT2017gfo file)
    16 Oct        publication, with the final sky map, 28 deg2
    """
    jr = store.journal
    clock.set("2017-08-17T12:41:20+00:00")
    handle_grb_notice(R.gbm_notice(), store)
    clock.set("2017-08-17T13:21:42+00:00")
    handle_gw_alert(gw170817_alert("PRELIMINARY", None, "2017-08-17T13:21:42Z"), store)
    clock.set("2017-08-17T17:54:51+00:00")
    handle_gw_alert(gw170817_alert("INITIAL", "bayestar", "2017-08-17T17:54:51Z"), store)

    disc = R.discovery_point()
    clock.set("2017-08-18T01:05:23+00:00")
    kn = OpticalTransient(name="SSS17a", ra=R.SSS17A["ra"], dec=R.SSS17A["dec"], pos_err_arcsec=0.5,
                          source="Swope (GCN 21529)", photometry=[disc])
    kn.redshift = Redshift(0.009783, 1e-5, "spectroscopic", "host NGC 4993 (NED)")
    kn.notes = "AT 2017gfo. host lookup: found."
    add_optical(kn, store, now_gps=float(Time("2017-08-18T01:05:23", scale="utc").gps))
    jr.log("host_lookup", "SSS17a: host NGC 4993, z = 0.009783", optical="SSS17a",
           host={"name": "NGC 4993"}, redshift=kn.redshift.__dict__)
    pending = R.followup_points()
    for te in sorted({p.mjd for p in pending}):
        clock.set(te)
        kn = next(c for c in store.transients() if c.name == "SSS17a")
        kn.photometry = [disc] + [p for p in pending if p.mjd <= te]
        add_optical(kn, store, now_gps=float(Time(te, format="mjd").gps))
    clock.set("2017-10-16T14:00:00+00:00")
    handle_gw_alert(gw170817_alert("UPDATE", "LALInference_v2", "2017-10-16T14:00:00Z"), store)


def gw170817_reprocessing(store: AlertStore, clock: Clock) -> None:
    """What RAVEN-KN itself did with the event, in November 2026: the
    decision from the optical data, and the targeted search of the public
    strain (``examples/gw170817_reverse/``, committed on 2026-11-13)."""
    jr = store.journal
    clock.set("2026-11-13T08:00:00+00:00")
    kn = next(c for c in store.transients() if c.name == "SSS17a")
    d = assess(kn, out_dir=store.root / "assess", ranges_mpc=R.O2_RANGES)
    log_decision(jr, d, None)
    assert d["trigger_targeted_search"], "the optical data should trigger the search"
    jr.log("search_request", "SSS17a: search request written, window "
           f"{d['t0_window']['gps_end'] - d['t0_window']['gps_start']:.0f} s (reprocessing of the public data)",
           optical="SSS17a",
           request=json.loads((store.root / "assess" / "SSS17a_targeted_search_request.json").read_text()))
    _sss17a_search(store, clock)


def _sss17a_search(store: AlertStore, clock: Clock) -> None:
    """The search of examples/gw170817_reverse (GWOSC data), with its 100 off-source trials."""
    run = "SSS17a_GRB_window_coherent"
    out = store.root / "runs" / "SSS17a" / run
    res = json.loads((REV / f"{run}_targeted_result.json").read_text())
    off = json.loads((REV / "offsource_merged.json").read_text())
    report = json.loads((REV / "report.json").read_text())
    for s in res["subsegments"]:
        s["offsource"] = off["offsource"]
        s["p_value_offsource"] = off["p_value"]
    res["loudest"]["p_value_offsource"] = off["p_value"]
    res["p_value"] = off["p_value"]
    clock.t += timedelta(minutes=20)
    rl = RunLog(out, "SSS17a", run, visibility="public", clock=clock, heartbeat_s=3600)
    rl.submitted("reprocessing of the public GWOSC data, coherent statistic, 100 off-source trials")
    clock.t += timedelta(minutes=5)
    with rl:
        rl.state["data_source"] = "public"
        rl.stage("searching", f"{run}: search on public data, window "
                 f"[{res['config']['gps_start']:.2f}, {res['config']['gps_end']:.2f}], coherent statistic",
                 data_source="public")
        clock.t += timedelta(seconds=report["search"]["runtime_s"])
        lo = res["loudest"]
        rl.subsegment(1, 1, {"rho_net": lo["loudest_rho_net"], "gps": lo["loudest_gps"],
                             "mchirp": lo["loudest_mchirp"], "p_value": res["p_value"]})
        (out / f"{run}_targeted_result.json").write_text(json.dumps(res))
        shutil.copy(REV / "snr_timeseries.json", out / f"{run}_snr_timeseries.json")
        rl.state["p_value"] = res["p_value"]
        rl.log("search_result", f"{run}: loudest {lo['loudest_rho_net']:.2f}, p-value "
               f"{res['p_value']:.3g}, {res['bank_size']} templates", p_value=res["p_value"],
               bank_size=res["bank_size"], loudest=lo, result_file=f"{run}_targeted_result.json")


def synthetic(store: AlertStore, clock: Clock) -> None:
    jr = store.journal
    SYN = "synthetic (logbook fixture)"
    # real notices from tests/data: Swift-BAT and its XRT position, IceCube bronze track
    clock.set("2026-10-19T08:00:00+00:00")
    for f in ("gcn.classic.voevent.SWIFT_BAT_GRB_POS_ACK.xml", "gcn.classic.voevent.SWIFT_XRT_POSITION.xml"):
        clock.t += timedelta(minutes=2)
        handle_grb_notice((LIVE / f).read_text(), store)
    nu = json.loads((LIVE / "gcn.notices.icecube.gold_bronze_track_alerts.001.json").read_text())
    clock.set(nu["trigger_time"].replace("Z", "+00:00"))
    clock.t += timedelta(seconds=40)
    handle_neutrino_notice(json.dumps(nu), store)
    t_nu = Time(nu["trigger_time"].replace("Z", ""), format="isot", scale="utc")

    # invented superevent 200 s after the neutrino, at its position: preliminary, then initial
    t_gw = t_nu + TimeDelta(200, format="sec")
    for k, (atype, sigma, far, dt) in enumerate((("PRELIMINARY", 7.0, 3e-8, 60), ("INITIAL", 3.5, 8e-9, 1500))):
        clock.set((t_gw + TimeDelta(dt, format="sec")).to_datetime(timezone.utc))
        handle_gw_alert(gw_alert("S261020fx", atype, t_gw.isot + "Z", clock.t.strftime("%Y-%m-%dT%H:%M:%SZ"),
                                 far, gaussian_moc(nu["ra"] + 2.0 - k, nu["dec"] - 1.5 + k, sigma)), store)

    # invented candidate found the next night, 0.3 deg from the neutrino
    mjd0 = t_nu.mjd + 0.9
    kn = OpticalTransient(name="FIX-KN-1", ra=nu["ra"] + 0.25, dec=nu["dec"] - 0.15, source=SYN, kn_score=0.55,
                     photometry=[PhotometryPoint(mjd0 - 1.6, "r", 21.8, upper_limit=True, survey=SYN),
                                 PhotometryPoint(mjd0, "r", 20.4, 0.08, survey=SYN)],
                     notes="synthetic candidate of the logbook fixture.")
    clock.set(mjd0 + 0.01)
    add_optical(kn, store, now_gps=float(Time(mjd0, format="mjd").gps))
    kn = next(c for c in store.transients() if c.name == "FIX-KN-1")
    jr.log("lookup_failed", "FIX-KN-1: host lookup failed (ConnectionError: NED unreachable)",
           level="warning", optical="FIX-KN-1", lookup="host", error="synthetic")
    kn.notes = (kn.notes + " host lookup: failed (ConnectionError).").strip()
    store.save(kn)
    prev = None
    for dt, mag in ((0.02, None), (1.0, 20.9), (2.1, 21.7)):
        clock.set(mjd0 + dt)
        if mag is not None:
            kn.photometry.append(PhotometryPoint(mjd0 + dt - 0.01, "r", mag, 0.1, survey=SYN))
        if dt >= 1.0 and kn.redshift is None:
            kn.redshift = Redshift(0.028, 1e-4, "spectroscopic", "host, " + SYN)
            kn.notes = kn.notes.replace(" host lookup: failed (ConnectionError).", " host lookup: found.")
            jr.log("host_lookup", "FIX-KN-1: host lookup, found", optical="FIX-KN-1",
                   redshift=kn.redshift.__dict__)
        add_optical(kn, store, now_gps=float(Time(mjd0 + dt, format="mjd").gps))
        d = assess(kn, out_dir=store.root / "assess")
        log_decision(jr, d, prev)
        prev = d
        if mag is None:
            # a second invented superevent, retracted
            clock.set("2026-10-21T14:10:00+00:00")
            handle_gw_alert(gw_alert("S261021fy", "PRELIMINARY", "2026-10-21T14:08:55Z", "2026-10-21T14:10:00Z",
                                     2e-7, gaussian_moc(40.0, -30.0, 12.0)), store)
            clock.set("2026-10-21T14:52:00+00:00")
            handle_gw_alert(gw_alert("S261021fy", "RETRACTION", None, "2026-10-21T14:52:00Z"), store)

    # runs: one failed on public data, one on LVK data with PE that fails the checks
    clock.set(mjd0 + 2.2)
    rl = RunLog(store.root / "runs" / "FIX-KN-1" / "FIX-KN-1_public", "FIX-KN-1", "FIX-KN-1_public",
                visibility="public", clock=clock, heartbeat_s=3600)
    rl.submitted("synthetic")
    try:
        with rl:
            rl.stage("searching", "FIX-KN-1_public: search on public data", data_source="public")
            clock.t += timedelta(seconds=40)
            raise RuntimeError("DataNotAvailable: the data of this window are not public yet "
                               "(H1, L1); use --source lvk")
    except RuntimeError:
        pass
    clock.t += timedelta(hours=1)
    run = "FIX-KN-1_lvk"
    rl = RunLog(store.root / "runs" / "FIX-KN-1" / run, "FIX-KN-1", run, visibility="lvk",
                clock=clock, heartbeat_s=3600)
    rl.submitted("synthetic; LVK frames")
    rng = np.random.default_rng(1)
    with rl:
        rl.stage("searching", f"{run}: search on lvk data", data_source="lvk")
        for i in range(1, 4):
            clock.t += timedelta(minutes=25)
            rl.subsegment(i, 3, {"rho_net": [6.1, 8.4, 6.6][i - 1], "gps": float(t_gw.gps) + 3 * i,
                                 "mchirp": 1.21, "p_value": [0.6, 0.03, 0.4][i - 1]})
        bg = np.round(5.2 + rng.gamma(2.0, 0.45, 200), 3).tolist()
        sub = {"gps_start": float(t_gw.gps) - 500, "gps_end": float(t_gw.gps) + 500, "ifos": ["H1", "L1"],
               "loudest_rho_net": 8.4, "loudest_gps": float(t_gw.gps) + 6, "loudest_mchirp": 1.21,
               "loudest_rho_ifo": {"H1": 5.6, "L1": 6.3}, "background": bg, "p_value": 0.03,
               "gates": {"H1": [], "L1": []}, "notes": "synthetic", "offsource": [],
               "p_value_offsource": None, "statistic": "coherent", "loudest_rho_coh": 8.6,
               "loudest_rho_null": 0.9}
        res = {"config": {"statistic": "coherent", "label": run, "ifos": ["H1", "L1"]}, "bank_size": 2148,
               "subsegments": [sub], "p_value": 0.03, "loudest": {k: v for k, v in sub.items()
                                                                  if k not in ("background", "offsource")},
               "data": {"source": "lvk", "note": "synthetic"}}
        (rl.dir / f"{run}_targeted_result.json").write_text(json.dumps(res))
        rl.state["p_value"] = 0.03
        rl.log("search_result", f"{run}: loudest 8.40, p-value 0.03, 2148 templates", p_value=0.03,
               bank_size=2148, loudest=res["loudest"], result_file=f"{run}_targeted_result.json")
        for which, (dmed, dlo, dhi, lbf, edge) in (("em", (118.0, 104.0, 133.0, 21.0, None)),
                                                    ("gw", (62.0, 60.5, 71.0, 96.0, "lower"))):
            rl.pe(which, "running")
            clock.t += timedelta(hours=2)
            n = 1500
            dist = (rng.normal(dmed, (dhi - dlo) / 3.3, n) if edge is None
                    else 60.0 + rng.exponential(dmed - 60.0, n))
            ang = np.abs(rng.normal(25 if which == "em" else 48, 8 if which == "em" else 22, n)) % 90
            summ = {"label": f"{run}_pe_{which}", "em_priors": which == "em", "ifos": ["H1", "L1"],
                    "candidate": {"mchirp": 1.21, "gps": float(t_gw.gps) + 6},
                    "maximum_likelihood": {"snr": 8.7 if which == "em" else 8.9},
                    "log_evidence": -40000.0, "log_bayes_factor": lbf, "n_samples": 6000,
                    "posterior": {"luminosity_distance": {"median": dmed, "lower": dlo, "upper": dhi}},
                    "distance_at_prior_edge": edge, "note": "synthetic"}
            (rl.dir / f"{run}_pe_{which}_pe_summary.json").write_text(json.dumps(summ))
            (rl.dir / f"{run}_pe_{which}_pe_samples.json").write_text(json.dumps(
                {"label": summ["label"], "em_priors": which == "em", "n_total": 6000,
                 "luminosity_distance": np.round(dist, 2).tolist(),
                 "viewing_angle_deg": np.round(ang, 2).tolist()}))
            rl.pe(which, "done", summary=summ, warnings=pe_warnings(summ, 8.6))


def circulars(store: AlertStore, clock: Clock) -> None:
    """Invented GCN Circulars about the alerts of :func:`synthetic`, in the
    format of ``gcn.circulars``: a burst followed up in X-rays and in the
    optical (upper limits, then a detection with a redshift), a neutrino
    with one optical search, a superevent with a counterpart report."""
    from raven_kn.circulars import handle_circular

    def send(number, when, subject, event, body):
        clock.set(when)
        t = datetime.fromisoformat(when)
        handle_circular({"circularId": number, "subject": subject, "eventId": event,
                         "submitter": "A. Observer at Fixture Observatory <observer@example.org>",
                         "createdOn": int(t.timestamp() * 1000), "submittedHow": "web",
                         "body": body + "\n\n[Invented for the logbook fixture.]"}, store)

    send(90001, "2026-10-17T16:21:00+00:00", "GRB 261017A: Swift detection of a burst", "GRB 261017A",
         "At 16:05:05 UT, the Swift Burst Alert Telescope (BAT) triggered and located GRB 261017A "
         "(trigger=1519536). Swift slewed immediately to the burst.")
    send(90002, "2026-10-17T18:02:00+00:00", "GRB 261017A: Swift-XRT refined analysis", "GRB 261017A",
         "We have analysed 4.2 ks of XRT data for GRB 261017A. The light curve can be modelled with a "
         "power-law decay with an index of 1.1.")
    send(90003, "2026-10-17T21:40:00+00:00", "GRB 261017A: NOT optical upper limits", "GRB 261017A",
         "We observed the field of GRB 261017A with the Nordic Optical Telescope, starting 4.9 h after "
         "the burst. No source is detected at the XRT position down to r > 23.1 (3 sigma).")
    send(90004, "2026-10-18T03:15:00+00:00", "GRB 261017A: GTC optical afterglow detection", "GRB 261017A",
         "We observed the field of GRB 261017A with the 10.4 m GTC, 11 h after the burst. We detect a "
         "faint source at the XRT position, with r = 24.3 +/- 0.2.")
    send(90005, "2026-10-18T22:30:00+00:00", "GRB 261017A: VLT/X-shooter spectroscopic redshift z = 0.42",
         "GRB 261017A", "We obtained a spectrum of the afterglow of GRB 261017A with X-shooter. From "
         "several emission lines of the host galaxy we derive a redshift z = 0.42.")
    send(90006, "2026-10-20T03:10:00+00:00", "IceCube-261020A - IceCube observation of a high-energy "
         "neutrino candidate track-like event", "IceCube-261020A",
         "IceCube detected a track-like event with a moderate probability of being of astrophysical origin.")
    send(90007, "2026-10-20T19:45:00+00:00", "IceCube-261020A: ZTF optical upper limits", "IceCube-261020A",
         "We observed the localization region of IceCube-261020A with ZTF. No new transient is detected "
         "down to r > 20.5 within the 90 % region.")
    send(90008, "2026-10-21T10:05:00+00:00", "LIGO/Virgo/KAGRA S261020fx: optical counterpart candidate "
         "FIX-KN-1", "S261020fx",
         "We report the discovery of a fast-fading optical transient, FIX-KN-1, in the localisation of "
         "S261020fx. On 2026-10-21T08:10:00 UT we measure r = 20.6 +/- 0.1 (AB). In a second image, "
         "taken on 2026-10-22T08:20:00 UT, the source has faded to r = 21.5 +/- 0.2, and i > 22.1.")


def listener_history(store: AlertStore) -> None:
    """Two starts of the GCN listener with 3 h without heartbeat between
    them, a Kafka error and a message that failed (synthetic)."""
    real_now = datetime.now(timezone.utc)
    sj = type(store.journal).for_store(store.root, "listen", clock=lambda: real_now - timedelta(hours=30))
    hb = Heartbeat(store.root, "listen", sj, ["igwn.gwalert", "gcn.classic.voevent.FERMI_GBM_GND_POS"])
    hb.message("igwn.gwalert", skipped=True)
    hb.kafka_error("KafkaError{code=_TRANSPORT,val=-195,str=\"kafka.gcn.nasa.gov:9092: Disconnected\"}")
    try:
        raise ValueError("not well-formed (invalid token): line 1, column 0")
    except ValueError as e:
        hb.failure("gcn.classic.voevent.SWIFT_XRT_POSITION", e)
    st = json.loads(hb.path.read_text())
    st["heartbeat"] = (real_now - timedelta(hours=3, minutes=4)).strftime("%Y-%m-%dT%H:%M:%S.000Z")
    hb.path.write_text(json.dumps(st))
    sj.clock = lambda: datetime.now(timezone.utc) - timedelta(minutes=4)
    Heartbeat(store.root, "listen", sj, st["topics"]).message("igwn.gwalert", skipped=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="store directory (emptied first)")
    args = ap.parse_args()
    out = Path(args.out)
    if out.exists():
        shutil.rmtree(out)
    clock = Clock()
    store = AlertStore(out, source="replay")
    store.journal.clock = clock
    gw170817(store, clock)
    synthetic(store, clock)
    circulars(store, clock)
    gw170817_reprocessing(store, clock)
    listener_history(store)
    for f in (out / "assess").glob("*_skymap.fits"):       # 95 MB each, not read by the service
        f.unlink()
    from raven_kn.journal import read_events
    ev, _, _ = read_events(out, limit=100000)
    kinds = {}
    for e in ev:
        kinds[e["kind"]] = kinds.get(e["kind"], 0) + 1
    print(json.dumps({"store": str(out), "entries": len(ev), "kinds": kinds}, indent=1))


if __name__ == "__main__":
    main()
