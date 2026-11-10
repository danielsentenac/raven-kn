"""GW170817 replayed in the reverse direction: GRB alert and optical follow-up
first, kilonova confirmation, GW search last.

What was known on 17-19 August 2017, in the order it arrived, without the
LIGO/Virgo alert:

1. Fermi-GBM trigger 524666471 at 12:41:06.47 UTC, ground position
   RA 176.8, Dec -39.8, 11.6 deg (1 sigma, statistical; GCN 21506, 21520).
   The listener stores the notice; no kilonova candidate exists yet.
2. SSS17a / AT2017gfo found by Swope at 2017-08-17 23:33 UTC, i = 16.0
   (GCN 21529), 24 deg from the GBM position. Added as a candidate, it is
   matched to the stored GRB notice: its merger-time window becomes RAVEN's
   [-5, +1] s around the GRB. Host NGC 4993 (NED, NED-D).
3. ``assess`` after every new epoch of photometry (the real AT2017gfo light
   curve, NMMA's reduced file): the search is held until the light curve
   shows a kilonova (absolute magnitude in the kilonova range and fast fading
   or reddening, from >= 2 epochs). No spectroscopic classification is used,
   so the light curve alone has to do it. BNS sensitivity of that day:
   horizons H1 107, L1 218, V1 58 Mpc (LIGO-Virgo 2017, PRL 119, 161101),
   i.e. ranges = horizon / 2.264.
4. Targeted search of the GWOSC strain in the 6 s GRB window, coherent
   statistic, background from off-source trials (6 s windows around the
   on-source, as PyGRB).

  python examples/gw170817_reverse.py --workdir examples/gw170817_reverse \\
      --n-offsource 100 --ca-bundle /etc/pki/tls/certs/ca-bundle.crt
"""

from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import numpy as np
from astropy.time import Time

from raven_kn.connection import KNCandidate, PhotometryPoint, Redshift, find_host
from raven_kn.listener import CandidateStore, add_candidate, handle_grb_notice
from raven_kn.search import assess

ROOT = Path(__file__).resolve().parents[1]
GBM_TEMPLATE = ROOT / "tests" / "data" / "gcn_live" / "gcn.classic.voevent.FERMI_GBM_GND_POS.msg"
GFO_DATA = ROOT / "examples" / "validation" / "AT2017gfo_reduced.dat"

GRB_ISO = "2017-08-17T12:41:06.47"
GBM = {"TrigID": "524666471", "ra": 176.8, "dec": -39.8, "err": 11.6}
SWOPE_ISO = "2017-08-17T23:33:17"
SSS17A = {"ra": 197.450371, "dec": -23.381486}          # 13:09:48.089 -23:22:53.35 (GCN 21529)
O2_HORIZONS = {"H1": 107.0, "L1": 218.0, "V1": 58.0}
O2_RANGES = {k: round(v / 2.264, 1) for k, v in O2_HORIZONS.items()}
GW170817_GPS = 1187008882.43                            # truth, used only in the report


def gbm_notice() -> str:
    """The real GRB 170817A ground position in the format of a live GBM
    ground-position VOEvent (the 2017 notice is not archived by GCN Kafka)."""
    xml = GBM_TEMPLATE.read_text()
    xml = re.sub(r'ivorn="[^"]*"', f'ivorn="ivo://nasa.gsfc.gcn/Fermi#GBM_Gnd_Pos_{GRB_ISO}_'
                 f'{GBM["TrigID"]}_1-000"', xml, count=1)
    xml = re.sub(r'(name="TrigID"\s+value=")[^"]*', rf'\g<1>{GBM["TrigID"]}', xml)
    xml = re.sub(r"<ISOTime>[^<]*</ISOTime>", f"<ISOTime>{GRB_ISO}</ISOTime>", xml)
    xml = re.sub(r"<C1>[^<]*</C1>", f"<C1>{GBM['ra']:.4f}</C1>", xml)
    xml = re.sub(r"<C2>[^<]*</C2>", f"<C2>{GBM['dec']:.4f}</C2>", xml)
    xml = re.sub(r"<Error2Radius>[^<]*</Error2Radius>",
                 f"<Error2Radius>{GBM['err']:.4f}</Error2Radius>", xml)
    return xml


def discovery_point() -> PhotometryPoint:
    return PhotometryPoint(mjd=float(Time(SWOPE_ISO, scale="utc").mjd), band="i", mag=16.0,
                           mag_err=0.1, survey="Swope (GCN 21529)")


def followup_points() -> list[PhotometryPoint]:
    """The real AT2017gfo light curve (NMMA example file), detections only."""
    rows = [l.split() for l in GFO_DATA.read_text().splitlines() if l.strip()]
    t = Time([r[0] for r in rows], format="isot", scale="utc").mjd
    return sorted((PhotometryPoint(float(ti), r[1], float(r[2]), float(r[3]), survey="literature")
                   for ti, r in zip(t, rows) if np.isfinite(float(r[3]))), key=lambda p: p.mjd)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", default="examples/gw170817_reverse")
    ap.add_argument("--n-offsource", type=int, default=100)
    ap.add_argument("--statistic", default="coherent", choices=["coherent", "network"])
    ap.add_argument("--mchirp", type=float, nargs=2, default=None,
                    help="chirp-mass range (default: the search request's, 1.0-1.5)")
    ap.add_argument("--ca-bundle", default=None)
    ap.add_argument("--skip-search", action="store_true")
    args = ap.parse_args()
    work = Path(args.workdir)
    work.mkdir(parents=True, exist_ok=True)
    t_grb = Time(GRB_ISO, scale="utc")
    report = {"timeline": []}

    def log(t_mjd, what, **kw):
        h = (t_mjd - t_grb.mjd) * 24
        print(f"[+{h:6.1f} h] {what}")
        report["timeline"].append({"hours_after_grb": round(h, 2), "event": what, **kw})

    # 1. GRB notice first: stored, nothing to attach to yet
    store = CandidateStore(work / "store")
    hits = handle_grb_notice(gbm_notice(), store)
    log(t_grb.mjd, f"Fermi-GBM notice {GBM['TrigID']} stored; candidates updated: {hits}")

    # 2. Swope discovery ~11 h later: matched to the stored GRB; host galaxy
    kn = KNCandidate(name="SSS17a", ra=SSS17A["ra"], dec=SSS17A["dec"], pos_err_arcsec=0.5,
                     source="Swope (GCN 21529)", photometry=[discovery_point()])
    t_disc = discovery_point().mjd
    add_candidate(kn, store, now_gps=Time(SWOPE_ISO, scale="utc").gps)
    kn = next(c for c in store.candidates() if c.name == "SSS17a")
    log(t_disc, f"SSS17a discovered (i = 16.0); GRB attached: {kn.grb is not None}",
        notes=kn.notes)
    try:
        host = find_host(kn)
        hs = f"host {host['host']['name']}, z = {kn.redshift.z:.4f}"
    except Exception as e:                      # offline: NGC 4993 from NED by hand
        kn.redshift = Redshift(0.009783, 1e-5, "spectroscopic", "host NGC 4993 (NED, offline)")
        host, hs = {"error": f"{type(e).__name__}: {e}"}, "host lookup failed: NGC 4993 by hand"
    store.save(kn)
    log(t_disc, hs, host=host)

    # 3. assess at discovery and after each new epoch, until confirmed
    pending = followup_points()
    epochs = [t_disc] + sorted({p.mjd for p in pending})
    decision, t_trigger = None, None
    for te in epochs:
        kn.photometry = [discovery_point()] + [p for p in pending if p.mjd <= te]
        decision = assess(kn, ranges_mpc=O2_RANGES)
        ev = decision["kilonova_evidence"]
        det = decision["detectability"]
        log(te, f"assess with {len(kn.photometry)} points: trigger {decision['trigger_targeted_search']}"
            f" (expected SNR {det['snr_median']:.0f}; kilonova: {ev['confirmed']}"
            f"{' by ' + ', '.join(ev['by']) if ev['by'] else ''}) - {ev['reasons'][-1]}",
            trigger=decision["trigger_targeted_search"], kilonova_evidence=ev)
        if decision["trigger_targeted_search"]:
            t_trigger = te
            break
    if t_trigger is None:
        print("never confirmed: no search")
        (work / "report.json").write_text(json.dumps(report, indent=1, default=str))
        return
    decision = assess(kn, out_dir=work / "assess", ranges_mpc=O2_RANGES)
    req = work / "assess" / f"{kn.name}_targeted_search_request.json"
    r = json.loads(req.read_text())
    inside = r["gps_start"] <= GW170817_GPS <= r["gps_end"]
    log(t_trigger, f"search request: window [{r['gps_start']:.2f}, {r['gps_end']:.2f}] "
        f"({r['gps_end'] - r['gps_start']:.0f} s), Mc {r['mchirp_range']}; contains GW170817: {inside}")
    report["window"] = {"gps_start": r["gps_start"], "gps_end": r["gps_end"],
                        "contains_gw170817": inside, "trigger_hours_after_grb":
                        (t_trigger - t_grb.mjd) * 24,
                        "expected_network_snr": decision["detectability"]["snr_median"],
                        "ranges_mpc": O2_RANGES}
    (work / "report.json").write_text(json.dumps(report, indent=1, default=str))
    if args.skip_search:
        return

    # 4. targeted search
    from raven_kn.targeted.search import TargetedSearchConfig, gwosc_fetcher, run
    over = dict(statistic=args.statistic, n_offsource=args.n_offsource,
                label=f"SSS17a_GRB_window_{args.statistic}")
    if args.mchirp:
        over["mchirp_min"], over["mchirp_max"] = args.mchirp
    cfg = TargetedSearchConfig.from_request(req, **over)
    t = time.time()
    res = run(cfg, gwosc_fetcher(cfg.sample_rate, args.ca_bundle), out_dir=work / "search")
    s = res.subsegments[0]
    off = sorted(s.offsource)
    print(f"[search] {res.bank_size} templates, {len(off)} off-source trials, {time.time() - t:.0f} s")
    print(f"    loudest {s.loudest_rho_net:.2f} at {s.loudest_gps:.3f} "
          f"({s.loudest_gps - GW170817_GPS:+.3f} s from GW170817), Mc {s.loudest_mchirp:.4f}, "
          f"per IFO {s.loudest_rho_ifo}")
    print(f"    off-source max {off[-1] if off else None}, p = {res.p_value}")
    report["search"] = {"bank_size": res.bank_size, "statistic": args.statistic,
                        "loudest": s.loudest_rho_net, "loudest_gps": s.loudest_gps,
                        "dt_from_gw170817_s": s.loudest_gps - GW170817_GPS,
                        "mchirp": s.loudest_mchirp, "rho_ifo": s.loudest_rho_ifo,
                        "n_offsource": len(off), "offsource_max": off[-1] if off else None,
                        "offsource_median": off[len(off) // 2] if off else None,
                        "p_value": res.p_value, "notes": s.notes, "gates": s.gates,
                        "runtime_s": time.time() - t}
    (work / "report.json").write_text(json.dumps(report, indent=1, default=str))


if __name__ == "__main__":
    main()
