"""Command line: ``raven-kn assess optical.json --out DIR``
and ``raven-kn coinc optical.json --superevent-far HZ --gw-skymap FILE``."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import os

from .connection import load_optical
from .search import assess, coinc_far
from .treatment.detectability import DEFAULT_RANGES_MPC
from .treatment.distance_inclination import posterior_for
from .treatment.t0 import estimate_window


def _ranges(s: str | None) -> dict[str, float]:
    if not s:
        return dict(DEFAULT_RANGES_MPC)
    return {k: float(v) for k, v in (kv.split("=") for kv in s.split(","))}


def _read_segments(path: str) -> list[tuple[float, float]]:
    out = []
    for line in open(path):
        f = line.split("#")[0].split()
        if len(f) >= 4:
            f = f[1:3]
        if len(f) >= 2:
            out.append((float(f[0]), float(f[1])))
    return out


def main(argv=None):
    p = argparse.ArgumentParser(prog="raven-kn")
    sub = p.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("assess", help="targeted-search trigger decision")
    a.add_argument("optical", help="JSON file of the optical transient")
    a.add_argument("--out", default=None)
    a.add_argument("--store", default=None,
                   help="write the decision in the logbook of this store (when it changed)")
    a.add_argument("--ranges", help="BNS ranges in Mpc, e.g. H1=150,L1=160,V1=55")
    a.add_argument("--snr-threshold", type=float, default=6.0)
    a.add_argument("--p-min", type=float, default=0.1)
    a.add_argument("--d-max-useful", type=float, default=500.0)
    a.add_argument("--chirp-mass", type=float, default=None,
                   help="source chirp mass in Msun to rescale the BNS ranges (e.g. 0.5 for a "
                        "sub-solar alert); default: BNS 1.2188")
    a.add_argument("--model-sys-frac", type=float, default=0.75,
                   help="model-systematic margin on a light-curve t0, as a fraction of the "
                        "extrapolation gap to the anchoring data (0 = statistical only)")
    a.add_argument("--no-kn-evidence", action="store_true",
                   help="trigger on detectability alone, without requiring evidence that the "
                        "transient is a kilonova (spectrum, light curve or classifier score)")

    c = sub.add_parser("coinc", help="joint FAR with a GW superevent")
    c.add_argument("optical", help="JSON file of the optical transient")
    c.add_argument("--superevent-id", default="S000000")
    c.add_argument("--superevent-far", type=float, required=True)
    c.add_argument("--gw-skymap")
    c.add_argument("--chirp-mass", type=float, default=None,
                   help="superevent chirp mass in Msun (distance-factor horizon); default BNS")

    n = sub.add_parser("nmma", help="NMMA light-curve fit -> t0 posterior")
    n.add_argument("optical", help="JSON file of the optical transient")
    n.add_argument("--workdir", required=True)
    n.add_argument("--model", default="Bu2019lm")
    n.add_argument("--model-prior", help="NMMA prior file for models without a built-in prior")
    n.add_argument("--svd-path")
    n.add_argument("--interpolation-type", default="tensorflow")
    n.add_argument("--sampler", default="pymultinest")
    n.add_argument("--max-age", type=float, default=7.0, help="days before first detection")
    n.add_argument("--fit-days", type=float, default=3.0,
                   help="fit only data up to this many days after the first detection "
                        "(<= 0: all data)")
    n.add_argument("--uninformed", action="store_true",
                   help="broad distance/inclination priors (fit result reusable as constraint)")
    n.add_argument("--dry-run", action="store_true", help="write inputs and print the command")
    n.add_argument("--import-only", action="store_true", help="read a finished run in --workdir")
    n.add_argument("--save", help="write the updated optical transient JSON here")
    n.add_argument("--executable", default=os.environ.get("RAVEN_KN_NMMA_EXE", "lightcurve-analysis"),
                   help="NMMA lightcurve-analysis path (default: $RAVEN_KN_NMMA_EXE or PATH), "
                        "e.g. ~/miniforge3/envs/nmma/bin/lightcurve-analysis")

    g = sub.add_parser("targeted", help="light targeted GW search (PyCBC) at a known position")
    g.add_argument("--request", help="<name>_targeted_search_request.json from assess")
    g.add_argument("--ra", type=float)
    g.add_argument("--dec", type=float)
    g.add_argument("--gps-start", type=float)
    g.add_argument("--gps-end", type=float)
    g.add_argument("--ifos", default="H1,L1,V1")
    g.add_argument("--mchirp-min", type=float, default=1.0)
    g.add_argument("--mchirp-max", type=float, default=1.5)
    g.add_argument("--mass-ratio", type=float, default=1.0,
                   help="m1/m2 of every template (default: equal masses)")
    g.add_argument("--spin1z", type=float, default=0.0, help="aligned spin of every template")
    g.add_argument("--spin2z", type=float, default=0.0)
    g.add_argument("--polarisation", choices=["free", "circular", "face_on", "face_away"],
                   help="coherent statistic: restrict to a circularly polarised signal (system seen "
                        "close to face-on and/or face-away, e.g. a GRB); default: from the request's "
                        "inclination prior, else free")
    g.add_argument("--distance-prior", action="store_true",
                   help="coherent statistic: rank candidates by the likelihood marginalised over the "
                        "request's EM distance and inclination prior (or --distance-range)")
    g.add_argument("--distance-range", help="EM luminosity distance interval in Mpc, as lo,hi")
    g.add_argument("--pe", choices=["em", "gw", "both"],
                   help="after the search, Bilby parameter estimation of the loudest candidate at the "
                        "EM sky position: em = EM distance and inclination priors, gw = distance "
                        "uniform in volume and isotropic orientation, both = the two, to compare")
    g.add_argument("--pe-result", metavar="FILE",
                   help="with --pe: take the candidate from this targeted result instead of searching")
    g.add_argument("--pe-nlive", type=int, default=500)
    g.add_argument("--pe-npool", type=int, default=1)
    g.add_argument("--pe-duration", type=float, default=64.0, help="analysed data length (s)")
    g.add_argument("--f-low", type=float, default=30.0)
    g.add_argument("--subsegment", type=float, default=1024.0)
    g.add_argument("--n-slides", type=int, default=200)
    g.add_argument("--n-offsource", type=int, default=0,
                   help="off-source trials (same length, outside the window) per sub-segment")
    g.add_argument("--statistic", choices=["network", "coherent"], default="network",
                   help="incoherent network SNR, or coherent SNR with null-stream cut (as PyGRB)")
    g.add_argument("--no-gating", action="store_true")
    g.add_argument("--source", default="auto", choices=["auto", "public", "lvk", "gwosc", "site", "replay"],
                   help="auto: public GWOSC data (datafind.gwosc.org frames, GWOSC timelines) if they "
                        "cover the window, else LVK data if configured and a SciToken is present; "
                        "public / lvk: force one; gwosc: GWOSC event files only, no segment check; "
                        "site: the online strain kept on the Virgo site; replay: the O4 replay kept "
                        "there, for tests")
    g.add_argument("--files", help="IFO=path[:channel],... (frame files need :channel); "
                                   "detectors without a file are read from GWOSC")
    g.add_argument("--lvk-preset", help="LVK frame types, channels and science flags of a run (O2)")
    g.add_argument("--frame-type", nargs="+", default=[], metavar="IFO=TYPE", help="LVK frame types")
    g.add_argument("--channel", nargs="+", default=[], metavar="IFO=CHANNEL", help="LVK h(t) channels")
    g.add_argument("--science-flag", nargs="+", default=[], metavar="IFO=FLAG",
                   help="LVK science flags, e.g. H1=H1:DMT-ANALYSIS_READY:1")
    g.add_argument("--veto-definer", help="LVK CBC veto-definer XML (CAT1 excluded, CAT2 reported)")
    g.add_argument("--datafind-host", help="LVK datafind server (default $GWDATAFIND_SERVER)")
    g.add_argument("--segment-url", default="https://segments.igwn.org")
    g.add_argument("--public-rate", type=int, choices=[4, 16],
                   help="public frames at 4 or 16 kHz (<IFO>_GWOSC_<run>_4KHZ_R1 / _16KHZ_R1)")
    g.add_argument("--write-manifest", metavar="FILE",
                   help="resolve segments and frame files for the window, write them to FILE and stop "
                        "(on a login node, where the token and datafind are)")
    g.add_argument("--extract-frames", metavar="DIR",
                   help="with --write-manifest: also write the strain of each detector around the "
                        "window to a small frame file in DIR (RAW frames, tape: too big to ship)")
    g.add_argument("--manifest", metavar="FILE",
                   help="read the data as resolved in FILE (--write-manifest): no network or credentials")
    g.add_argument("--frame-dir", help="with --manifest: directory holding the frame files (default: "
                                       "the current one, else their original paths)")
    g.add_argument("--ca-bundle", help="CA file if Python fails SSL to gwosc.org, "
                                       "e.g. /etc/pki/tls/certs/ca-bundle.crt")
    g.add_argument("--label", default=None)
    g.add_argument("--out", default="targeted_out")
    g.add_argument("--visibility", choices=["public", "lvk"], default=None,
                   help="logbook visibility of this run (default: public for public data, else lvk)")

    y = sub.add_parser("pygrb", help="PyGRB override configs (one per sub-segment) + bank")
    y.add_argument("--request", required=True, help="targeted search request JSON from assess")
    y.add_argument("--out", default="pygrb_configs")
    y.add_argument("--subsegment", type=float, default=1024.0)
    y.add_argument("--f-low", type=float, default=30.0)
    y.add_argument("--offsource-factor", type=float, default=6.0)
    y.add_argument("--base-config", nargs="+", default=["pygrb_base.ini"],
                   help="LVK PyGRB base configuration files, without the per-GRB file")
    y.add_argument("--event-template", help="per-GRB ini to adapt (e.g. LVK GRB170817A.ini); "
                                            "its [workflow] event options are replaced")
    y.add_argument("--bank-path", help="bank path as seen where the workflow runs")
    y.add_argument("--extra-override", nargs="*", default=[],
                   help="further section:option:value overrides")
    y.add_argument("--vetoes", help="text file of vetoed GPS intervals 'start end' (all detectors; "
                                    "segwizard columns 'id start end duration' also read); "
                                    "sub-segments are placed between them")
    y.add_argument("--veto-definer", help="veto-definer XML: its CAT flags (--veto-categories) are "
                                          "queried from DQSegDB and excluded")
    y.add_argument("--veto-categories", type=int, nargs="+", default=[1, 2])
    y.add_argument("--science-flag", nargs="+", default=[], metavar="IFO=FLAG",
                   help="exclude time outside science, e.g. H1=H1:DMT-ANALYSIS_READY:1")
    y.add_argument("--segment-url", default="https://segments.igwn.org")

    li = sub.add_parser("listen", help="GCN Kafka listener: GW alerts, GRB notices and neutrino alerts vs stored "
                                       "optical transients (credentials: $GCN_CLIENT_ID, "
                                       "$GCN_CLIENT_SECRET)")
    li.add_argument("--store", required=True, help="optical transient/alert/result directory")
    li.add_argument("--group-id", default="raven-kn",
                    help="Kafka consumer group: a restart resumes where it stopped ('' for none)")
    li.add_argument("--include-mock", action="store_true",
                    help="also process LVK mock superevents (MS...), e.g. to test the chain")
    li.add_argument("--grb-topics", nargs="*", default=None)
    li.add_argument("--nu-topics", nargs="*", default=None,
                    help="neutrino alert topics (default: IceCube gold/bronze tracks; none to disable)")
    li.add_argument("--no-circulars", action="store_true",
                    help="do not read the GCN Circulars (gcn.circulars)")
    li.add_argument("--domain", default=None, help="e.g. test.gcn.nasa.gov")
    li.add_argument("--use-combined-skymap", action="store_true",
                    help="when an alert carries a RAVEN coincidence (external_coinc), use its "
                         "combined GW x GRB map for the sky factor instead of the GW-only map")

    ci = sub.add_parser("circular", help="read GCN Circulars from the archive into the store: link them "
                                         "to its alerts and label them")
    ci.add_argument("ids", nargs="+", help="circular numbers or ranges, e.g. 46261 46270-46290")
    ci.add_argument("--store", required=True)

    ap = sub.add_parser("accept-photometry",
                        help="magnitudes reported in GCN Circulars about an optical transient: list them, "
                             "or copy them into its photometry once read (--accept)")
    ap.add_argument("name", help="name of the optical transient")
    ap.add_argument("--store", required=True)
    ap.add_argument("--circular", type=int, nargs="*", help="only those of these circulars")
    ap.add_argument("--accept", action="store_true", help="copy them (those that have a time)")

    ad = sub.add_parser("add", help="add an optical transient to the listener store and match it against "
                                    "GW alerts already received")
    ad.add_argument("optical", help="JSON file of the optical transient")
    ad.add_argument("--store", required=True)
    ad.add_argument("--find-host", action="store_true",
                    help="look up the host redshift/distance (NED, GLADE+) before matching")
    ad.add_argument("--use-combined-skymap", action="store_true")

    h = sub.add_parser("host", help="host galaxy, redshift and redshift-independent distance "
                                    "(NED, then GLADE+, then DES Y6; NED-D for nearby hosts)")
    h.add_argument("optical", help="JSON file of the optical transient")
    h.add_argument("--max-offset-kpc", type=float, default=50.0)
    h.add_argument("--overwrite", action="store_true",
                   help="replace a spectroscopic redshift already in the optical transient")
    h.add_argument("--no-des", action="store_true", help="skip the DES Y6 Gold photo-z fallback")
    h.add_argument("--save", help="write the updated optical transient JSON here")

    t = sub.add_parser("tns", help="TNS classification, redshift and spectra "
                                   "($TNS_BOT_ID, $TNS_BOT_NAME, $TNS_API_KEY)")
    t.add_argument("optical", help="JSON file of the optical transient")
    t.add_argument("--objname", required=True, help="e.g. AT2025ulz")
    t.add_argument("--save")

    z = sub.add_parser("ztffp", help="ZTF forced photometry from IPAC's ZFPS "
                                     "($ZTF_FPS_EMAIL, $ZTF_FPS_USERPASS)")
    z.add_argument("action", choices=["submit", "fetch"],
                   help="submit: queue the optical transient's position; fetch: merge the finished "
                        "light curve into the optical transient (results take minutes to days)")
    z.add_argument("optical", help="JSON file of the optical transient")
    z.add_argument("--days-before", type=float, default=60.0,
                   help="submit: JD range starts this long before the first detection (or now)")
    z.add_argument("--days-after", type=float, default=30.0,
                   help="submit: and ends this long after it (capped at now)")
    z.add_argument("--force", action="store_true",
                   help="submit even if the position was submitted in the last 90 days")
    z.add_argument("--save", help="fetch: write the updated optical transient JSON here")

    rl = sub.add_parser("runlog", help="mark a search or PE run in the logbook (e.g. submitted to a cluster)")
    rl.add_argument("out", help="output directory of the run")
    rl.add_argument("--request", required=True, help="the targeted-search request (gives the optical transient)")
    rl.add_argument("--label", default=None,
                    help="run name: the --label of the targeted command (default: the optical transient)")
    rl.add_argument("--stage", default="submitted", choices=["submitted"])
    rl.add_argument("--note", default="")
    rl.add_argument("--visibility", choices=["public", "lvk"], default="lvk")

    sv = sub.add_parser("serve", help="read-only HTTP service on a store, for the logbook application")
    sv.add_argument("--store", required=True)
    sv.add_argument("--host", default="127.0.0.1")
    sv.add_argument("--port", type=int, default=8765)
    sv.add_argument("--token-file", default=None,
                    help="file with the access token (else RAVEN_KN_TOKEN); required")
    sv.add_argument("--lvk-token-file", default=None,
                    help="file with the token that also shows runs on non-public data "
                         "(else RAVEN_KN_LVK_TOKEN)")
    sv.add_argument("--virgo-search-token-file", default=None,
                    help="file with the Virgo targeted-search token, required to start a targeted "
                         "search on the Virgo site (else RAVEN_KN_VIRGO_SEARCH_TOKEN; optional)")

    la = sub.add_parser("launcher", help="carry out the targeted searches asked for from the logbook "
                                         "application, on the data kept on the site")
    la.add_argument("--store", required=True)
    la.add_argument("--once", action="store_true", help="stop when no request is waiting")
    la.add_argument("--interval", type=float, default=10.0, help="seconds between two looks")

    fl = sub.add_parser("fink-listen", help="optical transients from a Fink Kafka livestream "
                                            "(credentials: finkctl auth register)")
    fl.add_argument("--survey", choices=["ztf", "lsst"], required=True)
    fl.add_argument("--store", required=True, help="the listener's store (shared with raven-kn listen)")
    fl.add_argument("--topics", nargs="*", default=None,
                    help="Fink topics (ZTF default: the three kilonova topics; Rubin: required, "
                         "see `finkctl topic list -survey lsst`)")
    fl.add_argument("--assess", action="store_true",
                    help="re-run assess on every update; decisions in <store>/assess/")
    fl.add_argument("--ranges", help="BNS ranges in Mpc for --assess, e.g. H1=150,L1=160,V1=55")
    fl.add_argument("--find-host", action="store_true",
                    help="look up the host galaxy once per new optical transient (needed for the "
                         "light-curve kilonova test)")
    fl.add_argument("--max-alerts", type=int, default=None)

    args = p.parse_args(argv)
    if args.cmd == "serve":
        from .serve import serve
        tok = (Path(args.token_file).read_text().strip() if args.token_file
               else os.environ.get("RAVEN_KN_TOKEN", ""))
        lvk = (Path(args.lvk_token_file).read_text().strip() if args.lvk_token_file
               else os.environ.get("RAVEN_KN_LVK_TOKEN", ""))
        if not tok and not lvk:
            p.error("an access token is required (--token-file or RAVEN_KN_TOKEN)")
        vts = (Path(args.virgo_search_token_file).read_text().strip() if args.virgo_search_token_file
               else os.environ.get("RAVEN_KN_VIRGO_SEARCH_TOKEN", ""))
        serve(args.store, args.host, args.port, tok or None, lvk or None, vts or None)
        return
    if args.cmd == "launcher":
        import logging
        from .launch import launcher
        from .listener import AlertStore
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
        launcher(AlertStore(args.store, source="launcher"), args.interval, args.once)
        return
    if args.cmd == "runlog":
        from .journal import RunLog
        req = json.loads(Path(args.request).read_text())
        cand = req.get("optical") or req.get("candidate", "unknown")
        RunLog(args.out, cand, args.label or cand,
               visibility=args.visibility).submitted(args.note)
        return
    if args.cmd == "fink-listen":
        import logging
        from .fink_stream import ZTF_KN_TOPICS, stream
        from .listener import AlertStore
        topics = args.topics if args.topics else (list(ZTF_KN_TOPICS) if args.survey == "ztf" else None)
        if not topics:
            p.error("--topics is required for --survey lsst (Fink has no Rubin kilonova topic yet)")
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
        n = stream(AlertStore(args.store, source=f"fink-{args.survey}"), args.survey, topics,
                   assess_kw={"ranges_mpc": _ranges(args.ranges)} if args.assess else None,
                   find_host=args.find_host, max_alerts=args.max_alerts,
                   on_result=lambda r: print(json.dumps(r, default=float), flush=True))
        print(json.dumps({"alerts_processed": n}))
        return
    if args.cmd == "ztffp":
        from astropy.time import Time
        from .connection import ZTFForcedPhotometryConnector, save_optical
        kn = load_optical(args.optical)
        zc = ZTFForcedPhotometryConnector(os.environ["ZTF_FPS_EMAIL"], os.environ["ZTF_FPS_USERPASS"])
        if args.action == "submit":
            now = float(Time.now().jd)
            fd = kn.first_detection
            ref = fd.mjd + 2400000.5 if fd is not None else now
            jd0, jd1 = ref - args.days_before, min(ref + args.days_after, now)
            codes = zc.submit([(kn.ra, kn.dec)], jd0, jd1, force=args.force)
            print(json.dumps({"submitted": kn.name, "ra": kn.ra, "dec": kn.dec, "jdstart": jd0,
                              "jdend": jd1, "http_status": codes}, indent=2))
        else:
            n = zc.add_to(kn)
            res = {"optical": kn.name, "points_added": n,
                   "status": "merged" if n else "no finished light curve for this position yet"}
            if args.save and n:
                save_optical(kn, args.save)
                res["saved"] = args.save
            print(json.dumps(res, indent=2))
        return
    if args.cmd in ("host", "tns"):
        from .connection import TNSConnector, find_host, save_optical
        kn = load_optical(args.optical)
        if args.cmd == "host":
            res = find_host(kn, max_offset_kpc=args.max_offset_kpc, overwrite=args.overwrite,
                            use_des=not args.no_des)
        else:
            res = TNSConnector(os.environ["TNS_BOT_ID"], os.environ["TNS_BOT_NAME"],
                               os.environ["TNS_API_KEY"]).add_to(kn, args.objname)
        res["redshift"] = kn.redshift.__dict__ if kn.redshift else None
        res["distance"] = kn.distance.__dict__ if kn.distance else None
        res["classification"] = kn.classification
        if args.save:
            save_optical(kn, args.save)
            res["saved"] = args.save
        print(json.dumps(res, indent=2, default=float))
        return
    if args.cmd == "accept-photometry":
        from .circulars import accept_photometry, reported_photometry
        from .listener import AlertStore
        store = AlertStore(args.store, source="cli")
        if args.accept:
            out = accept_photometry(store, args.name, args.circular)
        else:
            kn = next((t for t in store.transients() if t.name == args.name), None)
            if kn is None:
                raise SystemExit(f"no optical transient {args.name} in the store")
            out = [p for p in reported_photometry(store.root, kn)
                   if not args.circular or p["circular"] in args.circular]
        print(json.dumps(out, indent=2))
        return
    if args.cmd == "circular":
        from .circulars import expand_ids, fetch_circular, handle_circular
        from .listener import AlertStore
        store = AlertStore(args.store, source="cli")
        out = []
        for cid in expand_ids(args.ids):
            try:
                msg = fetch_circular(cid)
            except Exception as e:                           # a number not published, a network error
                out.append({"circular": cid, "error": f"{type(e).__name__}: {e}"})
                continue
            out.append({"circular": cid, "subject": msg.get("subject"),
                        "links": handle_circular(msg, store)})
        print(json.dumps(out, indent=2))
        return
    if args.cmd in ("listen", "add"):
        import logging
        from .listener import GRB_TOPICS, NU_TOPICS, AlertStore, add_optical, listen
        store = AlertStore(args.store, source="listen" if args.cmd == "listen" else "cli")
        if args.cmd == "add":
            kn = load_optical(args.optical)
            if args.find_host:
                from .connection import find_host
                find_host(kn)
            print(json.dumps(add_optical(kn, store,
                                           use_combined=args.use_combined_skymap),
                             indent=2, default=float))
            return
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
        listen(store, os.environ["GCN_CLIENT_ID"], os.environ["GCN_CLIENT_SECRET"],
               grb_topics=args.grb_topics if args.grb_topics is not None else GRB_TOPICS,
               nu_topics=args.nu_topics if args.nu_topics is not None else NU_TOPICS,
               domain=args.domain, use_combined=args.use_combined_skymap,
               on_result=lambda r: print(json.dumps(r, default=float), flush=True),
               group_id=args.group_id or None, include_mock=args.include_mock,
               circulars=not args.no_circulars)
        return
    if args.cmd == "pygrb":
        from .targeted.pygrb import write_pygrb_overrides
        from .targeted.search import TargetedSearchConfig
        cfg = TargetedSearchConfig.from_request(args.request, subsegment_s=args.subsegment,
                                                f_low=args.f_low)
        vetoes = _read_segments(args.vetoes) if args.vetoes else []
        if args.veto_definer or args.science_flag:
            from .targeted.segments import excluded_segments, write_segments
            vetoes += excluded_segments(
                cfg.ifos, cfg.gps_start, cfg.gps_end, veto_definer=args.veto_definer,
                categories=args.veto_categories, url=args.segment_url,
                science_flags=dict(f.split("=", 1) for f in args.science_flag))
            os.makedirs(args.out, exist_ok=True)
            write_segments(sorted(vetoes), os.path.join(args.out, "excluded_segments.txt"))
        files = write_pygrb_overrides(cfg, args.out, offsource_factor=args.offsource_factor,
                                      base_config=args.base_config,
                                      event_template=args.event_template,
                                      bank_path=args.bank_path,
                                      extra_overrides=args.extra_override,
                                      vetoes=vetoes)
        print(json.dumps({"n_subsegments": len(files), "out": args.out,
                          "excluded_s": sum(b - a for a, b in vetoes),
                          "run_script": f"{args.out}/run_pygrb.sh"}, indent=2))
        return
    if args.cmd == "targeted":
        from .targeted.search import (TargetedSearchConfig, file_fetcher, gwosc_fetcher, run)
        kw = dict(ifos=tuple(args.ifos.split(",")), mchirp_min=args.mchirp_min,
                  mchirp_max=args.mchirp_max, mass_ratio=args.mass_ratio, spin1z=args.spin1z,
                  spin2z=args.spin2z, f_low=args.f_low, subsegment_s=args.subsegment,
                  n_slides=args.n_slides, n_offsource=args.n_offsource, statistic=args.statistic,
                  gate_threshold=None if args.no_gating else 50.0)
        if args.distance_prior:
            kw["distance_prior"] = True
        if args.distance_range:
            kw["distance_range_mpc"] = tuple(float(x) for x in args.distance_range.split(","))
        if args.polarisation:
            kw["polarisation"] = args.polarisation
        if args.request:
            cfg = TargetedSearchConfig.from_request(args.request, **kw)
        else:
            cfg = TargetedSearchConfig(args.ra, args.dec, args.gps_start, args.gps_end, **kw)
        if args.label:
            cfg.label = args.label
        gw = gwosc_fetcher(cfg.sample_rate, args.ca_bundle)
        if args.files:
            files = dict(kv.split("=", 1) for kv in args.files.split(","))
            fetch = file_fetcher(files, cfg.sample_rate, fallback=gw)
        elif args.manifest:
            from .targeted.datasource import ManifestSource
            fetch = ManifestSource.load(args.manifest, cfg.sample_rate, args.frame_dir)
            print(f"data source: {fetch.name} (manifest {args.manifest})", flush=True)
        elif args.source in ("site", "replay"):
            from .targeted.datasource import SiteSource
            fetch = SiteSource.named("live" if args.source == "site" else "replay", cfg.sample_rate)
            print(f"data source: {fetch.name}", flush=True)
        elif args.source == "gwosc":
            fetch = gw
        else:
            from .targeted.datasource import LVKSource, PublicSource, select_source
            lvk = None
            if args.lvk_preset or args.frame_type:
                kw = {k: dict(x.split("=", 1) for x in v) for k, v in
                      (("frame_types", args.frame_type), ("channels", args.channel),
                       ("science_flags", args.science_flag)) if v}
                kw.update(sample_rate=cfg.sample_rate, veto_definer=args.veto_definer,
                          datafind_host=args.datafind_host, segment_url=args.segment_url)
                lvk = LVKSource.preset(args.lvk_preset, **kw) if args.lvk_preset else LVKSource(**kw)
            fetch = select_source(args.source, cfg.ifos, cfg.gps_start, cfg.gps_end,
                                  PublicSource(cfg.sample_rate, ca_bundle=args.ca_bundle,
                                               rate_khz=args.public_rate), lvk)
            print(f"data source: {fetch.name}", flush=True)
            if args.write_manifest:
                from .targeted.datasource import write_manifest
                man = write_manifest(fetch, cfg.ifos, cfg.gps_start, cfg.gps_end, args.write_manifest,
                                     extract_dir=args.extract_frames)
                print(json.dumps({"manifest": args.write_manifest, "source": man["source"],
                                  "span": man["span"],
                                  "frames": {i: len(m.get("frames", [])) for i, m in man["ifos"].items()},
                                  "errors": {i: m["error"] for i, m in man["ifos"].items() if "error" in m}},
                                 indent=2))
                return
        # logbook: <label>_state.json and <label>_journal.jsonl in the output directory (raven_kn.journal)
        from .journal import RunLog, pe_warnings
        src = getattr(fetch, "name", "gwosc" if fetch is gw else "files")
        vis = args.visibility or ("public" if src in ("public", "gwosc") else "lvk")
        req = json.loads(Path(args.request).read_text()) if args.request else {}
        cand = req.get("optical") or req.get("candidate") or cfg.label
        only_pe = bool(args.pe and args.pe_result)
        with RunLog(args.out, cand, cfg.label, visibility=vis,
                    kind="pe" if only_pe else "search") as rlog:
            rlog.state["data_source"] = src
            if only_pe:
                loudest = json.loads(Path(args.pe_result).read_text())["loudest"]
            else:
                rlog.stage("searching", f"{cfg.label}: search on {src} data, window "
                           f"[{cfg.gps_start:.2f}, {cfg.gps_end:.2f}], {cfg.statistic} statistic",
                           data_source=src)
                res = run(cfg, fetch, args.out, on_progress=lambda i, n, r: rlog.subsegment(
                    i, n, {"rho_net": r.loudest_rho_net, "gps": r.loudest_gps,
                           "mchirp": r.loudest_mchirp, "p_value": r.p_value}))
                loudest = res.loudest
                rlog.state["p_value"] = res.p_value
                rlog.log("search_result", f"{cfg.label}: loudest {loudest.get('loudest_rho_net', 0):.2f}"
                         f", p-value {res.p_value:.3g}, {res.bank_size} templates",
                         p_value=res.p_value, bank_size=res.bank_size, loudest=loudest,
                         result_file=f"{cfg.label}_targeted_result.json")
                print(json.dumps({"p_value": res.p_value,
                                  "p_value_extrapolated": res.loudest.get("p_value_extrapolated"),
                                  "bank_size": res.bank_size,
                                  "loudest": res.loudest}, indent=2, default=float))
            if args.pe and loudest.get("loudest_rho_net", 0) > 0:
                from .targeted.pe import run_pe
                ifos = list(loudest.get("loudest_rho_ifo") or cfg.ifos)
                for em in {"em": (True,), "gw": (False,), "both": (True, False)}[args.pe]:
                    which = "em" if em else "gw"
                    rlog.pe(which, "running")
                    out = run_pe(cfg, fetch, loudest["loudest_mchirp"], loudest["loudest_gps"], args.out,
                                 em=em, ifos=ifos, duration=args.pe_duration, nlive=args.pe_nlive,
                                 npool=args.pe_npool)
                    rlog.pe(which, "done", summary=out,
                            warnings=pe_warnings(out, loudest.get("loudest_rho_coh")
                                                 or loudest.get("loudest_rho_net")))
                    print(json.dumps(out, indent=2, default=float))
        return
    kn = load_optical(args.optical)
    if args.cmd == "nmma":
        from .connection import save_optical
        from .treatment import lightcurve_fit as lcf
        if args.import_only:
            lcf.import_result(kn, args.workdir)
            s = lcf.NMMASetup.load(args.workdir)
        else:
            kn, s = lcf.fit(kn, args.workdir, dry_run=args.dry_run, model=args.model,
                            model_prior=args.model_prior, svd_path=args.svd_path,
                            interpolation_type=args.interpolation_type, sampler=args.sampler,
                            max_age_d=args.max_age, informed=not args.uninformed,
                            fit_days=args.fit_days if args.fit_days > 0 else None,
                            executable=os.path.expanduser(args.executable))
        res = {"command": s.command_line, "trigger_mjd": s.trigger_mjd,
               "timeshift_bounds_d": s.timeshift_bounds, "informed_priors": s.informed,
               "n_detections": s.n_detections, "n_limits": s.n_limits, "notes": s.notes}
        if kn.t0_samples_gps:
            w = estimate_window(kn)
            res["t0_window"] = {"gps_start": w.gps_start, "gps_end": w.gps_end,
                                "duration_h": w.duration / 3600, "method": w.method,
                                "notes": w.notes}
        if args.save:
            save_optical(kn, args.save)
            res["saved"] = args.save
        print(json.dumps(res, indent=2, default=float))
        return
    if args.cmd == "assess":
        prev = None
        if args.out:
            pp = Path(args.out) / f"{kn.name.replace(' ', '_')}_decision.json"
            prev = json.loads(pp.read_text()) if pp.exists() else None
        res = assess(kn, args.out, _ranges(args.ranges), args.snr_threshold, args.p_min,
                     args.d_max_useful, model_sys_frac=args.model_sys_frac,
                     chirp_mass_msun=args.chirp_mass,
                     require_kn_evidence=not args.no_kn_evidence)
        if args.store:
            from .journal import Journal, log_decision
            log_decision(Journal.for_store(args.store, "cli"), res, prev)
    else:
        post = posterior_for(kn)
        win = estimate_window(kn)
        res = coinc_far(kn, {"superevent_id": args.superevent_id, "far": args.superevent_far},
                        win, args.gw_skymap, post, chirp_mass_msun=args.chirp_mass)
    print(json.dumps(res, indent=2, default=float))


if __name__ == "__main__":
    main()
