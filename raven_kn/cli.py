"""Command line: ``raven-kn assess candidate.json --out DIR``
and ``raven-kn coinc candidate.json --superevent-far HZ --gw-skymap FILE``."""

from __future__ import annotations

import argparse
import json
import os

from .connection import load_candidate
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
    a.add_argument("candidate")
    a.add_argument("--out", default=None)
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

    c = sub.add_parser("coinc", help="joint FAR with a GW superevent")
    c.add_argument("candidate")
    c.add_argument("--superevent-id", default="S000000")
    c.add_argument("--superevent-far", type=float, required=True)
    c.add_argument("--gw-skymap")
    c.add_argument("--chirp-mass", type=float, default=None,
                   help="superevent chirp mass in Msun (distance-factor horizon); default BNS")

    n = sub.add_parser("nmma", help="NMMA light-curve fit -> t0 posterior")
    n.add_argument("candidate")
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
    n.add_argument("--save", help="write the updated candidate JSON here")
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
    g.add_argument("--f-low", type=float, default=30.0)
    g.add_argument("--subsegment", type=float, default=1024.0)
    g.add_argument("--n-slides", type=int, default=200)
    g.add_argument("--n-offsource", type=int, default=0,
                   help="off-source trials (same length, outside the window) per sub-segment")
    g.add_argument("--statistic", choices=["network", "coherent"], default="network",
                   help="incoherent network SNR, or coherent SNR with null-stream cut (as PyGRB)")
    g.add_argument("--no-gating", action="store_true")
    g.add_argument("--source", default="gwosc", choices=["gwosc", "files"])
    g.add_argument("--files", help="IFO=path[:channel],... (frame files need :channel); "
                                   "detectors without a file are read from GWOSC")
    g.add_argument("--ca-bundle", help="CA file if Python fails SSL to gwosc.org, "
                                       "e.g. /etc/pki/tls/certs/ca-bundle.crt")
    g.add_argument("--label", default=None)
    g.add_argument("--out", default="targeted_out")

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
    y.add_argument("--segment-url", default="https://segments.ligo.org")

    li = sub.add_parser("listen", help="GCN Kafka listener: GW alerts and GRB notices vs stored "
                                       "kilonova candidates (credentials: $GCN_CLIENT_ID, "
                                       "$GCN_CLIENT_SECRET)")
    li.add_argument("--store", required=True, help="candidate/alert/result directory")
    li.add_argument("--group-id", default="raven-kn",
                    help="Kafka consumer group: a restart resumes where it stopped ('' for none)")
    li.add_argument("--include-mock", action="store_true",
                    help="also process LVK mock superevents (MS...), e.g. to test the chain")
    li.add_argument("--grb-topics", nargs="*", default=None)
    li.add_argument("--domain", default=None, help="e.g. test.gcn.nasa.gov")
    li.add_argument("--use-combined-skymap", action="store_true",
                    help="when an alert carries a RAVEN coincidence (external_coinc), use its "
                         "combined GW x GRB map for the sky factor instead of the GW-only map")

    ad = sub.add_parser("add", help="add a candidate to the listener store and match it against "
                                    "GW alerts already received")
    ad.add_argument("candidate")
    ad.add_argument("--store", required=True)
    ad.add_argument("--find-host", action="store_true",
                    help="look up the host redshift/distance (NED, GLADE+) before matching")
    ad.add_argument("--use-combined-skymap", action="store_true")

    h = sub.add_parser("host", help="host galaxy, redshift and redshift-independent distance "
                                    "(NED, then GLADE+, then DES Y6; NED-D for nearby hosts)")
    h.add_argument("candidate")
    h.add_argument("--max-offset-kpc", type=float, default=50.0)
    h.add_argument("--overwrite", action="store_true",
                   help="replace a spectroscopic redshift already in the candidate")
    h.add_argument("--no-des", action="store_true", help="skip the DES Y6 Gold photo-z fallback")
    h.add_argument("--save", help="write the updated candidate JSON here")

    t = sub.add_parser("tns", help="TNS classification, redshift and spectra "
                                   "($TNS_BOT_ID, $TNS_BOT_NAME, $TNS_API_KEY)")
    t.add_argument("candidate")
    t.add_argument("--objname", required=True, help="e.g. AT2025ulz")
    t.add_argument("--save")

    args = p.parse_args(argv)
    if args.cmd in ("host", "tns"):
        from .connection import TNSConnector, find_host, save_candidate
        kn = load_candidate(args.candidate)
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
            save_candidate(kn, args.save)
            res["saved"] = args.save
        print(json.dumps(res, indent=2, default=float))
        return
    if args.cmd in ("listen", "add"):
        import logging
        from .listener import GRB_TOPICS, CandidateStore, add_candidate, listen
        store = CandidateStore(args.store)
        if args.cmd == "add":
            kn = load_candidate(args.candidate)
            if args.find_host:
                from .connection import find_host
                find_host(kn)
            print(json.dumps(add_candidate(kn, store,
                                           use_combined=args.use_combined_skymap),
                             indent=2, default=float))
            return
        logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
        listen(store, os.environ["GCN_CLIENT_ID"], os.environ["GCN_CLIENT_SECRET"],
               grb_topics=args.grb_topics if args.grb_topics is not None else GRB_TOPICS,
               domain=args.domain, use_combined=args.use_combined_skymap,
               on_result=lambda r: print(json.dumps(r, default=float), flush=True),
               group_id=args.group_id or None, include_mock=args.include_mock)
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
                  mchirp_max=args.mchirp_max, f_low=args.f_low, subsegment_s=args.subsegment,
                  n_slides=args.n_slides, n_offsource=args.n_offsource, statistic=args.statistic,
                  gate_threshold=None if args.no_gating else 50.0)
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
        else:
            fetch = gw
        res = run(cfg, fetch, args.out)
        print(json.dumps({"p_value": res.p_value,
                          "p_value_extrapolated": res.loudest.get("p_value_extrapolated"),
                          "bank_size": res.bank_size,
                          "loudest": res.loudest}, indent=2, default=float))
        return
    kn = load_candidate(args.candidate)
    if args.cmd == "nmma":
        from .connection import save_candidate
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
            save_candidate(kn, args.save)
            res["saved"] = args.save
        print(json.dumps(res, indent=2, default=float))
        return
    if args.cmd == "assess":
        res = assess(kn, args.out, _ranges(args.ranges), args.snr_threshold, args.p_min,
                     args.d_max_useful, model_sys_frac=args.model_sys_frac,
                     chirp_mass_msun=args.chirp_mass)
    else:
        post = posterior_for(kn)
        win = estimate_window(kn)
        res = coinc_far(kn, {"superevent_id": args.superevent_id, "far": args.superevent_far},
                        win, args.gw_skymap, post, chirp_mass_msun=args.chirp_mass)
    print(json.dumps(res, indent=2, default=float))


if __name__ == "__main__":
    main()
