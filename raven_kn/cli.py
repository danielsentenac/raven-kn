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
    y.add_argument("--base-config", default="pygrb_base.ini",
                   help="LVK PyGRB base configuration the overrides are layered on")

    args = p.parse_args(argv)
    if args.cmd == "pygrb":
        from .targeted.pygrb import write_pygrb_overrides
        from .targeted.search import TargetedSearchConfig
        cfg = TargetedSearchConfig.from_request(args.request, subsegment_s=args.subsegment,
                                                f_low=args.f_low)
        files = write_pygrb_overrides(cfg, args.out, offsource_factor=args.offsource_factor,
                                      base_config=args.base_config)
        print(json.dumps({"n_subsegments": len(files), "out": args.out,
                          "run_script": f"{args.out}/run_pygrb.sh"}, indent=2))
        return
    if args.cmd == "targeted":
        from .targeted.search import (TargetedSearchConfig, file_fetcher, gwosc_fetcher, run)
        kw = dict(ifos=tuple(args.ifos.split(",")), mchirp_min=args.mchirp_min,
                  mchirp_max=args.mchirp_max, f_low=args.f_low, subsegment_s=args.subsegment,
                  n_slides=args.n_slides, gate_threshold=None if args.no_gating else 50.0)
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
