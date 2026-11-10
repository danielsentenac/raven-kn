"""Off-source trials of the GW170817 reverse replay, in chunks (cluster jobs),
and their merge into a p-value.

The on-source search (``examples/gw170817_reverse.py``) gives the loudest
statistic in the 6 s GRB window. Its background is the same search in
off-source windows of the same length around it (PyGRB-style). Each job runs
trials ``first .. first + n - 1`` of the sequence of
:func:`raven_kn.targeted.search.offsource_windows`, so chunks add up to the
same windows as one long run.

  # one chunk (HTCondor job)
  python gw170817_reverse_offsource.py run --request SSS17a_targeted_search_request.json \\
      --first 0 --n 5 --frames H1=H.gwf:H1:GWOSC-4KHZ_R1_STRAIN L1=... V1=... --out chunk_000.json
  # merge
  python gw170817_reverse_offsource.py merge --onsource report.json chunk_*.json --out offsource.json
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import replace
from pathlib import Path


def run_chunk(args) -> None:
    from raven_kn.targeted.search import (TargetedSearchConfig, build_bank, file_fetcher,
                                          gwosc_fetcher, offsource_windows, search_subsegment)
    over = {"statistic": args.statistic, "n_slides": 0}
    if args.mchirp:
        over["mchirp_min"], over["mchirp_max"] = args.mchirp
    cfg = TargetedSearchConfig.from_request(args.request, **over)
    length = cfg.gps_end - cfg.gps_start
    wins = offsource_windows(cfg, length, args.first + args.n)[args.first:]
    t = time.time()
    bank = build_bank(cfg.mchirp_min, cfg.mchirp_max, cfg.f_low, cfg.min_match,
                      cfg.sample_rate, cfg.max_templates)
    files = dict(f.split("=", 1) for f in args.frames) if args.frames else {}
    fetch = (file_fetcher(files, cfg.sample_rate) if files
             else gwosc_fetcher(cfg.sample_rate, args.ca_bundle))
    trials = []
    for k, (a, b) in enumerate(wins, start=args.first):
        r = search_subsegment(cfg, bank, fetch, a, b)
        trials.append({"index": k, "gps_start": a, "gps_end": b, "ifos": r.ifos,
                       "loudest": r.loudest_rho_net, "loudest_gps": r.loudest_gps,
                       "rho_ifo": r.loudest_rho_ifo, "gates": r.gates})
        print(f"trial {k} [{a:.2f}, {b:.2f}] {r.loudest_rho_net:.2f} {r.ifos} {r.notes[:80]}",
              flush=True)
    Path(args.out).write_text(json.dumps({"statistic": cfg.statistic, "bank_size": len(bank),
                                          "mchirp": [cfg.mchirp_min, cfg.mchirp_max],
                                          "trials": trials, "runtime_s": time.time() - t},
                                         indent=1, default=float))


def merge(args) -> None:
    from raven_kn.treatment.far import targeted_pvalue
    on = json.loads(Path(args.onsource).read_text())["search"]
    trials = []
    for f in args.chunks:
        trials += json.loads(Path(f).read_text())["trials"]
    idx = [t["index"] for t in trials]
    assert len(idx) == len(set(idx)), "overlapping chunks"
    bkg = sorted(t["loudest"] for t in trials if len(t["ifos"]) >= 2)
    p = targeted_pvalue(on["loudest"], bkg)
    out = {"onsource_loudest": on["loudest"], "onsource_gps": on["loudest_gps"],
           "n_trials": len(bkg), "missing_indices": sorted(set(range(max(idx) + 1)) - set(idx)),
           "offsource_max": bkg[-1], "offsource_median": bkg[len(bkg) // 2],
           "offsource_90pc": bkg[int(0.9 * (len(bkg) - 1))],
           "n_louder": sum(b >= on["loudest"] for b in bkg), "p_value": p,
           "p_value_note": f"p = (n_louder + 1) / (n_trials + 1); floor 1/{len(bkg) + 1}",
           "offsource": bkg}
    Path(args.out).write_text(json.dumps(out, indent=1, default=float))
    print(f"on-source {on['loudest']:.2f}; {len(bkg)} off-source trials: max {bkg[-1]:.2f}, "
          f"median {out['offsource_median']:.2f}; p = {p:.4f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--request", required=True)
    r.add_argument("--first", type=int, required=True)
    r.add_argument("--n", type=int, required=True)
    r.add_argument("--statistic", default="coherent")
    r.add_argument("--mchirp", type=float, nargs=2, default=None)
    r.add_argument("--frames", nargs="*", help="IFO=path:channel or IFO=file.hdf5")
    r.add_argument("--ca-bundle", default=None)
    r.add_argument("--out", required=True)
    m = sub.add_parser("merge")
    m.add_argument("--onsource", required=True, help="report.json of gw170817_reverse.py")
    m.add_argument("chunks", nargs="+")
    m.add_argument("--out", required=True)
    args = ap.parse_args()
    (run_chunk if args.cmd == "run" else merge)(args)


if __name__ == "__main__":
    main()
