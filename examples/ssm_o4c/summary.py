import glob, json
ev = json.load(open("events.json"))
print(f"{'run':38s} {'stat':>6s} {'rew':>6s} {'null':>5s} {'chi2':>5s}  {'H1':>5s} {'L1':>5s} {'V1':>5s}  {'Mc':>8s} {'dt_ms':>7s} {'off max':>7s} {'med':>5s} {'p_off':>6s}  gates")
for f in sorted(glob.glob("*_targeted_result.json")):
    r = json.load(open(f)); l = r["loudest"]; s = r["subsegments"][0]
    sid = f.split("/")[-1].split("_")[0]; off = sorted(s.get("offsource") or [0])
    ri = l.get("loudest_rho_ifo") or {}
    g = {k: len(v) for k, v in (l.get("gates") or {}).items() if v}
    print(f"{f.replace('_targeted_result.json',''):38s} {l.get('loudest_rho_coh') or l.get('loudest_rho_network') or l['loudest_rho_net']:6.2f} {l['loudest_rho_net']:6.2f} "
          f"{(l.get('loudest_rho_null') or 0):5.2f} {(l.get('loudest_chisq_r') or 0):5.2f}  {ri.get('H1',0):5.2f} {ri.get('L1',0):5.2f} {ri.get('V1',0):5.2f}  "
          f"{l['loudest_mchirp']:8.5f} {(l['loudest_gps']-ev[sid]['t0'])*1e3:7.1f} {off[-1]:7.2f} {off[len(off)//2]:5.2f} {(l.get('p_value_offsource') or 1):6.3f}  {g or ''}")
