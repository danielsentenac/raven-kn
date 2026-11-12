"""Science segments and CAT1/CAT2 vetoes from DQSegDB (IGWN SciToken) for the
LVK-data PyGRB test: what makes the collaboration's h(t) frames analysable."""
import argparse
import json

from gwpy.segments import DataQualityFlag

from raven_kn.targeted.segments import veto_definer_segments

ap = argparse.ArgumentParser()
ap.add_argument("--start", type=float, required=True)
ap.add_argument("--end", type=float, required=True)
ap.add_argument("--veto-definer", required=True)
ap.add_argument("--science", nargs="+", required=True, help="IFO=FLAG:VERSION")
ap.add_argument("--url", default="https://segments.igwn.org")
ap.add_argument("--out", required=True)
a = ap.parse_args()
out = {"gps": [a.start, a.end], "url": a.url, "veto_definer": a.veto_definer.split("/")[-1], "ifos": {}}
for spec in a.science:
    ifo, flag = spec.split("=", 1)
    q = DataQualityFlag.query(flag, a.start, a.end, url=a.url)
    out["ifos"][ifo] = {
        "science_flag": flag,
        "science": [[float(x), float(y)] for x, y in q.active],
        **{f"cat{c}": veto_definer_segments(a.veto_definer, [ifo], a.start, a.end, categories=(c,), url=a.url)
           for c in (1, 2)}}
    print(ifo, out["ifos"][ifo])
json.dump(out, open(a.out, "w"), indent=1)
