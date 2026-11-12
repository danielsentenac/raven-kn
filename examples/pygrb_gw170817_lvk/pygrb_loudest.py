"""Loudest zero-lag trigger of each pycbc_multi_inspiral output, and the gates
it applied: python pygrb_loudest.py out/INSPIRAL_BANK*.hdf > pygrb_loudest.json"""
import json
import sys

import h5py
import numpy as np

out = {"per_bank": []}
for path in sorted(sys.argv[1:]):
    f = h5py.File(path, "r")
    n = f["network"]
    k = int(np.argmax(n["reweighted_snr"][:]))
    b = h5py.File("bank/" + path.split("_")[-1].replace(".hdf", "") + ".hdf", "r")
    tid = int(n["template_id"][k])
    m1, m2 = float(b["mass1"][tid]), float(b["mass2"][tid])
    row = {"file": path, "n_triggers": int(len(n["event_id"])),
           "coherent_snr": float(n["coherent_snr"][k]), "reweighted_snr": float(n["reweighted_snr"][k]),
           "null_snr": float(n["null_snr"][k]), "network_chisq": float(n["my_network_chisq"][k]),
           "end_time_gc": float(n["end_time_gc"][k]), "mass1": m1, "mass2": m2,
           "mchirp": (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2}
    for ifo in ("H1", "L1", "V1"):
        ev = int(n[f"{ifo}_event_id"][k])
        j = int(np.flatnonzero(f[ifo]["event_id"][:] == ev)[0])
        row[f"{ifo}_snr"] = float(f[ifo]["snr"][j])
        row[f"{ifo}_gates"] = [float(t) for t in f[ifo]["gating/auto/time"][:]]
    out["per_bank"].append(row)
out["loudest"] = max(out["per_bank"], key=lambda r: r["reweighted_snr"])
print(json.dumps(out, indent=1))
