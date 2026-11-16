"""Step 1 (CIT login node): science segments and strain around each event.
raw : LIGO RAW frames (CAL-DELTAL_EXTERNAL_DQ converted to strain) + Virgo HoftOnline
hoft: LIGO H1/L1_HOFT_C00 (GDS-CALIB_STRAIN_CLEAN) + Virgo HoftOnline
Strain of the analysable stretch within +-MARGIN s is written to small frames."""
import json, sys, time
from raven_kn.targeted.datasource import LVKSource, LVK_PRESETS, write_manifest
MARGIN = 512.0
ev = json.load(open("events.json"))
V = {"frame": "HoftOnline", "chan": "V1:Hrec_hoft_16384Hz", "flag": "V1:ITF_SCIENCE:1"}
flags = {"H1": "H1:DMT-ANALYSIS_READY:1", "L1": "L1:DMT-ANALYSIS_READY:1", "V1": V["flag"]}
src = {
 "raw": dict(frame_types={"H1": "H1_R", "L1": "L1_R", "V1": V["frame"]},
             channels={"H1": "H1:CAL-DELTAL_EXTERNAL_DQ", "L1": "L1:CAL-DELTAL_EXTERNAL_DQ", "V1": V["chan"]},
             science_flags=flags, raw_calibration={"H1": {}, "L1": {}}),
 "hoft": dict(frame_types={"H1": "H1_HOFT_C00", "L1": "L1_HOFT_C00", "V1": V["frame"]},
              channels={"H1": "H1:GDS-CALIB_STRAIN_CLEAN", "L1": "L1:GDS-CALIB_STRAIN_CLEAN", "V1": V["chan"]},
              science_flags=flags),
}
for sid in (sys.argv[1:] or ev):
    e = ev[sid]
    a, b = e["t0"] + e["window"][0], e["t0"] + e["window"][1]
    for kind, kw in src.items():
        t = time.time()
        s = LVKSource(sample_rate=2048, datafind_host="datafind.ldas.cit:80", **kw)
        man = write_manifest(s, e["ifos"].split(","), a, b, f"manifest_{sid}_{kind}.json", margin=MARGIN,
                             extract_dir=f"frames_{sid}_{kind}")
        print(sid, kind, "%.0f s" % (time.time() - t),
              {i: (m.get("error") or [m["analysable"], len(m.get("original_frames", []))]) for i, m in man["ifos"].items()}, flush=True)
