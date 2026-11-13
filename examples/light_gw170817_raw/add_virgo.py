"""Three-detector manifest: H1 and L1 strain derived from the RAW frames at CIT
(manifest_raw.json) plus Virgo's h(t) cut out at CC-IN2P3 (virgo_lyon.sh).
Virgo's science segments and vetoes come from DQSegDB (IGWN SciToken):

  python add_virgo.py VETO_DEFINER.xml
"""
import json
import shutil
import sys
from pathlib import Path

from raven_kn.targeted.datasource import LVKSource, add_to_manifest

entry = json.loads(Path("virgo_entry.json").read_text())
x, y = entry["extracted"]
name = entry["frames"][0].rsplit("/", 1)[-1]
Path("frames_raw3").mkdir(exist_ok=True)
for f in list(Path("frames_raw").glob("*.gwf")) + [Path("frames_virgo") / name]:
    shutil.copy(f, Path("frames_raw3") / f.name)
entry["frames"] = ["file://" + str((Path("frames_raw3") / name).resolve())]
shutil.copy("manifest_raw.json", "manifest_raw3.json")
src = LVKSource.preset("O2", veto_definer=sys.argv[1])
man = add_to_manifest("manifest_raw3.json", "V1", entry, src.analysable("V1", x, y), src.vetoed("V1", x, y),
                      note="Virgo: reconstructed h(t) from the CC-IN2P3 tape archive (its raw frames hold no strain)")
print(json.dumps({i: {k: m[k] for k in ("channel", "analysable", "cat2", "original_channel")}
                  for i, m in man["ifos"].items()}, indent=1))
