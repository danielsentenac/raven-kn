"""GW170817 GRB window, coherent statistic, one polarisation mode ($1), GWOSC frames."""
import sys, time
from raven_kn.targeted.search import TargetedSearchConfig, file_fetcher, run
mode = sys.argv[1]
files = {i: f"{i[0]}-{i}_GWOSC_4KHZ_R1-1187006835-4096.gwf:{i}:GWOSC-4KHZ_R1_STRAIN" for i in ("H1", "L1", "V1")}
cfg = TargetedSearchConfig.from_request("SSS17a_targeted_search_request.json", statistic="coherent",
                                        n_offsource=100, mchirp_min=1.19, mchirp_max=1.205, label="pol_" + mode)
print("from request:", cfg.polarisation, cfg.viewing_angle_max_deg, flush=True)
cfg.polarisation = mode
t = time.time()
r = run(cfg, file_fetcher(files, cfg.sample_rate), out_dir=".", progress=True)
print(f"{time.time() - t:.0f} s", flush=True)
