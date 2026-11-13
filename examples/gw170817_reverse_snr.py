"""Matched-filter SNR time series of the GW170817 reverse-replay search, for the loudest template.

Same conditioning, PSD, template and coherent combination as the search
(``raven_kn.targeted.search``), for the one template that the full-bank search
found loudest. Writes the per-detector and coherent SNR over the 6 s GRB
window (peak in 10 ms bins) and at full rate around the peak.

  python examples/gw170817_reverse_snr.py \\
      --request examples/gw170817_reverse/SSS17a_targeted_search_request.json \\
      --result examples/gw170817_reverse/SSS17a_GRB_window_coherent_targeted_result.json \\
      --frames H1=H-H1_GWOSC_4KHZ_R1-1187006835-4096.gwf:H1:GWOSC-4KHZ_R1_STRAIN ... \\
      --out examples/gw170817_reverse/snr_timeseries.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from raven_kn.targeted import search as S

GW170817_GPS = 1187008882.43


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--request", required=True)
    ap.add_argument("--result", required=True, help="targeted result of the full-bank search")
    ap.add_argument("--frames", nargs="*", default=[], help="IFO=path:channel (default: GWOSC)")
    ap.add_argument("--ca-bundle", default=None)
    ap.add_argument("--zoom-s", type=float, default=0.06)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    import pycbc.psd
    from pycbc.filter import matched_filter_core
    from pycbc.waveform import get_fd_waveform

    cfg = S.TargetedSearchConfig.from_request(args.request, statistic="coherent")
    loud = json.loads(Path(args.result).read_text())["loudest"]
    start, end, mc = loud["gps_start"], loud["gps_end"], loud["loudest_mchirp"]
    files = dict(f.split("=", 1) for f in args.frames)
    fetch = (S.file_fetcher(files, cfg.sample_rate) if files
             else S.gwosc_fetcher(cfg.sample_rate, args.ca_bundle))

    # data stretch as in search_subsegment, sized for the whole bank
    fs = cfg.sample_rate
    lead = S.template_duration(cfg.mchirp_min, cfg.f_low) + 4 * cfg.pad_s
    ifos, z, sig, gates = [], [], [], {}
    L = int(round((end - start) * fs))
    delays = S._delays(cfg.ifos, cfg.ra_deg, cfg.dec_deg, 0.5 * (start + end))
    m = S._mass_equal(mc)
    for ifo in cfg.ifos:
        ts = fetch(ifo, start - lead, end + 4 * cfg.pad_s)
        if len(ts) % 2:
            ts = ts[:len(ts) - 1]
        ts, gates[ifo] = S._condition(ts, cfg.f_low, cfg.gate_threshold, cfg.gate_width_s,
                                      cfg.gate_taper_s)
        ts = S._taper_ends(ts, cfg.pad_s)
        psd = pycbc.psd.welch(ts, seg_len=int(16 * fs), seg_stride=int(8 * fs), avg_method="median")
        st = ts.to_frequencyseries()
        psd = pycbc.psd.interpolate(psd, st.delta_f)
        if len(psd) != len(st):
            psd.resize(len(st))
        psd = pycbc.psd.inverse_spectrum_truncation(psd, int(4 * fs), low_frequency_cutoff=cfg.f_low)
        hp, _ = get_fd_waveform(approximant="TaylorF2", mass1=m, mass2=m, delta_f=st.delta_f,
                                f_lower=cfg.f_low)
        hp.resize(len(st))
        x = (start + delays[ifo] - float(st.start_time)) * fs
        i0 = int(round(x))
        hp = hp * np.exp(-2j * np.pi * hp.sample_frequencies.numpy() * (x - i0) / fs)
        snr_u, _, norm = matched_filter_core(hp, st, psd=psd, low_frequency_cutoff=cfg.f_low)
        z.append(snr_u.numpy()[i0:i0 + L] * norm)       # aligned to geocentre time
        sig.append(4 * st.delta_f / norm)
        ifos.append(ifo)

    ant = S._antenna(ifos, cfg.ra_deg, cfg.dec_deg, 0.5 * (start + end))
    P = S.coherent_projector(sig, [ant[i][0] for i in ifos], [ant[i][1] for i in ifos])
    Z = np.stack(z)
    coh2 = np.real(np.einsum("ic,ij,jc->c", Z.conj(), P, Z))
    coh = np.sqrt(np.clip(coh2, 0, None))
    null = np.sqrt(np.clip(np.sum(np.abs(Z) ** 2, axis=0) - coh2, 0, None))
    t = start + np.arange(L) / fs
    ipk = int(np.argmax(coh))

    nb = max(int(round(0.01 * fs)), 1)                  # 10 ms bins, peak value

    def pool(a):
        n = len(a) // nb
        return a[:n * nb].reshape(n, nb).max(axis=1)

    k = int(round(args.zoom_s * fs))
    zs = slice(max(ipk - k, 0), ipk + k + 1)
    series = {"coherent": coh, "null": null, **{i: np.abs(Z[d]) for d, i in enumerate(ifos)}}
    out = {
        "mchirp": mc, "ifos": ifos, "gps_start": start, "gps_end": end, "sample_rate": fs,
        "gw170817_gps": GW170817_GPS, "gates": gates,
        "peak": {"gps": float(t[ipk]), "coherent": float(coh[ipk]), "null": float(null[ipk]),
                 **{i: float(abs(Z[d, ipk])) for d, i in enumerate(ifos)}},
        "search": {"reweighted": loud["loudest_rho_net"], "gps": loud["loudest_gps"]},
        "window": {"t": (start + (np.arange(L // nb) + 0.5) * nb / fs).round(4).tolist(),
                   **{n: np.round(pool(a), 3).tolist() for n, a in series.items()}},
        "zoom": {"t": t[zs].round(6).tolist(),
                 **{n: np.round(a[zs], 3).tolist() for n, a in series.items()}},
    }
    Path(args.out).write_text(json.dumps(out, default=float))
    print(json.dumps(out["peak"]), "reweighted in the search:", loud["loudest_rho_net"])


if __name__ == "__main__":
    main()
