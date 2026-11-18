"""Public GW170817 sky maps, reduced to small multi-order files for the logbook fixture.

Source: LIGO DCC G1701985 (public), https://dcc.ligo.org/public/0146/G1701985/001/
  bayestar.fits.gz          rapid HLV localisation, circulated on 2017-08-17 (GCN 21513)
  LALInference_v2.fits.gz   localisation published with the discovery (2017-10-16)

Each map is brought to nside 256 (probability summed, distance layers
probability-weighted); pixels outside the 99.99 % region are merged to
nside 16. Writes GW170817_<name>.multiorder.fits beside this script.

  python examples/logbook_fixture/fetch_gw170817_maps.py
"""
import sys
import urllib.request
from pathlib import Path

import healpy as hp
import numpy as np
from astropy.table import Table

HERE = Path(__file__).resolve().parent
URL = "https://dcc.ligo.org/public/0146/G1701985/001/"
FINE, COARSE = 8, 4


def reduce(path: Path) -> Table:
    prob, mu, sigma, norm = hp.read_map(path, field=(0, 1, 2, 3), nest=True)
    order = hp.nside2order(hp.get_nside(prob))
    k = 4 ** (order - FINE)

    def down(layer, kk):                              # probability-weighted mean of the children
        w = prob.reshape(-1, kk) if kk == k else None
        v = np.where(np.isfinite(layer), layer, 0.0).reshape(-1, kk)
        ws = w.sum(axis=1)
        return np.where(ws > 0, (v * w).sum(axis=1) / np.where(ws > 0, ws, 1), np.inf)
    p8 = prob.reshape(-1, k).sum(axis=1)
    layers8 = [down(x, k) for x in (mu, sigma, norm)]
    srt = np.argsort(-p8)
    cred = np.empty_like(p8)
    cred[srt] = np.cumsum(p8[srt])
    keep = cred <= 0.9999
    kc = 4 ** (FINE - COARSE)
    parent_keep = keep.reshape(-1, kc).any(axis=1)     # a coarse pixel is split if any child is kept
    uniq, dens, cols = [], [], [[], [], []]
    a8, a4 = hp.nside2pixarea(1 << FINE), hp.nside2pixarea(1 << COARSE)
    for ip in range(12 * 4 ** COARSE):
        ch = slice(ip * kc, (ip + 1) * kc)
        if parent_keep[ip]:
            uniq.extend(4 * 4 ** FINE + np.arange(ip * kc, (ip + 1) * kc))
            dens.extend(p8[ch] / a8)
            for c, l in zip(cols, layers8):
                c.extend(l[ch])
        else:
            uniq.append(4 * 4 ** COARSE + ip)
            w = p8[ch]
            dens.append(w.sum() / a4)
            for c, l in zip(cols, layers8):
                good = np.isfinite(l[ch]) & (w > 0)
                c.append(float((l[ch][good] * w[good]).sum() / w[good].sum()) if good.any() else np.inf)
    t = Table({"UNIQ": np.array(uniq, dtype=np.int64), "PROBDENSITY": np.array(dens),
               "DISTMU": np.array(cols[0]), "DISTSIGMA": np.array(cols[1]), "DISTNORM": np.array(cols[2])})
    t.meta["ORIGIN"] = "LIGO DCC G1701985, " + path.name
    return t


def main() -> None:
    for name in ("bayestar", "LALInference_v2"):
        raw = HERE / f"_{name}.fits.gz"
        if not raw.exists():
            urllib.request.urlretrieve(URL + f"{name}.fits.gz", raw)
        t = reduce(raw)
        out = HERE / f"GW170817_{name}.multiorder.fits"
        t.write(out, overwrite=True)
        area = 4 * 180 ** 2 / np.pi
        print(f"{out.name}: {len(t)} pixels, {out.stat().st_size / 1e3:.0f} kB, "
              f"total probability {float(np.sum(t['PROBDENSITY'] * np.pi / (3 * 4.0 ** (np.floor(np.log2(t['UNIQ'] / 4) / 2))))):.6f}")
        raw.unlink()


if __name__ == "__main__":
    sys.exit(main())
