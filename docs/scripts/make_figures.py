"""Figures for the documentation site (docs/img/), from the committed example results.

Run from the repository root:  .venv/bin/python docs/scripts/make_figures.py
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "img"

# dark chart surface, matching the workflow diagram
SURFACE = "#161b22"
TEXT = "#e6edf3"
MUTED = "#9da7b3"
GRID = "#30363d"
SERIES = ["#3987e5", "#d95926", "#199e70"]   # validated categorical slots 1-3 (dark)

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
    "text.color": TEXT, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 11,
    "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "legend.frameon": False, "lines.linewidth": 2, "axes.axisbelow": True,
})

SURVEYS = [("runs_inj", "ZTF-like"), ("runs_inj_rubin", "Rubin-like"), ("runs_inj_mixed", "ZTF + Rubin")]


def _rows(run):
    rows = json.loads((ROOT / "examples" / "injections" / run / "report.json").read_text())
    return [r for r in rows if r["fit_model"] in ("Bu2019lm", "Ka2017")]


def margin_coverage():
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    x = np.linspace(0, 1.2, 241)
    for (run, label), c in zip(SURVEYS, SERIES):
        f = np.array([r["f_needed"] for r in _rows(run)])
        y = [(f <= xi).mean() * 100 for xi in x]
        ax.plot(x, y, color=c, label=f"{label} ({len(f)} fits)")
    ax.axvline(0.75, color=MUTED, ls="--", lw=1)
    ax.text(0.765, 72, "default 0.75", color=MUTED, fontsize=10)
    ax.set_xlabel("model-systematic margin (fraction of the extrapolation gap)")
    ax.set_ylabel("fits whose window contains the merger (%)")
    ax.set_ylim(60, 101)
    ax.set_xlim(0, 1.2)
    ax.set_title("Merger inside the t0 window vs. margin, synthetic kilonovae")
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "injections_margin_coverage.png", dpi=150)
    plt.close(fig)


def window_widths():
    fig, ax = plt.subplots(figsize=(7.5, 3.9))
    labels = [l for _, l in SURVEYS]
    y = np.arange(len(labels))
    h = 0.36
    for k, (cov, name, c) in enumerate([(True, "field observed before the merger", SERIES[0]),
                                        (False, "late discovery", SERIES[1])]):
        med = []
        for run, _ in SURVEYS:
            w = [r["width_h_f0.75"] for r in _rows(run) if r["pre_merger_coverage"] == cov]
            med.append(np.median(w))
        pos = y + (k - 0.5) * (h + 0.04)
        ax.barh(pos, med, height=h, color=c, label=name)
        for p, m in zip(pos, med):
            ax.text(m + 2, p, f"{m:.0f} h", va="center", color=TEXT, fontsize=10)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("median t0 window (hours, margin 0.75)")
    ax.set_xlim(0, 170)
    ax.grid(axis="y", visible=False)
    ax.set_title("Width of the t0 window by survey", pad=28)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=10, borderaxespad=0.2)
    fig.tight_layout()
    fig.savefig(OUT / "injections_window_width.png", dpi=150)
    plt.close(fig)


def gw170817_background():
    d = ROOT / "examples" / "targeted_gw170817"
    fig, axes = plt.subplots(2, 1, figsize=(7.5, 4.8), sharex=True)
    for ax, stat, c in zip(axes, ["network", "coherent"], SERIES[:2]):
        r = json.loads((d / f"GW170817_eventframes_{stat}_targeted_result.json").read_text())
        s = r["subsegments"][0]
        bkg = np.asarray(s["background"])
        ax.hist(bkg, bins=np.arange(5, 33, 0.25), color=c, edgecolor=SURFACE, linewidth=0.5)
        rho = s["loudest_rho_net"]
        ax.axvline(rho, color=TEXT, lw=1.5)
        ax.text(rho - 0.4, ax.get_ylim()[1] * 0.75, f"GW170817  {rho:.1f}", ha="right", color=TEXT)
        ax.text(bkg.max() + 0.4, ax.get_ylim()[1] * 0.45,
                f"{len(bkg)} time slides\nmax {bkg.max():.1f}", color=MUTED, fontsize=10)
        ax.set_ylabel("slides")
        ax.set_title(f"{stat} statistic", fontsize=11)
    axes[-1].set_xlabel("loudest SNR in the window (newSNR, PyGRB gate)")
    fig.suptitle("GW170817: on-source loudest vs. time-slide background", x=0.02, ha="left",
                 fontweight="bold", color=TEXT)
    fig.tight_layout()
    fig.savefig(OUT / "gw170817_background.png", dpi=150)
    plt.close(fig)


def copy_fits():
    pairs = {
        "examples/validation/runs_fit3_Bu2019lm/A/nmma_out/AT2017gfo_A_bestfit_lightcurves.png":
            "at2017gfo_bu2019lm_bestfit.png",
        "examples/validation_grb211211a/runs_fit3_Bu2019lm/C/nmma_out/GRB211211A_C_bestfit_lightcurves.png":
            "grb211211a_bu2019lm_bestfit.png",
    }
    for src, dst in pairs.items():
        shutil.copy(ROOT / src, OUT / dst)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    margin_coverage()
    window_widths()
    gw170817_background()
    copy_fits()
    print("written to", OUT)
