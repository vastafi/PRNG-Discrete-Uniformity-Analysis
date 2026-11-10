#!/usr/bin/env python3
"""
make_figures.py

Draws Figure 1 of the paper (per-run composite index of every source, with the
null band) from the results.json written by reproduce_all.py. It does no
statistics of its own: every number plotted is read from that file.
Output: figures/figure1_per_run.pdf and a 600 dpi PNG.

Usage
-----
  python3 scripts/reproduce_all.py --data data --out results/main
  python3 scripts/make_figures.py  --results results/main/results.json --out figures
"""

import argparse, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DPI = 600
TIER_COLOUR = {1: "#1b7837", 2: "#5e81ac", 3: "#b2182b", 4: "#7b3294", 5: "#555555"}
TIER_HATCH = {1: "", 2: "//", 3: "xx", 4: "..", 5: "\\\\"}

WRAP = {"Python secrets": "Python\nsecrets", "Python numpy": "Python\nNumPy",
        "Java ThreadLocalRandom": "Java\nThreadLocalRandom",
        "Java SecureRandom": "Java\nSecureRandom", "pi": r"$\pi$"}

CRIT_LABEL = {"dmean": "mean", "dvar": "variance", "dskew": "skewness",
              "dexc": "exc. kurtosis", "dH": "entropy", "chi": r"$\chi^2$ GoF",
              "dks": "KS distance", "chiser": r"$\chi^2$ serial", "r1": "lag-1 autocorr."}


def label(name):
    return WRAP.get(name, name)


def style():
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.bbox": "tight"})


def save(fig, out, stem, pdf=True):
    fig.savefig(os.path.join(out, f"{stem}.png"), dpi=DPI)
    if pdf:
        fig.savefig(os.path.join(out, f"{stem}.pdf"))
    plt.close(fig)
    print("  wrote", stem)




def figure_per_run(res, out, pdf):
    per_run, cis = res["per_run_crqi"], res["per_run_ci"]
    idx = {k: v["index_twoside"] for k, v in res["scores"].items()}
    tiers, (lo, hi) = res["tiers"], res["null_band_95"]
    order = sorted(idx, key=lambda k: -idx[k])

    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.axhspan(lo, hi, color="0.92", zorder=0)
    for i, s in enumerate(order):
        v = np.asarray(per_run[s], float)
        c = TIER_COLOUR.get(tiers[s], "0.5")
        ax.plot([i, i], [v.min(), v.max()], color="0.6", lw=0.9, zorder=2)
        ax.scatter([i] * v.size, v, s=20, color=c, alpha=0.85, zorder=3)
        ax.plot([i - 0.3, i + 0.3], [v.mean()] * 2, color="k", lw=1.5, zorder=4)
        ax.fill_between([i - 0.18, i + 0.18], cis[s][0], cis[s][1],
                        color=c, alpha=0.18, zorder=1)
    ax.set_xticks(range(len(order)), [label(s) for s in order], rotation=45,
                  ha="right", fontsize=8.5)
    ax.set_ylabel("per-run composite index", fontsize=9)
    ax.tick_params(axis="y", labelsize=8.5)
    ax.set_ylim(0, 1)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Line2D([], [], ls="", marker="o", ms=5, color=c, label="run"),
               Line2D([], [], color="k", lw=1.5, label="mean of six runs"),
               Patch(color=c, alpha=0.25, label="95% CI of the mean"),
               Line2D([], [], color="0.6", lw=0.9, label="range of runs"),
               Patch(color="0.92", label="95% null band")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.22), ncol=3,
              frameon=False, fontsize=9.5, handlelength=1.4, columnspacing=1.0)
    save(fig, out, "figure1_per_run", pdf)






def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results/main/results.json")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--no-pdf", action="store_true")
    args = ap.parse_args()
    with open(args.results) as fh:
        res = json.load(fh)
    os.makedirs(args.out, exist_ok=True)
    style()
    pdf = not args.no_pdf
    figure_per_run(res, args.out, pdf)
    print(f"\nfigures written to {args.out}")


if __name__ == "__main__":
    main()