#!/usr/bin/env python3
"""
make_figures.py  (version 2)

Draws Figures 1 to 4 from the results.json produced by reproduce_all.py v2.
It does no statistics of its own: every number plotted is read from that file.

Changes with respect to version 1
---------------------------------
* No hard-coded TOP_CUT.  Version 1 decided the "top tier" with a constant
  0.60 that appears nowhere in the paper; the tiers are now read from
  results.json, where they come from the overlap of the per-run confidence
  intervals.
* Figure 2 now shows the null band and the per-run confidence interval, which
  is what makes the within-source versus between-source comparison legible.
* Two new figures:
    Figure 3 - the score of every source on every criterion, which shows at a
               glance which criterion drives each result;
    Figure 4 - the Spearman correlation of the nine criteria under H0, which
               documents the redundancy of the criterion set, in particular
               the perfect correlation between the entropy deficit and the
               goodness-of-fit chi-square.
* Vector PDF by default plus 600 dpi PNG, and a colour scheme that stays
  legible in greyscale; tier is encoded by colour AND by hatching.

Usage
-----
  python3 reproduce_all.py --data data --out results
  python3 make_figures.py  --results results/results.json --out figures
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


def figure1(res, out, pdf):
    idx = {k: v["index_twoside"] for k, v in res["scores"].items()}
    tiers = res["tiers"]
    thr, (lo, hi) = res["threshold_5pct"], res["null_band_95"]
    order = sorted(idx, key=lambda k: -idx[k])

    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    ax.axhspan(lo, hi, color="0.90", zorder=0)
    for i, s in enumerate(order):
        t = tiers[s]
        ax.bar(i, idx[s], width=0.66, color=TIER_COLOUR.get(t, "0.5"),
               hatch=TIER_HATCH.get(t, ""), edgecolor="white", linewidth=0.6, zorder=3)
    ax.axhline(thr, ls="--", c="k", lw=1.0, zorder=4)
    ax.axhline(hi, ls=":", c="0.35", lw=1.0, zorder=4)
    ax.set_xticks(range(len(order)), [label(s) for s in order], rotation=45,
                  ha="right", fontsize=7.5)
    ax.set_ylabel("composite index (two-sided calibration)")
    ax.set_ylim(0, 1)
    sidak = res.get("threshold_sidak")
    extra_h, extra_l = [], []
    if sidak is not None:
        ax.axhline(sidak, ls="-.", c="0.25", lw=1.0, zorder=4)
        extra_h = [plt.Line2D([], [], ls="-.", c="0.25")]
        extra_l = [f"Šidák-corrected alarm ({sidak:.3f})"]
    # a single tier is not a classification, so it is not shown in the legend
    ax.legend([plt.Line2D([], [], ls="--", c="k")] + extra_h
              + [plt.Line2D([], [], ls=":", c="0.35"),
                 plt.Rectangle((0, 0), 1, 1, color="0.90")],
              [f"5% alarm ({thr:.3f})"] + extra_l + ["97.5% null bound", "95% null band"],
              fontsize=6.5, ncol=2, loc="upper right", frameon=False)
    save(fig, out, "figure1", pdf)


def figure2(res, out, pdf):
    per_run, cis = res["per_run_crqi"], res["per_run_ci"]
    idx = {k: v["index_twoside"] for k, v in res["scores"].items()}
    tiers, (lo, hi) = res["tiers"], res["null_band_95"]
    order = sorted(idx, key=lambda k: -idx[k])

    fig, ax = plt.subplots(figsize=(6.6, 3.6))
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
                  ha="right", fontsize=7.5)
    ax.set_ylabel("per-run composite index")
    ax.set_ylim(0, 1)
    w = res.get("within_source_range_mean")
    b = res.get("between_source_range")
    if w and b:
        ax.set_title(f"mean within-source range {w:.2f} vs between-source range {b:.2f}",
                     fontsize=8)
    save(fig, out, "figure2", pdf)


def figure3(res, out, pdf):
    crits = res["diagnostics"]["criterion_correlations"]["criteria"]
    idx = {k: v["index_twoside"] for k, v in res["scores"].items()}
    order = sorted(idx, key=lambda k: -idx[k])
    M = np.array([[res["scores"][s]["t_twoside"][c] for c in crits] for s in order])

    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    im = ax.imshow(M, cmap="RdYlGn", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(crits)), [CRIT_LABEL.get(c, c) for c in crits],
                  rotation=40, ha="right")
    ax.set_yticks(range(len(order)), [s.replace("\n", " ") for s in order], fontsize=7.5)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=5.5,
                    color="0.1")
    fig.colorbar(im, ax=ax, shrink=0.9, label="two-sided typicality score")
    save(fig, out, "figure3", pdf)


def figure4(res, out, pdf):
    cc = res["diagnostics"]["criterion_correlations"]
    crits, M = cc["criteria"], np.array(cc["spearman"])
    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1)
    lab = [CRIT_LABEL.get(c, c) for c in crits]
    ax.set_xticks(range(len(crits)), lab, rotation=40, ha="right")
    ax.set_yticks(range(len(crits)), lab)
    for i in range(len(crits)):
        for j in range(len(crits)):
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=5.5,
                    color="0.1")
    fig.colorbar(im, ax=ax, shrink=0.85, label=r"Spearman correlation under $H_0$")
    save(fig, out, "figure4", pdf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results/results.json")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--no-pdf", action="store_true")
    args = ap.parse_args()
    with open(args.results) as fh:
        res = json.load(fh)
    os.makedirs(args.out, exist_ok=True)
    style()
    pdf = not args.no_pdf
    figure1(res, args.out, pdf)
    figure2(res, args.out, pdf)
    figure3(res, args.out, pdf)
    figure4(res, args.out, pdf)
    print(f"\nfigures written to {args.out}")


if __name__ == "__main__":
    main()