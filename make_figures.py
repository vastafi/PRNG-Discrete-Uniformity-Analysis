#!/usr/bin/env python3
"""
Draws Figure 1 and Figure 2 of

  "Concordance with the Discrete Uniform Distribution of Pseudorandom Number
   Generators in the Most Commonly Used Computer Applications"
  V. Astafi, A. Leahu, D. Ciorba

from the results.json produced by reproduce_all.py. It does no statistics of
its own: every number it plots is read from that file, so the figures always
agree with the tables.

Usage
-----
  python3 reproduce_all.py --data data --out results      # first
  python3 make_figures.py  --results results/results.json --out figures

Output: figures/figure1.png and figures/figure2.png at 300 dpi, the resolution
required for submission. Add --pdf for vector versions as well.

Requirements: numpy, matplotlib.
Runtime: a couple of seconds.
"""

import argparse, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# Appearance
# ----------------------------------------------------------------------------

TOP_CUT = 0.60        # a source is "top tier" at or above this index value;
                      # "bottom tier" below the 5% alarm threshold; middle otherwise
DPI     = 300
FIGSIZE = (8.2, 4.3)

COLOURS = {"top": "#2e7d32", "middle": "#7b83c4", "bottom": "#b3282d"}

# two-line labels keep the x axis readable
WRAP = {
    "Python secrets":         "Python\nsecrets",
    "Python numpy":           "Python\nnumpy",
    "Java ThreadLocalRandom": "Java\nThreadLocalRandom",
    "Java SecureRandom":      "Java\nSecureRandom",
    "pi":                     r"$\pi$",
}

def label(name):
    return WRAP.get(name, name)

def tier_of(value, threshold):
    if value >= TOP_CUT:
        return "top"
    if value < threshold:
        return "bottom"
    return "middle"

# ----------------------------------------------------------------------------
# Figures
# ----------------------------------------------------------------------------

def figure1(res, path):
    """Composite index per source, coloured by tier, with the alarm threshold
    and the 95% null band."""
    scores = {k: v["crqi"] for k, v in res["scores"].items()}
    thr    = res["threshold_5pct"]
    lo, hi = res["null_band_95"]
    order  = sorted(scores, key=lambda k: -scores[k])
    values = [scores[k] for k in order]
    tiers  = [tier_of(v, thr) for v in values]

    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    ax.axhspan(lo, hi, color="0.90", zorder=0)
    ax.bar(range(len(order)), values,
           color=[COLOURS[t] for t in tiers], zorder=3, width=0.68)
    ax.axhline(thr, ls="--", c="k", lw=1.1, zorder=4)

    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([label(k) for k in order], rotation=45, ha="right", fontsize=7.5)
    ax.set_ylabel("Composite index (CRQI)")
    ax.set_ylim(0, 1.0)

    handles = [plt.Rectangle((0, 0), 1, 1, color=COLOURS[t]) for t in ("top", "middle", "bottom")]
    ax.legend(handles + [plt.Line2D([], [], ls="--", c="k"),
                         plt.Rectangle((0, 0), 1, 1, color="0.90")],
              ["top tier", "middle tier", "bottom tier",
               f"5% alarm ({thr:.3f})", "95% null band"],
              fontsize=7, loc="upper right", ncol=2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path

def figure2(res, path):
    """Per-run composite index: one dot per run, a bar at the source mean and a
    vertical line spanning the within-source range."""
    scores  = {k: v["crqi"] for k, v in res["scores"].items()}
    per_run = res["per_run_crqi"]
    thr     = res["threshold_5pct"]
    order   = sorted(scores, key=lambda k: -scores[k])

    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    for i, name in enumerate(order):
        v = per_run[name]
        colour = COLOURS[tier_of(scores[name], thr)]
        ax.scatter([i] * len(v), v, s=22, color=colour, zorder=3, alpha=0.85)
        ax.plot([i - 0.28, i + 0.28], [np.mean(v)] * 2, color="k", lw=1.4, zorder=4)
        ax.plot([i, i], [min(v), max(v)], color="0.55", lw=0.9, zorder=2)

    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([label(k) for k in order], rotation=45, ha="right", fontsize=7.5)
    ax.set_ylabel("Per-run composite index")
    ax.legend([plt.Line2D([], [], color="k", lw=1.4)], ["source mean"], fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path

# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results/results.json",
                    help="the file written by reproduce_all.py")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--pdf", action="store_true", help="also write vector PDF versions")
    args = ap.parse_args()

    with open(args.results) as fh:
        res = json.load(fh)
    os.makedirs(args.out, exist_ok=True)

    print(f"{len(res['scores'])} source(s); alarm threshold {res['threshold_5pct']:.3f}")
    for ext in (["png", "pdf"] if args.pdf else ["png"]):
        p1 = figure1(res, os.path.join(args.out, f"figure1.{ext}"))
        p2 = figure2(res, os.path.join(args.out, f"figure2.{ext}"))
        print("  wrote", p1)
        print("  wrote", p2)

    within  = res.get("within_source_range_mean")
    between = res.get("between_source_range")
    if within is not None and not np.isnan(within):
        print(f"\nFigure 2 caption figures: mean within-source range {within:.2f}, "
              f"between-source range of means {between:.2f}")


if __name__ == "__main__":
    main()