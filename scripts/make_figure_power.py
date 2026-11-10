"""Power curves of six aggregation rules against two kinds of defect (paper figure).
Reads results/extra/pc_tilt.json and results/extra/pc_block.json (produced by power_curves.py).
Usage: python3 scripts/make_figure_power.py [--extra results/extra] [--out figures/figure2_power.pdf]"""
import argparse, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--extra", default="results/extra")
ap.add_argument("--out", default="figures/figure2_power.pdf")
a = ap.parse_args()

RULES = [("additive", "additive index", "#2a78d6", "o", "-"),
         ("geometric", "geometric", "#eb6834", "s", "--"),
         ("minimum", "minimum", "#1baf7a", "^", "-."),
         ("fisher", "Fisher", "#eda100", "D", ":"),
         ("X2_alone", "GoF alone", "#e87ba4", "v", "--"),
         ("serial_alone", "serial alone", "#008300", "P", "-.")]
PANELS = [("pc_tilt.json", "(a)", r"tilt of the marginal, $\delta$"),
          ("pc_block.json", "(b)", r"fraction of blocks reordered, $f$")]

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
for ax, (fn, letter, xlab) in zip(axes, PANELS):
    d = json.load(open(os.path.join(a.extra, fn)))
    grid = d["grid"]
    for key, lab, col, mk, ls in RULES:
        y = [d["rates"][str(g)][key] for g in grid]
        ax.plot(grid, y, color=col, marker=mk, ms=5, lw=1.6, ls=ls, label=lab)
    ax.axhline(0.05, color="#888888", lw=0.8, ls=":")
    ax.set_xscale("log"); ax.set_xticks(grid); ax.set_xticklabels([f"{g:g}" for g in grid])
    ax.minorticks_off(); ax.set_ylim(-0.02, 1.03); ax.set_xlabel(xlab)
    ax.text(-0.02, 1.04, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="right", va="bottom")
    ax.grid(axis="y", color="#e6e6e6", lw=0.6)
axes[0].set_ylabel("detection rate (200 streams)")
axes[0].text(0.00105, 0.075, "5% level", fontsize=7, color="#666666")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=6, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=(0, 0.08, 1, 1))
os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
fig.savefig(a.out, bbox_inches="tight"); fig.savefig(a.out.replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
print("saved", a.out)
