"""Quick analyses from stored results: G statistic, Pearson-G agreement,
ANOVA without pi, sensitivity of the required number of runs."""
import argparse, json, math
import numpy as np
from scipy.stats import chi2, norm, f as fdist, spearmanr
from common import stats_from_counts, build_null, SEED, HX

ap = argparse.ArgumentParser()
ap.add_argument("--results", default="results/results.json")
ap.add_argument("--out", default="results_extra/quick.json")
args = ap.parse_args()
d = json.load(open(args.results))
N = 15_688_140
out = {}

# --- G for every source, pooled and per run
rows = {}
for s, o in d["observed"].items():
    st = stats_from_counts(o["counts"])
    G, X = float(st["G"][0]), float(st["chi"][0])
    rows[s] = {"X2": X, "G": G, "pG": float(chi2.sf(G, 9)), "pX": float(chi2.sf(X, 9)),
               "G_from_entropy": float(2 * N * math.log(2) * (HX - (HX - st["dH"][0]))),
               "rel_diff": (G - X) / X}
out["G_sources"] = rows
order_G = sorted(rows, key=lambda s: rows[s]["G"])
order_X = sorted(rows, key=lambda s: rows[s]["X2"])
out["same_order_G_X"] = order_G == order_X
out["max_abs_rel_diff_sources"] = max(abs(r["rel_diff"]) for r in rows.values())

# --- Pearson vs G over the multinomial null at the pooled N and at the per-run N
for lab, n in [("pooled", N), ("run", 2_614_690)]:
    nl = build_null(n, np.random.default_rng(SEED))
    rel = np.abs(nl["G"] / nl["chi"] - 1)
    ident = np.max(np.abs((2 * n * math.log(2)) * nl["dH"] / nl["G"] - 1))
    out[f"null_{lab}"] = {"spearman_G_X": float(spearmanr(nl["G"], nl["chi"]).statistic),
                          "max_rel_GX": float(rel.max()), "median_rel_GX": float(np.median(rel)),
                          "max_rel_identity_dH_G": float(ident),
                          "frac_rank_swaps_pairs": None}

# --- ANOVA on per-run indices, with and without pi (and without Random.org)
def anova(groups):
    k = len(groups); n_i = np.array([len(g) for g in groups]); allv = np.concatenate(groups)
    grand = allv.mean()
    msb = np.sum(n_i * (np.array([np.mean(g) for g in groups]) - grand) ** 2) / (k - 1)
    msw = np.sum([np.sum((np.asarray(g) - np.mean(g)) ** 2) for g in groups]) / (n_i.sum() - k)
    F = msb / msw
    return {"F": float(F), "df1": int(k - 1), "df2": int(n_i.sum() - k),
            "p": float(fdist.sf(F, k - 1, n_i.sum() - k)), "sigma_w": float(math.sqrt(msw)),
            "range_means": float(np.ptp([np.mean(g) for g in groups])),
            "mean_within_range": float(np.mean([np.ptp(g) for g in groups]))}

pr = d["per_run_crqi"]
out["anova_all10"] = anova([pr[s] for s in pr])
out["anova_no_pi"] = anova([pr[s] for s in pr if s != "pi"])
out["anova_prng8"] = anova([pr[s] for s in pr if s not in ("pi", "Random.org")])
out["pi_window_sd"] = float(np.std(pr["pi"], ddof=1))
out["pi_window_range"] = float(np.ptp(pr["pi"]))

# --- required runs, with a 95% CI on sigma and a sensitivity grid
def n_req(sig, delta, alpha):
    return math.ceil(2 * (norm.ppf(1 - alpha / 2) + norm.ppf(0.8)) ** 2 * sig ** 2 / delta ** 2)

sens = {}
for lab, sig, df in [("anova_no_pi", out["anova_no_pi"]["sigma_w"], out["anova_no_pi"]["df2"]),
                     ("multiseed", 0.129, 1194)]:
    lo = sig * math.sqrt(df / chi2.ppf(0.975, df)); hi = sig * math.sqrt(df / chi2.ppf(0.025, df))
    sens[lab] = {"sigma": sig, "ci": [lo, hi],
                 "n_single": [n_req(x, 0.1, 0.05) for x in (lo, sig, hi)],
                 "n_bonf": [n_req(x, 0.1, 0.05 / 45) for x in (lo, sig, hi)]}
out["sigma_sensitivity"] = sens
sg = out["anova_no_pi"]["sigma_w"]
out["power_table"] = {str(D): [n_req(sg, D, 0.05), n_req(sg, D, 0.05 / 45)]
                      for D in (0.05, 0.10, 0.15, 0.20, 0.25)}
json.dump(out, open(args.out, "w"), indent=1)
print(json.dumps(out, indent=1))
