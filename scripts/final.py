"""Summaries: corrected correlation table, G-based index variant, positive
control with Wilson intervals, Monte Carlo check of the power formula."""
import argparse, json, math
import numpy as np
from scipy.stats import spearmanr, binomtest, ttest_ind, skew, norm
from common import (SEED, COUNT, SEQ, SIGNED, WEIGHTS, build_null, build_seq_null, scores, index_value,
                    two_sided, load_joint)

ap = argparse.ArgumentParser()
ap.add_argument("--results", default="results/main/results.json")
ap.add_argument("--multiseed", default="results/multiseed/multiseed.json")
ap.add_argument("--pc", nargs="+", required=True, help="pc1000.py outputs, in order")
ap.add_argument("--joint-run", nargs="+", required=True)
ap.add_argument("--joint-pooled", nargs="+", required=True, help="joint_null.py outputs at the pooled length, in order")
ap.add_argument("--out", default="results/extra/final.json")
args = ap.parse_args()

out = {}
N = 15_688_140
null = build_null(N, np.random.default_rng(SEED))
seq = build_seq_null(np.random.default_rng(SEED + 3))

# 1. correlation table: count block on all 120,000 multinomial draws; every entry
#    involving a pair criterion from complete sequences (joint null), |r1| as scored
cols = COUNT
M = np.eye(9)
names = COUNT + SEQ
C = spearmanr(np.column_stack([null[k] for k in COUNT])).statistic
M[:7, :7] = C
for lab, files in [("pooled", args.joint_pooled), ("run", args.joint_run)]:
    J = load_joint(files); n = len(J["chi"])
    J["r1"] = np.abs(J["r1"])
    S = spearmanr(np.column_stack([J[k] for k in names])).statistic
    out[f"corr_{lab}_n"] = n
    out[f"corr_{lab}_cross"] = {a: {b: float(S[names.index(a), names.index(b)]) for b in SEQ} for a in COUNT}
    out[f"corr_{lab}_r1_chiser"] = float(S[7, 8])
    out[f"corr_{lab}_G_X"] = float(spearmanr(J["G"], J["chi"]).statistic)
    if lab == "pooled":
        M[7:, :] = S[7:, :]; M[:, 7:] = S[:, 7:]
out["corr_table"] = {"names": names, "M": M.tolist()}
out["corr_dH_chi_full"] = float(C[COUNT.index("dH"), COUNT.index("chi")])

# 2. index with G in place of Pearson's statistic
d = json.load(open(args.results))
from common import stats_from_counts
diff = {}
for s, o in d["observed"].items():
    obs = {k: o[k] for k in COUNT + list(SIGNED.values())}; obs.update(r1=o["r1"], chiser=o["chiser"])
    t = scores(null, seq, obs, N); base = index_value(t)
    tG = dict(t); tG["chi"] = two_sided(float(stats_from_counts(o["counts"])["G"][0]), null["G"])
    diff[s] = {"pearson": base, "G": index_value(tG), "stored": d["scores"][s].get("index_twoside")}
out["index_G_variant"] = diff
out["index_G_max_abs_diff"] = max(abs(v["pearson"] - v["G"]) for v in diff.values())

# 3. positive control, 1000 replicates
parts = [json.load(open(p)) for p in args.pc]
A = parts[0]
assert all(abs(p["thr"] - A["thr"]) < 1e-12 for p in parts)
pc = {"thr": A["thr"], "thr_geo": A["thr_geo"]}
for s in ["clean", "A", "B", "C"]:
    add = np.array(sum((p["res"][s]["add"] for p in parts), []))
    geo = np.array(sum((p["res"][s]["geo"] for p in parts), []))
    R = add.size
    row = {"R": int(R), "mean_add": float(add.mean()), "sd_add": float(add.std(ddof=1)),
           "mean_geo": float(geo.mean())}
    for lab, v, th in [("add", add, A["thr"]), ("geo", geo, A["thr_geo"])]:
        k = int((v < th).sum()); ci = binomtest(k, R).proportion_ci(method="wilson")
        row[f"rate_{lab}"] = k / R; row[f"ci_{lab}"] = [ci.low, ci.high]
    pc[s] = row
out["positive_control"] = pc

# 4. Monte Carlo check of the normal-approximation sample size, on the
#    multi-seed index values (centred by generator, bounded on [0,1])
ms = json.load(open(args.multiseed))
vals = [np.asarray(v, float) for v in ms["index_by_seed"].values()]
resid = np.concatenate([v - v.mean() for v in vals]) + 0.5
out["multiseed_shape"] = {"min": float(resid.min()), "max": float(resid.max()),
                          "skew": float(skew(resid)), "sd": float(resid.std(ddof=1))}
rng = np.random.default_rng(SEED + 11)
chk = {}
for n_runs, alpha in [(26, 0.05), (56, 0.05 / 45), (21, 0.05), (45, 0.05 / 45)]:
    hits = 0; R = 4000
    for _ in range(R):
        a = rng.choice(resid, n_runs); b = np.clip(rng.choice(resid, n_runs) + 0.1, 0, 1)
        hits += ttest_ind(a, b).pvalue < alpha
    chk[f"n{n_runs}_a{alpha:.5f}"] = hits / R
out["power_check_delta0.1"] = chk
json.dump(out, open(args.out, "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "corr_table"}, indent=1))
print(np.round(M, 3))
