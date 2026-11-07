"""Null distribution of the disagreement between the one-sided and the two-sided
orderings of ten ideal sources, and where the observed disagreement falls in it.
Usage: reversal_null.py [--results results/results.json] [--out results_extra/reversal_null.json]"""
import argparse, json
import numpy as np
from scipy.stats import spearmanr
from common import SEED, COUNT, SIGNED, WEIGHTS, build_null, build_seq_null, _mid

ap = argparse.ArgumentParser()
ap.add_argument("--results", default="results/results.json")
ap.add_argument("--out", default="results_extra/reversal_null.json")
a = ap.parse_args()

# observed: one-sided (version-2) index vs corrected two-sided index on the ten sources
d = json.load(open(a.results)); S = d["sources"]
one = np.array([d["scores"][s]["crqi_v1_oneside"] for s in S])
two = np.array([d["scores"][s]["index_twoside"] for s in S])
rho_obs = float(spearmanr(one, two).statistic)
rank_two_of_first = int((-two).argsort().argsort()[np.argmax(one)]) + 1

# null: per-draw one-sided and two-sided scores from mid-ranks of the null simulation
N = 15_688_140
null = build_null(N, np.random.default_rng(SEED))
seq = build_seq_null(np.random.default_rng(SEED + 3))
B = len(null["chi"])
rng = np.random.default_rng([SEED, 4242])
U, T = {}, {}
for k in COUNT:
    U[k] = _mid(null[k])                                      # one-sided: upper tail
    h = _mid(null[SIGNED.get(k, k)]); T[k] = np.minimum(1, 2 * np.minimum(h, 1 - h))
pick = rng.integers(0, len(seq["r1"]), size=B)
U["r1"] = T["r1"] = _mid(seq["r1"])[pick]
U["chiser"] = _mid(seq["chiser"])[pick]
h = _mid(seq["chiser"])[pick]; T["chiser"] = np.minimum(1, 2 * np.minimum(h, 1 - h))
idx_u = sum(WEIGHTS[k] * U[k] for k in WEIGHTS)
idx_t = sum(WEIGHTS[k] * T[k] for k in WEIGHTS)

R, m = 20000, len(S)
rho, rank_first = np.empty(R), np.empty(R, dtype=int)
for i in range(R):
    s = rng.choice(B, m, replace=False)
    rho[i] = spearmanr(idx_u[s], idx_t[s]).statistic
    rank_first[i] = (-idx_t[s]).argsort().argsort()[np.argmax(idx_u[s])] + 1
out = {"R": R, "m": m, "observed_rho": rho_obs, "observed_rank_of_onesided_first": rank_two_of_first,
       "rho_median": float(np.median(rho)), "rho_q05_q95": [float(x) for x in np.quantile(rho, [0.05, 0.95])],
       "P_rho_le_observed": float((rho <= rho_obs).mean()),
       "P_first_falls_to_observed_rank_or_lower": float((rank_first >= rank_two_of_first).mean()),
       "P_first_becomes_last": float((rank_first == m).mean())}
json.dump(out, open(a.out, "w"), indent=1)
print(json.dumps(out, indent=1))
