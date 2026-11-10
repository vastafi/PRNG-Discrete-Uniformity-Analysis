"""Which fluctuations does each tail convention reward?

Ten ideal sources are drawn at random from the null simulation, 20,000 times, and
ordered by the one-sided and by the two-sided index. For the source placed first
under each convention we record where its goodness-of-fit statistic X^2 falls in
its own null distribution. An ideal source picked at random falls in the lower 5%
tail (under-dispersion) with probability 0.05; any excess for the winner is what
the convention rewards.
Usage: winner_tail.py [--out results/extra/winner_tail.json]"""
import argparse, json
import numpy as np
from common import SEED, COUNT, SIGNED, WEIGHTS, build_null, build_seq_null, _mid

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="results/extra/winner_tail.json")
a = ap.parse_args()

N = 15_688_140
null = build_null(N, np.random.default_rng(SEED))
seq = build_seq_null(np.random.default_rng(SEED + 3))
B = len(null["chi"])
rng = np.random.default_rng([SEED, 5151])
U, T = {}, {}
for k in COUNT:
    U[k] = _mid(null[k])
    h = _mid(null[SIGNED.get(k, k)]); T[k] = np.minimum(1, 2 * np.minimum(h, 1 - h))
pick = rng.integers(0, len(seq["r1"]), size=B)
U["r1"] = T["r1"] = _mid(seq["r1"])[pick]
U["chiser"] = _mid(seq["chiser"])[pick]
h = _mid(seq["chiser"])[pick]; T["chiser"] = np.minimum(1, 2 * np.minimum(h, 1 - h))
idx = {"one_sided": sum(WEIGHTS[k] * U[k] for k in WEIGHTS),
       "two_sided": sum(WEIGHTS[k] * T[k] for k in WEIGHTS)}
left = 1.0 - _mid(null["chi"])          # left-tail probability of X^2 in its null

R, m = 20000, 10
res = {"R": R, "m": m}
for name, ix in idx.items():
    lt = np.empty(R)
    for i in range(R):
        s = rng.choice(B, m, replace=False)
        lt[i] = left[s[np.argmax(ix[s])]]
    res[name] = {"P_winner_X2_lower_5pct": float((lt < 0.05).mean()),
                 "P_winner_X2_lower_1pct": float((lt < 0.01).mean()),
                 "P_winner_X2_upper_5pct": float((lt > 0.95).mean()),
                 "median_left_tail_prob_of_winner": float(np.median(lt))}
json.dump(res, open(a.out, "w"), indent=1)
print(json.dumps(res, indent=1))
