"""Proposition 2 (power ceiling of the additive index): evaluates P0(I_A < c), where I_A is the
weighted sum of the seven count-based two-sided scores under the multinomial null at the pooled length.
Usage: python3 ceiling_bound.py   (writes results_extra/ceiling_bound.json)"""
import json, os
import numpy as np
from common import SEED, COUNT, SIGNED, WEIGHTS, build_null, build_seq_null, _mid

N = 15_688_140
null = build_null(N, np.random.default_rng(SEED))
seq = build_seq_null(np.random.default_rng(SEED + 3))
B = len(null["chi"]); rng = np.random.default_rng([SEED, 777])
T = {}
for k in COUNT:
    h = _mid(null[SIGNED.get(k, k)]); T[k] = np.minimum(1, 2 * np.minimum(h, 1 - h))
pick = rng.integers(0, len(seq["r1"]), size=B)
T["r1"] = _mid(seq["r1"])[pick]
h = _mid(seq["chiser"])[pick]; T["chiser"] = np.minimum(1, 2 * np.minimum(h, 1 - h))
IA = sum(WEIGHTS[k] * T[k] for k in COUNT)
I = IA + WEIGHTS["r1"] * T["r1"] + WEIGHTS["chiser"] * T["chiser"]
c = float(np.quantile(I, 0.05))
p = float((IA < 0.287).mean())
out = {"B": int(B), "weight_count_block": sum(WEIGHTS[k] for k in COUNT), "threshold_5pct_this_run": c,
       "threshold_used": 0.287, "P0_IA_below_threshold": p, "mc_se": float(np.sqrt(p * (1 - p) / B))}
os.makedirs("results_extra", exist_ok=True)
json.dump(out, open("results_extra/ceiling_bound.json", "w"), indent=1)
print(json.dumps(out, indent=1))
