"""Detection rate vs defect size for six detection rules.
Usage: power_curves.py --defect tilt|block --reps R --out f.json"""
import argparse, json, time
import numpy as np
from scipy.stats import chi2
from common import (SEED, COUNT, SEQ, SIGNED, WEIGHTS, build_null, build_seq_null, scores,
                    index_value, index_geom, _mid, tilt, block_sort, criteria_of)

ap = argparse.ArgumentParser()
ap.add_argument("--defect", choices=["tilt", "block"], required=True)
ap.add_argument("--reps", type=int, default=200); ap.add_argument("--out", required=True)
a = ap.parse_args(); N = 15_688_140
null = build_null(N, np.random.default_rng(SEED))
seq = build_seq_null(np.random.default_rng(SEED + 3))
B = len(null["chi"])
T = {}
for k in COUNT:
    h = _mid(null[SIGNED.get(k, k)]); T[k] = np.minimum(1, 2 * np.minimum(h, 1 - h))
pick = np.random.default_rng(SEED + 2).integers(0, len(seq["r1"]), size=B)
T["r1"] = _mid(seq["r1"])[pick]
h = _mid(seq["chiser"])[pick]; T["chiser"] = np.minimum(1, 2 * np.minimum(h, 1 - h))
K = list(WEIGHTS)
M = np.column_stack([T[k] for k in K])
fl = 1.0 / B
stat_null = {"additive": M @ np.array([WEIGHTS[k] for k in K]),
             "geometric": np.exp(np.log(np.maximum(M, fl)) @ np.array([WEIGHTS[k] for k in K])),
             "minimum": M.min(axis=1),
             "fisher": -(-2 * np.log(np.maximum(M, fl)).sum(axis=1))}   # low = alarm
thr = {r: float(np.quantile(v, 0.05)) for r, v in stat_null.items()}
grid = [0.001, 0.002, 0.003, 0.004, 0.006, 0.008] if a.defect == "tilt" else [0.002, 0.005, 0.01, 0.02, 0.05, 0.1]
res = {"thr": thr, "grid": grid, "rates": {}}
t0 = time.time()
for g in grid:
    hits = {r: 0 for r in list(thr) + ["X2_alone", "serial_alone"]}
    for rep in range(a.reps):
        rng = np.random.default_rng([SEED, 31, int(g * 1e6), rep, 0 if a.defect == "tilt" else 1])
        x = rng.integers(0, 10, N, dtype=np.int8)
        x = tilt(x, g, rng) if a.defect == "tilt" else block_sort(x, 4, g, rng)
        o = criteria_of(x); t = scores(null, seq, o, N)
        v = np.array([t[k] for k in K])
        st = {"additive": index_value(t), "geometric": index_geom(t), "minimum": v.min(),
              "fisher": -(-2 * np.log(np.maximum(v, fl)).sum())}
        for r in thr:
            hits[r] += st[r] < thr[r]
        hits["X2_alone"] += chi2.sf(o["chi"], 9) < 0.05
        hits["serial_alone"] += chi2.sf(o["chiser"], 81) < 0.05
    res["rates"][str(g)] = {r: hits[r] / a.reps for r in hits}
    print(g, res["rates"][str(g)], round(time.time() - t0), flush=True)
    json.dump(res, open(a.out, "w"), indent=1)
print("done", flush=True)
