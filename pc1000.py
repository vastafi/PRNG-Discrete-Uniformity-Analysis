"""Positive control with many replicates . Usage: pc1000.py --reps R --start S --out f.json"""
import argparse, json, time
import numpy as np
from common import (SEED, build_null, build_seq_null, scores, index_value, index_geom,
                    null_index, tilt, block_sort, criteria_of)

ap = argparse.ArgumentParser()
ap.add_argument("--reps", type=int, default=1000); ap.add_argument("--start", type=int, default=0)
ap.add_argument("--N", type=int, default=15_688_140); ap.add_argument("--out", required=True)
a = ap.parse_args(); N = a.N
null = build_null(N, np.random.default_rng(SEED))
seq = build_seq_null(np.random.default_rng(SEED + 3))
thr = float(np.quantile(null_index(null, seq, np.random.default_rng(SEED + 2)), 0.05))
thr_g = float(np.quantile(null_index(null, seq, np.random.default_rng(SEED + 2), form="geometric"), 0.05))
print("thresholds", thr, thr_g, flush=True)
names = ["clean", "A", "B", "C"]
res = {n: {"add": [], "geo": [], "t": []} for n in names}
t0 = time.time()
for rep in range(a.start, a.start + a.reps):
    rng = np.random.default_rng([SEED, 777, rep])
    ref = rng.integers(0, 10, N, dtype=np.int8)
    st = {"clean": ref, "A": tilt(ref, 0.01, rng), "B": block_sort(ref, 4, 1.0, rng),
          "C": block_sort(tilt(ref, 0.001, rng), 4, 0.1, rng)}
    for n in names:
        t = scores(null, seq, criteria_of(st[n]), N)
        res[n]["add"].append(index_value(t)); res[n]["geo"].append(index_geom(t)); res[n]["t"].append(t)
    if (rep - a.start) % 10 == 0:
        print(rep, round(time.time() - t0), flush=True)
        json.dump({"thr": thr, "thr_geo": thr_g, "res": res}, open(a.out, "w"))
json.dump({"thr": thr, "thr_geo": thr_g, "res": res}, open(a.out, "w"))
print("done", flush=True)
