"""Joint null of all nine criteria (plus G) from complete uniform sequences.
Usage: joint_null.py --N n --B b --out f.npz"""
import argparse, time
import numpy as np
from common import SEED, criteria_of

ap = argparse.ArgumentParser()
ap.add_argument("--N", type=int); ap.add_argument("--B", type=int); ap.add_argument("--out"); ap.add_argument("--tag", type=int, default=0)
a = ap.parse_args()
rng = np.random.default_rng([SEED, 999, a.N] + ([a.tag] if a.tag else []))
keys = ["dmean", "dvar", "dskew", "dexc", "smean", "svar", "sskew", "sexc",
        "dH", "chi", "G", "dks", "r1", "chiser"]
out = {k: np.empty(a.B) for k in keys}
t0 = time.time()
for b in range(a.B):
    o = criteria_of(rng.integers(0, 10, a.N, dtype=np.int8))
    for k in keys:
        out[k][b] = o[k]
    if b % 100 == 0:
        print(b, round(time.time() - t0), flush=True)
        np.savez(a.out, n_done=b + 1, **{k: v[:b + 1] for k, v in out.items()})
np.savez(a.out, n_done=a.B, **out)
print("done", flush=True)
