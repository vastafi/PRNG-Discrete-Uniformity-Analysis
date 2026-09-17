#!/usr/bin/env python3
"""
multiseed.py

Turns a statement about one generated stream into a statement about a
generator.  For every algorithmic source the script draws S independent seeds,
scores each seeded stream with the nine criteria, and reports the distribution
of the composite index per generator together with a formal comparison.

This is the experiment that decides whether the ranking of the paper is a
property of the generators or of the particular runs that were analysed.  Two
numbers come out of it:

    sigma_within   spread of the index across seeds of the same generator
    sigma_between  spread of the per-generator means

and the F test of the "generator" effect.  If the effect is not significant,
that is the headline result of the paper and it is a much stronger result than
a ranking.

Generators covered natively: the numpy bit generators, Python's random module
and Python's secrets module.  Sources that live outside Python (R, Java, C++,
Excel, Mathematica) must be exported by their own environment; put them in a
folder with one file per seed and pass --external FOLDER.  The folder layout is
the same as for reproduce_all.py: one subfolder per source, one file per seed.

Usage
-----
    python3 multiseed.py --seeds 200 --n 1000000 --out results_multiseed
    python3 multiseed.py --seeds 200 --n 1000000 --external data_seeds
"""

import argparse, json, os
import numpy as np
from scipy.stats import f as fdist

from reproduce_all import (stats_from_counts, lag1, chi2_serial, build_null,
                           build_seq_null, percentile_scores, index_value,
                           discover, load_digits, SEED)

NATIVE = ["PCG64", "MT19937", "Philox", "SFC64", "python_random", "python_secrets"]


def native_digits(name, n, seed):
    if name in ("PCG64", "MT19937", "Philox", "SFC64"):
        rng = np.random.Generator(getattr(np.random, name)(seed))
        return rng.integers(0, 10, n, dtype=np.int8)
    if name == "python_random":
        import random
        r = random.Random(seed)
        return np.fromiter((r.randrange(10) for _ in range(n)), dtype=np.int8, count=n)
    if name == "python_secrets":
        import secrets
        return np.frombuffer(np.frombuffer(secrets.token_bytes(4 * n), dtype=np.uint32)
                             .astype(np.int64) % 10, dtype=np.int64).astype(np.int8)
    raise ValueError(name)


def score_stream(x, null, seq, n):
    c = np.bincount(x.astype(np.int64), minlength=10)
    o = {k: float(v[0]) for k, v in stats_from_counts(c).items()}
    o.update(r1=lag1(x), chiser=chi2_serial(x))
    _, t = percentile_scores(null, seq, o, n)
    return index_value(t)


def anova(groups):
    g = [np.asarray(v, float) for v in groups.values() if len(v) > 1]
    k = len(g)
    n_i = np.array([x.size for x in g])
    grand = np.concatenate(g).mean()
    msb = float(np.sum(n_i * (np.array([x.mean() for x in g]) - grand) ** 2)) / (k - 1)
    msw = float(np.sum([np.sum((x - x.mean()) ** 2) for x in g])) / (n_i.sum() - k)
    n0 = (n_i.sum() - np.sum(n_i ** 2) / n_i.sum()) / (k - 1)
    sig_b = max(0.0, (msb - msw) / n0)
    return {"F": msb / msw, "p_value": float(fdist.sf(msb / msw, k - 1, n_i.sum() - k)),
            "sigma_within": float(np.sqrt(msw)), "sigma_between": float(np.sqrt(sig_b)),
            "ICC": sig_b / (sig_b + msw) if sig_b + msw > 0 else 0.0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=200)
    ap.add_argument("--n", type=int, default=1_000_000, help="digits per seeded stream")
    ap.add_argument("--null", type=int, default=20_000)
    ap.add_argument("--seq-null", type=int, default=100)
    ap.add_argument("--external", default=None, help="folder with externally generated seeds")
    ap.add_argument("--out", default="results_multiseed")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    print(f"calibrating the null once, at n = {args.n:,} ...")
    null = build_null(args.n, np.random.default_rng(SEED + 20), args.null)
    seq = build_seq_null(args.n, np.random.default_rng(SEED + 21), args.seq_null,
                         n_cal=min(args.n, 1_000_000))

    results = {}
    for name in NATIVE:
        vals = []
        for s in range(args.seeds):
            x = native_digits(name, args.n, SEED + 1000 * s)
            vals.append(score_stream(x, null, seq, args.n))
        results[name] = vals
        v = np.asarray(vals)
        print(f"  {name:16s} mean {v.mean():.3f}  sd {v.std(ddof=1):.3f}  "
              f"min {v.min():.3f}  max {v.max():.3f}")

    if args.external:
        for src, files in discover(args.external).items():
            vals = []
            for p in files:
                x = load_digits(p)
                vals.append(score_stream(x[:args.n], null, seq, min(args.n, x.size)))
            results[src] = vals
            v = np.asarray(vals)
            print(f"  {src:16s} mean {v.mean():.3f}  sd {v.std(ddof=1):.3f}  ({len(files)} seeds)")

    a = anova(results)
    print(f"\ngenerator effect: F = {a['F']:.2f}, p = {a['p_value']:.4g}")
    print(f"sigma_within = {a['sigma_within']:.4f}   sigma_between = {a['sigma_between']:.4f}"
          f"   ICC = {a['ICC']:.3f}")
    if a["p_value"] > 0.05:
        print("=> at this sample size the generators are statistically indistinguishable;\n"
              "   report that as the result rather than a ranking.")

    with open(os.path.join(args.out, "multiseed.json"), "w") as fh:
        json.dump({"n": args.n, "seeds": args.seeds, "index_by_seed": results,
                   "anova": a}, fh, indent=1)
    print(f"written to {args.out}")


if __name__ == "__main__":
    main()
