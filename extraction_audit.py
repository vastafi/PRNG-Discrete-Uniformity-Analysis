#!/usr/bin/env python3
"""
extraction_audit.py

Answers the question a referee will ask first: are you testing the generator,
or the generator plus the code that turns its output into decimal digits?

The script holds the ALGORITHM fixed and varies only the extraction path, then
scores every path with the same nine criteria used in the paper.  If the paths
differ by as much as the ten sources differ from each other, then the ranking
of the paper is partly a ranking of extraction code, and that has to be said in
the manuscript.

Five extraction paths, all applied to numpy's PCG64 with the same seed:

  A  integers(0, 10)            direct uniform integer draw, no bias
  B  floor(10 * random())       float in [0,1) scaled and truncated
  C  integers(0, 2**31) % 10    modulo reduction, the classic biased path
  D  digit k of random()        one decimal digit of the float's expansion
  E  digits of "%.15f"          all 15 printed decimals, as a spreadsheet export

Paths D and E are the ones that resemble how digits are obtained from
environments that only expose a floating point uniform, which is the case for
several of the ten sources in the paper.

Usage
-----
    python3 extraction_audit.py --n 2000000 --out results_extraction
    python3 extraction_audit.py --n 2000000 --generator MT19937
"""

import argparse, json, os
import numpy as np

from reproduce_all import (stats_from_counts, lag1, chi2_serial, build_null,
                           build_seq_null, percentile_scores, index_value,
                           ALL_CRITERIA, SEED)


def make_digits(path, n, rng):
    if path == "A_integers":
        return rng.integers(0, 10, n, dtype=np.int8)
    if path == "B_float_scaled":
        return np.floor(10 * rng.random(n)).astype(np.int8)
    if path == "C_modulo":
        return (rng.integers(0, 2 ** 31, n) % 10).astype(np.int8)
    if path == "D_float_digit":
        u = rng.random(n)
        return np.floor(u * 10 ** 6 % 10).astype(np.int8)     # the 6th decimal
    if path == "E_printed_digits":
        m = max(1, n // 15)
        u = rng.random(m)
        s = "".join(f"{v:.15f}"[2:] for v in u)
        return np.frombuffer(s.encode(), dtype=np.uint8)[:n].astype(np.int8) - ord("0")
    raise ValueError(path)


PATHS = ["A_integers", "B_float_scaled", "C_modulo", "D_float_digit", "E_printed_digits"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=2_000_000)
    ap.add_argument("--generator", default="PCG64",
                    choices=["PCG64", "MT19937", "Philox", "SFC64"])
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--null", type=int, default=20_000)
    ap.add_argument("--seq-null", type=int, default=100)
    ap.add_argument("--out", default="results_extraction")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    bg = getattr(np.random, args.generator)
    print(f"generator {args.generator}, n = {args.n:,} digits per extraction path")

    print("calibrating nulls ...")
    null = build_null(args.n, np.random.default_rng(args.seed + 10), args.null)
    seq = build_seq_null(args.n, np.random.default_rng(args.seed + 11), args.seq_null,
                         n_cal=min(args.n, 1_000_000))

    rows = {}
    for p in PATHS:
        rng = np.random.Generator(bg(args.seed))
        x = make_digits(p, args.n, rng)
        c = np.bincount(x.astype(np.int64), minlength=10)
        o = {k: float(v[0]) for k, v in stats_from_counts(c).items()}
        o.update(r1=lag1(x), chiser=chi2_serial(x))
        _, t = percentile_scores(null, seq, o, args.n)
        rows[p] = {"index_twoside": index_value(t),
                   "scores": t,
                   "chi2_gof": o["chi"], "ks": o["dks"], "r1": o["r1"]}
        print(f"  {p:18s} index = {rows[p]['index_twoside']:.3f}   "
              f"chi2 = {o['chi']:8.2f}   D = {o['dks']:.6f}")

    vals = [r["index_twoside"] for r in rows.values()]
    spread = max(vals) - min(vals)
    print(f"\nspread of the index across extraction paths, same algorithm: {spread:.3f}")
    print("compare this with the between-source range reported in results.json;")
    print("if it is of the same order, the paper cannot attribute differences to")
    print("the generators without documenting the extraction path of each source.")

    with open(os.path.join(args.out, "extraction_audit.json"), "w") as fh:
        json.dump({"generator": args.generator, "n": args.n, "rows": rows,
                   "index_spread_across_paths": spread}, fh, indent=1)
    print(f"written to {args.out}")


if __name__ == "__main__":
    main()
