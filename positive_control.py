#!/usr/bin/env python3
"""
positive_control.py

Re-runs the positive control of Section 3.8 under the two-sided calibration.

A composite index is only useful if it reacts to real departures from
uniformity. Three streams are built from a clean reference stream and
deliberately corrupted, then scored against the same null reference used
throughout Section 3:

  A  marginal bias only      digit probabilities tilted linearly,
                             p_d proportional to 1 + delta * (d - 4.5) / 4.5,
                             serial behaviour untouched
  B  serial dependence only  an exact permutation of the clean stream, so the
                             digit counts and therefore all seven count-based
                             criteria are identical to the reference by
                             construction, while local reordering inside blocks
                             of length m induces lag-1 dependence
  C  both defects, weaker    the tilt at delta/10 and the block reordering
                             applied to a fraction of the stream

Both forms of the index are reported: the additive one used in the paper and
the weighted geometric one, which does not let seven typical criteria mask one
atypical criterion. That comparison is the point of re-running the control.

The corruption parameters are printed and written to the output file. They must
be reported in the paper: a positive control with undisclosed magnitudes says
nothing.

Usage
-----
    python3 positive_control.py --data data --source "Python secrets" \\
        --out results_control
    python3 positive_control.py --data data --delta 0.01 --block 4 --fraction 0.1
"""

import argparse, json, os
import numpy as np

from reproduce_all import (discover, load_digits, stats_from_counts, lag1,
                           chi2_serial, build_null, build_seq_null,
                           percentile_scores, index_value, index_geom,
                           null_index, ALL_CRITERIA, WEIGHTS, SEED)


def tilt(x, delta, rng):
    """Linearly tilted marginal, serial behaviour untouched.

    Each digit is independently resampled from the tilted law, so the counts
    move but nothing in the ordering does.
    """
    d = np.arange(10)
    p = (1.0 + delta * (d - 4.5) / 4.5) / 10.0
    p = np.clip(p, 1e-9, None)
    p /= p.sum()
    return rng.choice(10, size=x.size, p=p).astype(np.int8)


def block_sort(x, m, fraction=1.0, rng=None):
    """Exact permutation of x: sort inside blocks of length m.

    Counts are preserved exactly, so every count-based criterion is unchanged;
    the ordering inside blocks creates positive lag-1 dependence.
    """
    y = x.copy()
    n = (y.size // m) * m
    blocks = y[:n].reshape(-1, m)
    if fraction >= 1.0:
        blocks.sort(axis=1)
    else:
        k = int(blocks.shape[0] * fraction)
        idx = rng.choice(blocks.shape[0], size=k, replace=False)
        blocks[idx] = np.sort(blocks[idx], axis=1)
    y[:n] = blocks.reshape(-1)
    return y


def criteria_of(x):
    c = np.bincount(x.astype(np.int64), minlength=10)
    o = {k: float(v[0]) for k, v in stats_from_counts(c).items()}
    o.update(r1=lag1(x), chiser=chi2_serial(x))
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--source", default="Python secrets",
                    help="which source provides the clean reference stream; "
                         "use 'synthetic' for a freshly generated ideal stream")
    ap.add_argument("--reps", type=int, default=1,
                    help="repeat the control on independent clean streams and "
                         "report detection rates instead of single values; "
                         "forces synthetic clean streams when greater than 1")
    ap.add_argument("--delta", type=float, default=0.01,
                    help="tilt of the marginal for stream A")
    ap.add_argument("--block", type=int, default=4,
                    help="block length reordered in stream B")
    ap.add_argument("--fraction", type=float, default=0.1,
                    help="fraction of blocks reordered in stream C")
    ap.add_argument("--null", type=int, default=120_000)
    ap.add_argument("--seq-null", type=int, default=200)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--out", default="results_control")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    rng = np.random.default_rng(args.seed + 40)
    synthetic = args.source == "synthetic" or args.reps > 1
    if synthetic:
        sources = discover(args.data)
        first = next(iter(sources.values()))
        N = int(sum(load_digits(p).size for p in first))
        clean = rng.integers(0, 10, N, dtype=np.int8)
        print(f"clean reference: synthetic ideal stream, N = {N:,}")
    else:
        sources = discover(args.data)
        if args.source not in sources:
            raise SystemExit(f"source not found: {args.source}; available: {list(sources)}")
        clean = np.concatenate([load_digits(p) for p in sources[args.source]])
        N = clean.size
        print(f"clean reference: {args.source}, N = {N:,}")
    print(f"parameters: delta = {args.delta}, block = {args.block}, fraction = {args.fraction}")

    print("calibrating nulls ...")
    null = build_null(N, np.random.default_rng(args.seed), args.null)
    seq = build_seq_null(N, np.random.default_rng(args.seed + 3), args.seq_null,
                         n_cal=min(N, 3_000_000))
    cn = null_index(null, seq, np.random.default_rng(args.seed + 2))
    cg = null_index(null, seq, np.random.default_rng(args.seed + 2), form="geometric")
    thr = float(np.quantile(cn, 0.05))
    thr_geo = float(np.quantile(cg, 0.05))
    print(f"alarm threshold at 5%: additive {thr:.3f}, geometric {thr_geo:.3f}\n")

    def build(ref):
        return {"clean": ref,
                "A_marginal_bias": tilt(ref, args.delta, rng),
                "B_serial_only": block_sort(ref, args.block, 1.0, rng),
                "C_both_weaker": block_sort(tilt(ref, args.delta / 10.0, rng),
                                            args.block, args.fraction, rng)}

    names = ["clean", "A_marginal_bias", "B_serial_only", "C_both_weaker"]
    acc = {n: {"add": [], "geo": [], "alarm_add": 0, "alarm_geo": 0, "scores": None}
           for n in names}
    counts_ok = True
    for rep in range(args.reps):
        ref = clean if rep == 0 else rng.integers(0, 10, N, dtype=np.int8)
        streams = build(ref)
        counts_ok &= bool(np.array_equal(
            np.bincount(streams["B_serial_only"].astype(np.int64), minlength=10),
            np.bincount(ref.astype(np.int64), minlength=10)))
        for n in names:
            o = criteria_of(streams[n])
            _, t = percentile_scores(null, seq, o, N)
            a, g = index_value(t), index_geom(t)
            acc[n]["add"].append(a); acc[n]["geo"].append(g)
            acc[n]["alarm_add"] += int(a < thr); acc[n]["alarm_geo"] += int(g < thr_geo)
            if rep == 0:
                acc[n]["scores"] = t
        if args.reps > 1:
            print(f"  replicate {rep + 1}/{args.reps}", end="\r", flush=True)
    print(" " * 30, end="\r")
    print(f"stream B preserves the digit counts exactly: {counts_ok}\n")

    R = args.reps
    if R == 1:
        for n in names:
            a, g = acc[n]["add"][0], acc[n]["geo"][0]
            print(f"  {n:18s} additive {a:.3f} {'ALARM' if a < thr else '     '}"
                  f"   geometric {g:.3f} {'ALARM' if g < thr_geo else ''}")
    else:
        print(f"  {'stream':18s}{'additive mean':>15s}{'detected':>10s}"
              f"{'geometric mean':>16s}{'detected':>10s}")
        for n in names:
            print(f"  {n:18s}{np.mean(acc[n]['add']):15.3f}"
                  f"{acc[n]['alarm_add'] / R:9.0%}"
                  f"{np.mean(acc[n]['geo']):16.3f}{acc[n]['alarm_geo'] / R:9.0%}")
        print("\n  the 'clean' row is the false alarm rate and should sit near 5%;")
        print("  the other rows are the power of the index against each defect.")

    print(f"\n{'stream':18s}" + "".join(f"{c:>9s}" for c in ALL_CRITERIA))
    for n in names:
        print(f"{n:18s}" + "".join(f"{acc[n]['scores'][c]:9.3f}" for c in ALL_CRITERIA))

    rows = {n: {"index_additive_mean": float(np.mean(acc[n]["add"])),
                "index_geometric_mean": float(np.mean(acc[n]["geo"])),
                "detection_rate_additive": acc[n]["alarm_add"] / R,
                "detection_rate_geometric": acc[n]["alarm_geo"] / R,
                "index_additive_all": acc[n]["add"],
                "index_geometric_all": acc[n]["geo"],
                "scores_first_replicate": acc[n]["scores"]} for n in names}
    same = counts_ok

    out = {"source": args.source, "N": int(N), "delta": args.delta, "block": args.block,
           "fraction": args.fraction, "threshold_5pct": thr, "threshold_5pct_geometric": thr_geo, "reps": args.reps,
           "counts_preserved_in_B": bool(same), "weights": WEIGHTS, "streams": rows}
    with open(os.path.join(args.out, "positive_control.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwritten to {args.out}")


if __name__ == "__main__":
    main()
