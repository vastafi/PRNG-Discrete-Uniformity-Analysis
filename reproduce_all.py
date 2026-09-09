#!/usr/bin/env python3
"""
Reproduction script for:
  "Concordance with the Discrete Uniform Distribution of Pseudorandom Number
   Generators in the Most Commonly Used Computer Applications"
  V. Astafi, A. Leahu, D. Ciorba

Recomputes every statistic in the paper from the raw digit files.

The script DISCOVERS the data layout by itself. Point --data at any folder whose
subfolders each hold the runs of one source; file and folder names do not matter:

    data/                            data/
      java_threadlocal/                Java thread local/
        threadlocal_data.csv             threadlocal_data.csv
        threadlocal_data1.csv            threadlocal_data1.csv
        ...                              ...
      pi/                              Pi/
      python_secrets/                  Python secrets/

Within a folder, the file without a trailing number is run 1 and the numbered
files follow in order. Any number of sources and any number of runs per source
will work; the paper uses ten sources of six runs each.

Usage
-----
    python3 reproduce_all.py --data data --out results

Requirements: numpy, scipy.
Runtime: roughly 1.5 minutes per source on a current laptop.
"""

import argparse, json, os, re, sys
import numpy as np
from scipy.stats import chi2, norm

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

N_NULL = 120_000          # null replicates
SEED   = 20260908
DIGITS = np.arange(10)

# theoretical characteristics of X ~ U{0,...,9}
EX, VARX, ASX = 4.5, 8.25, 0.0
EXCX = -6 * (10**2 + 1) / (5 * (10**2 - 1))     # = -1.2242424...
HX   = np.log2(10)                              # = 3.3219280...

# nine criteria: six main (w = 0.15) and three shape (w = 1/30)
WEIGHTS = {"dmean": 0.15, "dvar": 0.15, "chi": 0.15, "chiser": 0.15,
           "dks": 0.15, "r1": 0.15,
           "dskew": 1/30, "dexc": 1/30, "dH": 1/30}
COUNT_CRITERIA = ["dmean", "dvar", "dskew", "dexc", "dH", "chi", "dks"]

# folder name (normalised) -> label used in the paper
LABELS = {
    "cpp": "C++", "c": "C++", "cplusplus": "C++", "cpprandom": "C++",
    "excel": "Excel", "microsoftexcel": "Excel", "excelnum": "Excel",
    "javasecurerandom": "Java SecureRandom", "javasecure": "Java SecureRandom",
    "securerandom": "Java SecureRandom",
    "javathreadlocal": "Java ThreadLocalRandom",
    "javathreadlocalrandom": "Java ThreadLocalRandom", "threadlocal": "Java ThreadLocalRandom",
    "mathematica": "Mathematica", "wolframmathematica": "Mathematica",
    "pi": "pi", "pidigits": "pi",
    "pythonnumpy": "Python numpy", "numpy": "Python numpy", "pythonnumpyrandom": "Python numpy",
    "pythonsecrets": "Python secrets", "secrets": "Python secrets",
    "r": "R", "rnumbers": "R",
    "randomorg": "Random.org", "random": "Random.org",
}

def label_for(folder_name):
    key = re.sub(r"[^a-z0-9]", "", folder_name.lower())
    return LABELS.get(key, folder_name)

# ----------------------------------------------------------------------------
# Discovery and loading
# ----------------------------------------------------------------------------

def run_sort_key(filename):
    """Order runs: the un-numbered file first, then 1, 2, 3, ..."""
    stem = os.path.splitext(filename)[0]
    m = re.match(r"^(.*?)(\d*)$", stem)
    return (m.group(1), int(m.group(2)) if m.group(2) else 0)

def discover(data_dir):
    """Map each subfolder of data_dir to its ordered list of run files."""
    if not os.path.isdir(data_dir):
        sys.exit(f"error: --data path is not a folder: {data_dir}")
    found = {}
    for entry in sorted(os.listdir(data_dir)):
        folder = os.path.join(data_dir, entry)
        if not os.path.isdir(folder) or entry.startswith("."):
            continue
        files = [f for f in os.listdir(folder)
                 if f.lower().endswith((".csv", ".txt", ".dat")) and not f.startswith(".")]
        if not files:
            continue
        found[label_for(entry)] = [os.path.join(folder, f)
                                   for f in sorted(files, key=run_sort_key)]
    if not found:
        sys.exit(f"error: no source subfolders with data files under {data_dir}")
    return found

def load_digits(path):
    """Read one run. Accepts one digit per line, or comma/space separated."""
    with open(path, "r") as fh:
        head = fh.read(4096)
    sep = "," if head.count(",") > head.count("\n") else None
    a = np.loadtxt(path, dtype=np.int64, delimiter=sep).ravel()
    bad = a[(a < 0) | (a > 9)]
    if bad.size:
        sys.exit(f"error: {path} contains values outside 0-9 (e.g. {bad[0]})")
    return a.astype(np.int8)

# ----------------------------------------------------------------------------
# Statistics
# ----------------------------------------------------------------------------

def stats_from_counts(C):
    """C has shape (B, 10). Returns a dict of arrays of length B."""
    C = np.atleast_2d(np.asarray(C, dtype=np.float64))
    n = C.sum(axis=1)
    f = C / n[:, None]
    m = (f * DIGITS).sum(axis=1)                                         # equation (1)
    m2 = ((f * (DIGITS - m[:, None]) ** 2).sum(axis=1)) * n / (n - 1)    # equation (2)
    sd = np.sqrt(m2)
    S3 = (C * (DIGITS - m[:, None]) ** 3).sum(axis=1)
    S4 = (C * (DIGITS - m[:, None]) ** 4).sum(axis=1)
    g1 = (n / ((n - 1) * (n - 2))) * S3 / sd ** 3                        # equation (3)
    g2 = (n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3)) * S4 / sd ** 4 \
         - 3 * (n - 1) ** 2 / ((n - 2) * (n - 3))                        # equation (4)
    L = np.where(f > 0, f * np.log2(np.where(f > 0, f, 1.0)), 0.0)
    H = -L.sum(axis=1)                                                   # equation (5)
    E = n / 10.0
    chisq = ((C - E[:, None]) ** 2 / E[:, None]).sum(axis=1)
    D = np.max(np.abs(np.cumsum(C, axis=1) / n[:, None] - (DIGITS + 1) / 10.0), axis=1)
    return {"mean": m, "var": m2, "skew": g1, "exc": g2, "H": H,
            "dmean": np.abs(m - EX), "dvar": np.abs(m2 - VARX),
            "dskew": np.abs(g1 - ASX), "dexc": np.abs(g2 - EXCX),
            "dH": np.abs(H - HX), "chi": chisq, "dks": D}

def lag1(x):                                                             # equation (6)
    x = x.astype(np.float64)
    d = x - x.mean()
    return float((d[:-1] * d[1:]).sum() / (d ** 2).sum())

def chi2_serial(x):
    """10 x 10 contingency table over consecutive pairs, df = 81."""
    a = x[:-1].astype(np.int64)
    b = x[1:].astype(np.int64)
    T = np.bincount(a * 10 + b, minlength=100).reshape(10, 10).astype(np.float64)
    r = T.sum(axis=1, keepdims=True)
    c = T.sum(axis=0, keepdims=True)
    E = r @ c / T.sum()
    return float(((T - E) ** 2 / E).sum())

def chi2_between(a, b):
    """Independence between two distinct runs of the same source, df = 81."""
    m = min(a.size, b.size)
    T = np.bincount(a[:m].astype(np.int64) * 10 + b[:m].astype(np.int64),
                    minlength=100).reshape(10, 10).astype(np.float64)
    r = T.sum(axis=1, keepdims=True)
    c = T.sum(axis=0, keepdims=True)
    E = r @ c / T.sum()
    return float(((T - E) ** 2 / E).sum())

# ----------------------------------------------------------------------------
# Null calibration
# ----------------------------------------------------------------------------

def build_null(N, rng, B=N_NULL):
    """B multinomial count vectors at sample size N, drawn jointly so that the
    dependence among the seven count-based criteria is preserved."""
    C = np.empty((B, 10), dtype=np.int64)
    for i in range(0, B, 10_000):
        C[i:i + 10_000] = rng.multinomial(N, np.full(10, 0.1), size=min(10_000, B - i))
    return stats_from_counts(C)

def percentile_scores(null, observed, N):
    """u = proportion of null draws worse than the observed value (equation 9)."""
    u = {k: float((null[k] > observed[k]).mean()) for k in COUNT_CRITERIA}
    u["r1"] = float(2 * norm.sf(abs(observed["r1"] + 1.0 / N) * np.sqrt(N)))
    u["chiser"] = float(chi2.sf(observed["chiser"], 81))
    return u

def crqi(u, weights=WEIGHTS):                                            # equation (10)
    return sum(weights[k] * u[k] for k in weights)

def null_crqi(null, weights=WEIGHTS):
    """Null distribution of the index itself, for the alarm threshold.

    Uses its own random stream (SEED + 2) rather than a shared generator, so
    that the thresholds do not depend on how much randomness earlier steps
    happened to consume; the values are then identical on every run."""
    rng = np.random.default_rng(SEED + 2)
    B = len(null["chi"])
    U = {}
    for k in COUNT_CRITERIA:
        order = null[k].argsort()
        ranks = np.empty(B)
        ranks[order] = np.arange(B)
        U[k] = 1.0 - ranks / (B - 1)          # large statistic -> small u
    U["r1"] = rng.random(B)                    # asymptotically independent
    U["chiser"] = rng.random(B)
    return sum(weights[k] * U[k] for k in weights)

# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="folder holding one subfolder per source")
    ap.add_argument("--out", default="./results")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    sources = discover(args.data)
    print(f"found {len(sources)} source(s) under {args.data}:")
    for name, files in sources.items():
        print(f"  {name:24s} {len(files)} run(s)  [{os.path.basename(files[0])} ...]")
    print()

    # ---- load, pool and describe ----
    runs, obs = {}, {}
    for name, files in sources.items():
        runs[name] = [load_digits(p) for p in files]
        sizes = {a.size for a in runs[name]}
        if len(sizes) > 1:
            print(f"  note: {name} has runs of differing length {sorted(sizes)}")
        pooled = np.concatenate(runs[name])
        counts = np.bincount(pooled.astype(np.int64), minlength=10)
        s = {k: float(v[0]) for k, v in stats_from_counts(counts).items()}
        s["N"] = int(pooled.size)
        s["n_runs"] = len(runs[name])
        s["counts"] = counts.tolist()
        s["r1"] = lag1(pooled)
        s["chiser"] = chi2_serial(pooled)
        if len(runs[name]) >= 2:
            s["chi_between_run1_run2"] = chi2_between(runs[name][0], runs[name][1])
            s["p_between"] = float(chi2.sf(s["chi_between_run1_run2"], 81))
        s["p_gof"] = float(chi2.sf(s["chi"], 9))          # df = k - 1 = 9
        s["p_serial"] = float(chi2.sf(s["chiser"], 81))
        s["per_run"] = []
        for a in runs[name]:
            c1 = np.bincount(a.astype(np.int64), minlength=10)
            t = {k: float(v[0]) for k, v in stats_from_counts(c1).items()}
            t["n"] = int(a.size)
            t["r1"] = lag1(a)
            t["chiser"] = chi2_serial(a)
            t["p_gof"] = float(chi2.sf(t["chi"], 9))
            s["per_run"].append(t)
        obs[name] = s
        print(f"  computed {name}", flush=True)

    Ns = sorted({s["N"] for s in obs.values()})
    if len(Ns) > 1:
        print(f"\nnote: sources have different pooled sizes {Ns};")
        print("      each is calibrated against a null distribution at its own size.")

    # ---- composite index ----
    print("\ncalibrating null distributions ...", flush=True)
    nulls = {N: build_null(N, np.random.default_rng(SEED)) for N in Ns}
    scores = {name: None for name in sources}
    for name, s in obs.items():
        u = percentile_scores(nulls[s["N"]], s, s["N"])
        scores[name] = {"u": u, "crqi": crqi(u)}

    cn = null_crqi(nulls[Ns[0]])
    m = len(sources)
    thr5 = float(np.quantile(cn, 0.05))
    sidak_alpha = 1 - (1 - 0.05) ** (1 / m)                              # equation (11)
    thr_sidak = float(np.quantile(cn, sidak_alpha))
    band = (float(np.quantile(cn, 0.025)), float(np.quantile(cn, 0.975)))

    # ---- sensitivity to the weights ----
    SCHEMES = {
        "original": WEIGHTS,
        "equal": {k: 1/9 for k in WEIGHTS},
        "gof_heavy": {"chi": .20, "dks": .20, "chiser": .20,
                      "dmean": .13, "dvar": .13, "r1": .13,
                      "dskew": 1/300, "dexc": 1/300, "dH": 1/300},
        "six_criteria": {"dmean": 1/6, "dvar": 1/6, "chi": 1/6, "chiser": 1/6,
                         "dks": 1/6, "r1": 1/6, "dskew": 0.0, "dexc": 0.0, "dH": 0.0},
    }
    rank_range = {}
    for name in sources:
        rr = []
        for w in SCHEMES.values():
            vals = {n2: crqi(scores[n2]["u"], w) for n2 in sources}
            rr.append(1 + sorted(vals, key=lambda x: -vals[x]).index(name))
        rank_range[name] = (min(rr), max(rr))

    # ---- rank-sum baseline ----
    key_of = {"r1": lambda s: abs(s["r1"]), "chiser": lambda s: s["chiser"]}
    ranksum = {n2: 0 for n2 in sources}
    for crit in ["dmean", "dvar", "dskew", "dexc", "dH", "chi", "dks", "chiser", "r1"]:
        f = key_of.get(crit, lambda s, c=crit: s[c])
        for i, n2 in enumerate(sorted(sources, key=lambda x: f(obs[x]))):
            ranksum[n2] += i + 1

    # ---- per-run index, recalibrated at the per-run size ----
    run_sizes = sorted({t["n"] for s in obs.values() for t in s["per_run"]})
    run_nulls = {n: build_null(n, np.random.default_rng(SEED + 1)) for n in run_sizes}
    per_run_crqi = {name: [crqi(percentile_scores(run_nulls[t["n"]], t, t["n"]))
                           for t in s["per_run"]] for name, s in obs.items()}
    spreads = [max(v) - min(v) for v in per_run_crqi.values() if len(v) > 1]
    within = float(np.mean(spreads)) if spreads else float("nan")
    means = [float(np.mean(v)) for v in per_run_crqi.values()]
    between = float(max(means) - min(means)) if len(means) > 1 else float("nan")

    # ---- output ----
    out = {"sources": {k: [os.path.basename(p) for p in v] for k, v in sources.items()},
           "n_null": N_NULL, "seed": SEED, "m_sources": m,
           "observed": obs, "scores": scores, "rank_range": rank_range,
           "rank_sum": ranksum, "per_run_crqi": per_run_crqi,
           "threshold_5pct": thr5, "threshold_sidak": thr_sidak, "null_band_95": band,
           "within_source_range_mean": within, "between_source_range": between}
    path = os.path.join(args.out, "results.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)

    # ---- report ----
    print(f"\nalarm threshold 5% = {thr5:.3f}   Sidak (m = {m}) = {thr_sidak:.3f}")
    print(f"95% null band = [{band[0]:.3f}, {band[1]:.3f}]")
    if not np.isnan(within):
        print(f"mean within-source range = {within:.3f}   "
              f"between-source range of means = {between:.3f}")

    print(f"\n{'source':24s}{'N':>12s}{'chi2 gof':>10s}{'p':>8s}{'D_KS':>10s}"
          f"{'chi2 ser':>10s}{'p':>8s}{'CRQI':>8s}{'S':>5s}")
    for name in sorted(sources, key=lambda x: -scores[x]["crqi"]):
        s = obs[name]
        print(f"{name:24s}{s['N']:12,d}{s['chi']:10.3f}{s['p_gof']:8.4f}{s['dks']:10.6f}"
              f"{s['chiser']:10.2f}{s['p_serial']:8.4f}{scores[name]['crqi']:8.3f}"
              f"{ranksum[name]:5d}")
    print(f"\nwritten to {path}")


if __name__ == "__main__":
    main()