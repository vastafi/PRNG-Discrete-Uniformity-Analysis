#!/usr/bin/env python3
"""
reproduce_all.py  (version 2)

Recomputes every statistic of

  "A Multi-Criterion Statistical Assessment of Discrete Uniformity in Random
   and Pseudorandom Number Sources"
  V. Astafi, A. Leahu, D. Ciorba

from the raw digit files.  Drop-in replacement for version 1: same data
discovery, same CLI, same results.json keys plus new ones.

What changed with respect to version 1
--------------------------------------
(a) Two-sided calibration.  v1 scored every criterion with
        u = P(T_null > T_obs),
    which returns u -> 1 when the observed statistic is unusually SMALL.  For
    chi-square, KS or the entropy deficit an unusually small value means
    under-dispersion, not quality, so v1 rewards agreement that is too good to
    be true.  v2 adds
        t = 2 * min( P(T_null > T_obs), P(T_null < T_obs) ),
    which is ~U(0,1) under H0 and flags both tails.  Both scores are reported.

(b) The serial chi-square and the lag-1 autocorrelation are calibrated by
    Monte Carlo on simulated digit streams instead of by their asymptotic laws.
    The overlapping-pairs statistic is slightly under-dispersed relative to
    chi2(81) (simulated variance ~152 against 162), so the asymptotic p-value
    used in v1 is mildly conservative.

(c) v1 built the null distribution of the index with
        U["r1"] = rng.random(B);  U["chiser"] = rng.random(B)
    i.e. it ASSUMED those two scores are exactly uniform under H0.  For r1 that
    is nearly true, for the serial statistic it is not (see (b)).  v2 uses the
    simulated null of both statistics.

(d) Tiers are no longer defined by a hard-coded constant (v1 used TOP_CUT =
    0.60 inside make_figures.py, a number that appears nowhere in the paper).
    v2 groups sources by overlap of the 95% confidence interval of their
    per-run index, so the tier structure is a statement about resolution.

(e) Diagnostics the referees will ask for, computed and written out:
    - the identity  H - H_hat = chi2 / (2 N ln 2), which makes the entropy
      criterion a deterministic function of the goodness-of-fit criterion;
    - the Spearman correlation matrix of the nine criteria under H0;
    - the within-source / between-source variance decomposition of the index;
    - the joint null probability of the ThreadLocalRandom configuration.

(f) Every table is written to CSV, so no rank is transcribed by hand.

Usage
-----
    python3 reproduce_all.py --data data --out results
    python3 reproduce_all.py --data data --out results --seq-null 200
"""

import argparse, csv, json, os, re, sys
import numpy as np
from scipy.stats import chi2, norm, spearmanr, t as tdist, f as fdist

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

N_NULL = 120_000          # null replicates, count-based criteria
N_SEQ_NULL = 200          # null replicates, sequence-based criteria (expensive)
SEED = 20260908
DIGITS = np.arange(10)

EX, VARX, ASX = 4.5, 8.25, 0.0
EXCX = -6 * (10 ** 2 + 1) / (5 * (10 ** 2 - 1))     # = -1.2242424...
HX = np.log2(10)                                    # = 3.3219280...

WEIGHTS = {"dmean": 0.15, "dvar": 0.15, "chi": 0.15, "chiser": 0.15,
           "dks": 0.15, "r1": 0.15,
           "dskew": 1 / 30, "dexc": 1 / 30, "dH": 1 / 30}

# the entropy criterion is a monotone transform of "chi" (see report_identity);
# this alternative drops it and redistributes its weight
WEIGHTS_NO_ENTROPY = {"dmean": 0.15, "dvar": 0.15, "chi": 0.18, "chiser": 0.18,
                      "dks": 0.15, "r1": 0.15,
                      "dskew": 0.02, "dexc": 0.02, "dH": 0.0}

COUNT_CRITERIA = ["dmean", "dvar", "dskew", "dexc", "dH", "chi", "dks"]
SEQ_CRITERIA = ["r1", "chiser"]
ALL_CRITERIA = COUNT_CRITERIA + SEQ_CRITERIA

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
# Discovery and loading  (unchanged from version 1)
# ----------------------------------------------------------------------------

def run_sort_key(filename):
    stem = os.path.splitext(filename)[0]
    m = re.match(r"^(.*?)(\d*)$", stem)
    return (m.group(1), int(m.group(2)) if m.group(2) else 0)


def discover(data_dir):
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
    with open(path, "r") as fh:
        head = fh.read(4096)
    sep = "," if head.count(",") > head.count("\n") else None
    a = np.loadtxt(path, dtype=np.int64, delimiter=sep).ravel()
    bad = a[(a < 0) | (a > 9)]
    if bad.size:
        sys.exit(f"error: {path} contains values outside 0-9 (e.g. {bad[0]})")
    return a.astype(np.int8)


# ----------------------------------------------------------------------------
# Statistics  (unchanged from version 1)
# ----------------------------------------------------------------------------

def stats_from_counts(C):
    C = np.atleast_2d(np.asarray(C, dtype=np.float64))
    n = C.sum(axis=1)
    f = C / n[:, None]
    m = (f * DIGITS).sum(axis=1)
    m2 = ((f * (DIGITS - m[:, None]) ** 2).sum(axis=1)) * n / (n - 1)
    sd = np.sqrt(m2)
    S3 = (C * (DIGITS - m[:, None]) ** 3).sum(axis=1)
    S4 = (C * (DIGITS - m[:, None]) ** 4).sum(axis=1)
    g1 = (n / ((n - 1) * (n - 2))) * S3 / sd ** 3
    g2 = (n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3)) * S4 / sd ** 4 \
         - 3 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    L = np.where(f > 0, f * np.log2(np.where(f > 0, f, 1.0)), 0.0)
    H = -L.sum(axis=1)
    E = n / 10.0
    chisq = ((C - E[:, None]) ** 2 / E[:, None]).sum(axis=1)
    D = np.max(np.abs(np.cumsum(C, axis=1) / n[:, None] - (DIGITS + 1) / 10.0), axis=1)
    return {"mean": m, "var": m2, "skew": g1, "exc": g2, "H": H,
            "dmean": np.abs(m - EX), "dvar": np.abs(m2 - VARX),
            "dskew": np.abs(g1 - ASX), "dexc": np.abs(g2 - EXCX),
            "dH": np.abs(H - HX), "chi": chisq, "dks": D}


def lag1(x):
    x = x.astype(np.float64)
    d = x - x.mean()
    return float((d[:-1] * d[1:]).sum() / (d ** 2).sum())


def chi2_serial(x):
    a = x[:-1].astype(np.int64)
    b = x[1:].astype(np.int64)
    T = np.bincount(a * 10 + b, minlength=100).reshape(10, 10).astype(np.float64)
    r = T.sum(axis=1, keepdims=True)
    c = T.sum(axis=0, keepdims=True)
    E = r @ c / T.sum()
    return float(((T - E) ** 2 / E).sum())


def chi2_between(a, b):
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
    """Multinomial null for the seven count-based criteria, drawn jointly."""
    C = np.empty((B, 10), dtype=np.int64)
    for i in range(0, B, 10_000):
        C[i:i + 10_000] = rng.multinomial(N, np.full(10, 0.1), size=min(10_000, B - i))
    return stats_from_counts(C)


def build_seq_null(N, rng, B=N_SEQ_NULL, n_cal=None):
    """Monte Carlo null for the serial chi-square and |r1|.

    The serial statistic is asymptotically free of N and |r1| is standardised
    by sqrt(N), so the null may be calibrated at n_cal <= N to save time; pass
    two values of n_cal to check that the result is stable.
    """
    n_cal = int(min(N, n_cal or N))
    s = np.empty(B)
    r = np.empty(B)
    for i in range(B):
        x = rng.integers(0, 10, n_cal, dtype=np.int8)
        s[i] = chi2_serial(x)
        r[i] = abs(lag1(x) + 1.0 / n_cal) * np.sqrt(n_cal)
        if (i + 1) % 25 == 0:
            print(f"    sequence null {i + 1}/{B}", end="\r", flush=True)
    print(" " * 40, end="\r")
    return {"chiser": s, "r1": r, "n_cal": n_cal}


def two_sided(obs, null):
    """t = 2 min(P(null > obs), P(null < obs)), the typicality score."""
    hi = float((null > obs).mean())
    lo = float((null < obs).mean())
    return float(min(1.0, 2.0 * min(hi, lo))), hi


def percentile_scores(null, seq_null, observed, N):
    """Returns (u_oneside, t_twoside) for the nine criteria.

    u_oneside reproduces version 1 exactly, so the published numbers can still
    be regenerated; t_twoside is the corrected score.
    """
    u, t = {}, {}
    for k in COUNT_CRITERIA:
        t[k], u[k] = two_sided(observed[k], null[k])
    # version-1 asymptotic scores, kept for comparison
    u["r1"] = float(2 * norm.sf(abs(observed["r1"] + 1.0 / N) * np.sqrt(N)))
    u["chiser"] = float(chi2.sf(observed["chiser"], 81))
    # version-2 Monte Carlo scores
    r1_std = abs(observed["r1"] + 1.0 / N) * np.sqrt(N)
    t["r1"] = float((seq_null["r1"] > r1_std).mean())          # already two-sided in |r1|
    t["chiser"], _ = two_sided(observed["chiser"], seq_null["chiser"])
    return u, t


def index_value(scores, weights=WEIGHTS):
    return float(sum(weights[k] * scores[k] for k in weights))


def index_geom(scores, weights=WEIGHTS, floor=1.0 / N_NULL):
    return float(np.exp(sum(weights[k] * np.log(max(scores[k], floor)) for k in weights)))


def null_index(null, seq_null, rng, weights=WEIGHTS, B=None, form="additive"):
    """Null distribution of the two-sided index.

    Count-based scores come from the ranks of the joint multinomial draws, so
    their mutual dependence is preserved; the two sequence-based scores are
    resampled from their own Monte Carlo null, which is legitimate because they
    are asymptotically independent of the digit counts.
    """
    B = B or len(null["chi"])
    T = {}
    for k in COUNT_CRITERIA:
        r = null[k].argsort().argsort()
        hi = 1.0 - (r + 0.5) / len(r)
        T[k] = np.minimum(1.0, 2.0 * np.minimum(hi, 1.0 - hi))
    for k in SEQ_CRITERIA:
        v = seq_null[k]
        r = v.argsort().argsort()
        hi = 1.0 - (r + 0.5) / len(r)
        tt = hi if k == "r1" else np.minimum(1.0, 2.0 * np.minimum(hi, 1.0 - hi))
        T[k] = rng.choice(tt, size=B, replace=True)
    if form == "geometric":
        floor = 1.0 / max(B, 1)
        return np.exp(sum(weights[k] * np.log(np.maximum(T[k], floor)) for k in weights))
    return sum(weights[k] * T[k] for k in weights)


# ----------------------------------------------------------------------------
# Diagnostics
# ----------------------------------------------------------------------------

def report_identity(null, N):
    """H - H_hat equals chi2 / (2 N ln 2) to second order; verify it."""
    deficit = HX - null["H"]
    predicted = null["chi"] / (2 * N * np.log(2))
    rho = float(spearmanr(null["dH"], null["chi"]).statistic)
    rel = float(np.max(np.abs(deficit / predicted - 1.0)))
    return {"spearman_dH_chi2": rho,
            "max_relative_error_of_identity": rel,
            "theoretical_entropy_bias": float(-(10 - 1) / (2 * N * np.log(2))),
            "null_mean_abs_entropy_deviation": float(null["dH"].mean()),
            "null_5_95_abs_entropy_deviation": [float(x) for x in
                                                np.percentile(null["dH"], [5, 95])]}


def criterion_correlations(null, seq_null, rng):
    B = min(len(null["chi"]), 20_000)
    cols = []
    for k in ALL_CRITERIA:
        v = null[k] if k in null else seq_null[k]
        cols.append(v[:B] if len(v) >= B else rng.choice(v, size=B, replace=True))
    M = spearmanr(np.column_stack(cols)).statistic
    return {"criteria": ALL_CRITERIA, "spearman": np.asarray(M).tolist()}


def variance_components(per_run):
    groups = [np.asarray(v, float) for v in per_run.values() if len(v) > 1]
    if len(groups) < 2:
        return {}
    k = len(groups)
    n_i = np.array([g.size for g in groups])
    grand = np.concatenate(groups).mean()
    msb = float(np.sum(n_i * (np.array([g.mean() for g in groups]) - grand) ** 2)) / (k - 1)
    msw = float(np.sum([np.sum((g - g.mean()) ** 2) for g in groups])) / (n_i.sum() - k)
    n0 = (n_i.sum() - np.sum(n_i ** 2) / n_i.sum()) / (k - 1)
    sig_b = max(0.0, (msb - msw) / n0)
    return {"F": msb / msw, "df_between": int(k - 1), "df_within": int(n_i.sum() - k),
            "p_value": float(fdist.sf(msb / msw, k - 1, n_i.sum() - k)),
            "sigma2_within": msw, "sigma2_between": sig_b,
            "ICC": sig_b / (sig_b + msw) if (sig_b + msw) > 0 else 0.0}


def joint_anomaly_probability(null, observed, rng):
    """Probability, under H0, of the ThreadLocalRandom configuration:
    low-order moments unusually good AND goodness-of-fit unusually bad."""
    def upper(k):
        return (null[k] > observed[k]).mean()
    um, uv, uc = upper("dmean"), upper("dvar"), upper("chi")
    sel = ((null["dmean"] < observed["dmean"]) &
           (null["dvar"] < observed["dvar"]) &
           (null["chi"] > observed["chi"]))
    p = float(sel.mean())
    return {"u_mean": float(um), "u_var": float(uv), "u_chi2": float(uc),
            "p_joint": p, "p_joint_over_10_sources": float(1 - (1 - p) ** 10)}


def tiers_from_cis(per_run, index_by_source):
    """Group sources whose 95% CI of the mean per-run index overlaps.

    Replaces the hard-coded TOP_CUT = 0.60 of version 1. The tiers are then a
    statement about what the data can resolve, not about a chosen constant.
    """
    ci = {}
    for s, v in per_run.items():
        v = np.asarray(v, float)
        if v.size > 1:
            h = tdist.ppf(0.975, v.size - 1) * v.std(ddof=1) / np.sqrt(v.size)
        else:
            h = 0.0
        ci[s] = (float(v.mean() - h), float(v.mean() + h))
    # order by the mean of the per-run indices, the quantity the CIs describe;
    # ordering by the pooled index while grouping by per-run CIs mixes two
    # different sample sizes and can put a source first and last at once
    mean_by_source = {s: float(np.mean(v)) for s, v in per_run.items()}
    order = sorted(mean_by_source, key=lambda s: -mean_by_source[s])
    tiers, current, lo_cur = {}, 1, ci[order[0]][0]
    for s in order:
        if ci[s][1] < lo_cur:          # no overlap with the current group
            current += 1
            lo_cur = ci[s][0]
        else:
            lo_cur = min(lo_cur, ci[s][0])
        tiers[s] = current
    return tiers, ci, mean_by_source


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--out", default="./results")
    ap.add_argument("--null", type=int, default=N_NULL)
    ap.add_argument("--seq-null", type=int, default=N_SEQ_NULL)
    ap.add_argument("--n-cal", type=int, default=None,
                    help="stream length for the sequence-based null (default: N)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    sources = discover(args.data)
    print(f"found {len(sources)} source(s) under {args.data}")
    for name, files in sources.items():
        print(f"  {name:24s} {len(files)} run(s)")

    runs, obs = {}, {}
    for name, files in sources.items():
        runs[name] = [load_digits(p) for p in files]
        pooled = np.concatenate(runs[name])
        counts = np.bincount(pooled.astype(np.int64), minlength=10)
        s = {k: float(v[0]) for k, v in stats_from_counts(counts).items()}
        s.update(N=int(pooled.size), n_runs=len(runs[name]), counts=counts.tolist(),
                 r1=lag1(pooled), chiser=chi2_serial(pooled))
        if len(runs[name]) >= 2:
            s["chi_between_run1_run2"] = chi2_between(runs[name][0], runs[name][1])
            s["p_between"] = float(chi2.sf(s["chi_between_run1_run2"], 81))
        s["p_gof"] = float(chi2.sf(s["chi"], 9))
        s["p_serial_asymptotic"] = float(chi2.sf(s["chiser"], 81))
        s["per_run"] = []
        for a in runs[name]:
            c1 = np.bincount(a.astype(np.int64), minlength=10)
            tt = {k: float(v[0]) for k, v in stats_from_counts(c1).items()}
            tt.update(n=int(a.size), r1=lag1(a), chiser=chi2_serial(a))
            s["per_run"].append(tt)
        obs[name] = s
        print(f"  computed {name}", flush=True)

    Ns = sorted({s["N"] for s in obs.values()})
    run_sizes = sorted({t["n"] for s in obs.values() for t in s["per_run"]})

    print("\ncalibrating count-based nulls ...", flush=True)
    nulls = {N: build_null(N, np.random.default_rng(SEED), args.null) for N in Ns}
    run_nulls = {n: build_null(n, np.random.default_rng(SEED + 1), args.null) for n in run_sizes}
    print("calibrating sequence-based nulls (this is the slow part) ...", flush=True)
    seq_nulls = {N: build_seq_null(N, np.random.default_rng(SEED + 3), args.seq_null, args.n_cal)
                 for N in Ns}
    seq_run_nulls = {n: build_seq_null(n, np.random.default_rng(SEED + 4), args.seq_null, args.n_cal)
                     for n in run_sizes}

    scores = {}
    for name, s in obs.items():
        u, t = percentile_scores(nulls[s["N"]], seq_nulls[s["N"]], s, s["N"])
        scores[name] = {
            "u_oneside": u, "t_twoside": t,
            "crqi_v1_oneside": index_value(u),
            "index_twoside": index_value(t),
            "index_twoside_geom": index_geom(t),
            "index_twoside_no_entropy": index_value(t, WEIGHTS_NO_ENTROPY),
        }

    per_run = {}
    for name, s in obs.items():
        vals = []
        for tt in s["per_run"]:
            _, t2 = percentile_scores(run_nulls[tt["n"]], seq_run_nulls[tt["n"]], tt, tt["n"])
            vals.append(index_value(t2))
        per_run[name] = vals

    rng = np.random.default_rng(SEED + 2)
    cn = null_index(nulls[Ns[0]], seq_nulls[Ns[0]], rng)
    m = len(sources)
    thr5 = float(np.quantile(cn, 0.05))
    sidak_alpha = 1 - (1 - 0.05) ** (1 / m)
    thr_sidak = float(np.quantile(cn, sidak_alpha))
    band = (float(np.quantile(cn, 0.025)), float(np.quantile(cn, 0.975)))

    idx = {k: v["index_twoside"] for k, v in scores.items()}
    tiers, cis, mean6 = tiers_from_cis(per_run, idx)

    SCHEMES = {
        "original": WEIGHTS,
        "equal": {k: 1 / 9 for k in WEIGHTS},
        "gof_heavy": {"chi": .20, "dks": .20, "chiser": .20, "dmean": .13,
                      "dvar": .13, "r1": .13, "dskew": 1 / 300, "dexc": 1 / 300,
                      "dH": 1 / 300},
        "six_criteria": {"dmean": 1 / 6, "dvar": 1 / 6, "chi": 1 / 6, "chiser": 1 / 6,
                         "dks": 1 / 6, "r1": 1 / 6, "dskew": 0.0, "dexc": 0.0, "dH": 0.0},
        "no_entropy": WEIGHTS_NO_ENTROPY,
    }
    rank_range = {}
    for name in sources:
        rr = []
        for w in SCHEMES.values():
            vals = {n2: index_value(scores[n2]["t_twoside"], w) for n2 in sources}
            rr.append(1 + sorted(vals, key=lambda x: -vals[x]).index(name))
        rank_range[name] = (min(rr), max(rr))

    key_of = {"r1": lambda s: abs(s["r1"]), "chiser": lambda s: s["chiser"]}
    ranksum = {n2: 0 for n2 in sources}
    for crit in ALL_CRITERIA:
        f = key_of.get(crit, lambda s, c=crit: s[c])
        for i, n2 in enumerate(sorted(sources, key=lambda x: f(obs[x]))):
            ranksum[n2] += i + 1

    diag = {
        "entropy_chi2_identity": report_identity(nulls[Ns[0]], Ns[0]),
        "criterion_correlations": criterion_correlations(nulls[Ns[0]], seq_nulls[Ns[0]], rng),
        "variance_components": variance_components(per_run),
    }
    tl = next((n for n in sources if "ThreadLocal" in n), None)
    if tl:
        diag["threadlocal_joint_anomaly"] = joint_anomaly_probability(nulls[obs[tl]["N"]], obs[tl], rng)

    spreads = [max(v) - min(v) for v in per_run.values() if len(v) > 1]
    within = float(np.mean(spreads)) if spreads else float("nan")
    means = [float(np.mean(v)) for v in per_run.values()]
    between = float(max(means) - min(means)) if len(means) > 1 else float("nan")

    out = {"sources": {k: [os.path.basename(p) for p in v] for k, v in sources.items()},
           "n_null": args.null, "n_seq_null": args.seq_null, "seed": SEED, "m_sources": m,
           "observed": obs, "scores": scores, "rank_range": rank_range, "rank_sum": ranksum,
           "per_run_crqi": per_run, "per_run_ci": cis, "per_run_mean": mean6, "tiers": tiers,
           "threshold_5pct": thr5, "threshold_sidak": thr_sidak, "null_band_95": band,
           "within_source_range_mean": within, "between_source_range": between,
           "diagnostics": diag}
    with open(os.path.join(args.out, "results.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    # ---- CSV tables, so that nothing is transcribed by hand ----
    order = sorted(sources, key=lambda s: -idx[s])
    with open(os.path.join(args.out, "table_criteria.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["source"] + [f"obs_{c}" for c in ALL_CRITERIA]
                   + [f"rank_{c}" for c in ALL_CRITERIA])
        ranks = {}
        for c in ALL_CRITERIA:
            f = key_of.get(c, lambda s, cc=c: s[cc])
            srt = sorted(sources, key=lambda x: f(obs[x]))
            ranks[c] = {s: i + 1 for i, s in enumerate(srt)}
        for s in order:
            f = lambda c: key_of.get(c, lambda ss, cc=c: ss[cc])(obs[s])
            w.writerow([s] + [f(c) for c in ALL_CRITERIA] + [ranks[c][s] for c in ALL_CRITERIA])

    with open(os.path.join(args.out, "table_index.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["source", "index_twoside", "per_run_mean", "index_twoside_geom",
                    "crqi_v1_oneside", "index_no_entropy", "tier", "ci_lo", "ci_hi",
                    "rank_range", "rank_sum"])
        for s in order:
            w.writerow([s, scores[s]["index_twoside"], mean6[s], scores[s]["index_twoside_geom"],
                        scores[s]["crqi_v1_oneside"], scores[s]["index_twoside_no_entropy"],
                        tiers[s], cis[s][0], cis[s][1],
                        f"{rank_range[s][0]}-{rank_range[s][1]}", ranksum[s]])

    with open(os.path.join(args.out, "table_scores.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["source"] + [f"t_{c}" for c in ALL_CRITERIA] + [f"u1_{c}" for c in ALL_CRITERIA])
        for s in order:
            w.writerow([s] + [scores[s]["t_twoside"][c] for c in ALL_CRITERIA]
                       + [scores[s]["u_oneside"][c] for c in ALL_CRITERIA])

    # ---- report ----
    d = diag["entropy_chi2_identity"]
    print(f"\nentropy/chi2 identity: Spearman = {d['spearman_dH_chi2']:.4f}, "
          f"max relative error {d['max_relative_error_of_identity']:.2e}")
    print(f"  => the entropy criterion carries no information beyond the goodness-of-fit chi2")
    print(f"entropy null centre |H_hat-H| = {d['null_mean_abs_entropy_deviation']:.3e} "
          f"(NOT zero); 5-95% = {d['null_5_95_abs_entropy_deviation']}")
    vc = diag["variance_components"]
    if vc:
        print(f"variance components: within = {vc['sigma2_within']:.5f}, "
              f"between = {vc['sigma2_between']:.5f}, ICC = {vc['ICC']:.3f}, "
              f"F = {vc['F']:.2f}, p = {vc['p_value']:.4f}")
    if tl:
        j = diag["threadlocal_joint_anomaly"]
        print(f"{tl}: joint null probability = {j['p_joint']:.2e} "
              f"({j['p_joint_over_10_sources']:.3f} over ten sources)")
    print(f"\nalarm threshold 5% = {thr5:.3f}   Sidak = {thr_sidak:.3f}   "
          f"95% null band = [{band[0]:.3f}, {band[1]:.3f}]")
    print(f"mean within-source range = {within:.3f}   between-source range = {between:.3f}")
    print(f"\n{'source':24s}{'index2':>9s}{'geom':>8s}{'CRQIv1':>9s}{'tier':>6s}{'S':>5s}")
    for s in order:
        print(f"{s:24s}{idx[s]:9.3f}{scores[s]['index_twoside_geom']:8.3f}"
              f"{scores[s]['crqi_v1_oneside']:9.3f}{tiers[s]:6d}{ranksum[s]:5d}")
    print(f"\nwritten to {args.out}")


if __name__ == "__main__":
    main()
