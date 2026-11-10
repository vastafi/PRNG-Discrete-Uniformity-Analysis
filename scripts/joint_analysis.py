"""Compare the block-product null of the index (used in the paper) with the
joint null obtained from complete sequences, and test block independence."""
import argparse, json, sys
import numpy as np
from scipy.stats import spearmanr, chi2_contingency, ks_2samp
from common import SEED, COUNT, SEQ, build_null, build_seq_null, scores, index_value, null_index, load_joint

def dcor(x, y):
    def cent(a):
        D = np.abs(a[:, None, :] - a[None, :, :]).sum(-1) if a.ndim == 2 else np.abs(a[:, None] - a[None, :])
        return D - D.mean(0) - D.mean(1)[:, None] + D.mean()
    A, B = cent(x), cent(y)
    v = (A * B).mean(); return float(np.sqrt(max(v, 0) / np.sqrt((A * A).mean() * (B * B).mean())))

def rank01(v):
    return (v.argsort().argsort() + 0.5) / len(v)

ap = argparse.ArgumentParser()
ap.add_argument("--joint-run", nargs="+", required=True)
ap.add_argument("--joint-pooled", nargs="+", required=True, help="in order")
ap.add_argument("--out", default="results/extra/joint_analysis.json")
args = ap.parse_args()
out = {}
seq = build_seq_null(np.random.default_rng(SEED + 3))
for lab, files, N in [("run", args.joint_run, 2_614_690), ("pooled", args.joint_pooled, 15_688_140)]:
    J = load_joint(files); n = len(J["chi"])
    J_abs = dict(J); J_abs["r1"] = np.abs(J["r1"])  # the autocorrelation is scored as |r1|
    null = build_null(N, np.random.default_rng(SEED))
    # index of every complete null sequence, scored exactly as an observed source
    idx = np.array([index_value(scores(null, seq, {k: J[k][b] for k in J}, N)) for b in range(n)])
    prod = null_index(null, seq, np.random.default_rng(SEED + 2))
    q = lambda v, p: float(np.quantile(v, p))
    boot = np.random.default_rng(1)
    bq = [np.quantile(boot.choice(idx, n), 0.05) for _ in range(2000)]
    # cross-block dependence
    cnt = np.column_stack([rank01(J[k]) for k in COUNT]); sq = np.column_stack([rank01(J_abs[k]) for k in SEQ])
    sp = [[float(spearmanr(J[a], J_abs[b]).statistic) for b in SEQ] for a in COUNT]
    pvals = []
    for a in COUNT:
        for b in SEQ:
            T = np.histogram2d(rank01(J[a]), rank01(J_abs[b]), bins=5, range=[[0, 1], [0, 1]])[0]
            pvals.append(float(chi2_contingency(T)[1]))
    m = min(n, 2000); dc = dcor(cnt[:m], sq[:m])
    perm = boot
    dperm = [dcor(cnt[:m], sq[:m][perm.permutation(m)]) for _ in range(99)]
    out[lab] = {"n_joint": n,
                "thr5_joint": q(idx, 0.05), "thr5_joint_boot_se": float(np.std(bq)),
                "thr5_product": q(prod, 0.05),
                "band_joint": [q(idx, 0.025), q(idx, 0.975)], "band_product": [q(prod, 0.025), q(prod, 0.975)],
                "mean_joint": float(idx.mean()), "mean_product": float(prod.mean()),
                "sd_joint": float(idx.std(ddof=1)), "sd_product": float(prod.std(ddof=1)),
                "ks_joint_vs_product_p": float(ks_2samp(idx, prod).pvalue),
                "false_alarm_joint_at_product_thr": float((idx < q(prod, 0.05)).mean()),
                "cross_spearman_max_abs": float(np.max(np.abs(sp))),
                "cross_chi2_min_p": min(pvals), "cross_chi2_bonf_min_p": min(1.0, min(pvals) * len(pvals)),
                "dcor_blocks": dc, "dcor_perm_p": float((1 + sum(x >= dc for x in dperm)) / 100)}
    if lab == "pooled" and n >= 4000:
        # second, independent subset of 2,000 sequences (paper, Section 3.1)
        sel = slice(2000, 4000)
        c2 = np.column_stack([rank01(J[k])[sel] for k in COUNT])
        s2 = np.column_stack([rank01(J_abs[k])[sel] for k in SEQ])
        d2 = dcor(c2, s2); r5 = np.random.default_rng(5)
        p2 = [dcor(c2, s2[r5.permutation(2000)]) for _ in range(199)]
        out[lab]["dcor_blocks_subset2"] = d2
        out[lab]["dcor_perm_p_subset2"] = float((1 + sum(x >= d2 for x in p2)) / 200)
    print(lab, json.dumps(out[lab], indent=1), flush=True)
json.dump(out, open(args.out, "w"), indent=1)
