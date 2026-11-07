"""Self-contained implementation of the scoring used in the paper:
nine criteria, multinomial null for the count block, sequence null for the pair
block, two-sided typicality scores, additive and geometric index."""
import numpy as np

SEED = 20260908
DIGITS = np.arange(10)
EX, VARX, ASX = 4.5, 8.25, 0.0
EXCX = -6 * (10 ** 2 + 1) / (5 * (10 ** 2 - 1))
HX = np.log2(10)
WEIGHTS = {"dmean": 0.15, "dvar": 0.15, "chi": 0.15, "chiser": 0.15, "dks": 0.15, "r1": 0.15,
           "dskew": 1 / 30, "dexc": 1 / 30, "dH": 1 / 30}
COUNT = ["dmean", "dvar", "dskew", "dexc", "dH", "chi", "dks"]
# The four moment criteria are scored two-sided on the SIGNED deviation
SIGNED = {"dmean": "smean", "dvar": "svar", "dskew": "sskew", "dexc": "sexc"}
SEQ = ["r1", "chiser"]


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
    g2 = (n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3)) * S4 / sd ** 4 - 3 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    L = np.where(f > 0, f * np.log2(np.where(f > 0, f, 1.0)), 0.0)
    H = -L.sum(axis=1)
    E = n / 10.0
    chisq = ((C - E[:, None]) ** 2 / E[:, None]).sum(axis=1)
    G = 2 * np.where(C > 0, C * np.log(np.where(C > 0, C, 1.0) / E[:, None]), 0.0).sum(axis=1)
    D = np.max(np.abs(np.cumsum(C, axis=1) / n[:, None] - (DIGITS + 1) / 10.0), axis=1)
    return {"dmean": np.abs(m - EX), "dvar": np.abs(m2 - VARX), "dskew": np.abs(g1 - ASX),
            "dexc": np.abs(g2 - EXCX),
            "smean": m - EX, "svar": m2 - VARX, "sskew": g1 - ASX, "sexc": g2 - EXCX,
            "dH": np.abs(H - HX), "chi": chisq, "G": G, "dks": D}


def lag1(x):
    d = x.astype(np.float64)
    d -= d.mean()
    return float((d[:-1] * d[1:]).sum() / (d * d).sum())


def chi2_serial(x):
    a = x[:-1].astype(np.int64)
    T = np.bincount(a * 10 + x[1:], minlength=100).reshape(10, 10).astype(np.float64)
    E = T.sum(axis=1, keepdims=True) @ T.sum(axis=0, keepdims=True) / T.sum()
    return float(((T - E) ** 2 / E).sum())


def criteria_of(x):
    c = np.bincount(x, minlength=10)
    o = {k: float(v[0]) for k, v in stats_from_counts(c).items()}
    o.update(r1=lag1(x), chiser=chi2_serial(x))
    return o


def build_null(N, rng, B=120_000):
    C = np.empty((B, 10), dtype=np.int64)
    for i in range(0, B, 10_000):
        C[i:i + 10_000] = rng.multinomial(N, np.full(10, 0.1), size=min(10_000, B - i))
    return stats_from_counts(C)


def build_seq_null(rng, B=4000, n_cal=1_000_000):
    s, r = np.empty(B), np.empty(B)
    for i in range(B):
        x = rng.integers(0, 10, n_cal, dtype=np.int8)
        s[i] = chi2_serial(x)
        r[i] = abs(lag1(x) + 1.0 / n_cal) * np.sqrt(n_cal)
    return {"chiser": s, "r1": r}


def two_sided(obs, null):
    return float(min(1.0, 2.0 * min((null > obs).mean(), (null < obs).mean())))


def scores(null, seq, o, N):
    t = {k: two_sided(o[SIGNED.get(k, k)], null[SIGNED.get(k, k)]) for k in COUNT}
    t["r1"] = float((seq["r1"] > abs(o["r1"] + 1.0 / N) * np.sqrt(N)).mean())
    t["chiser"] = two_sided(o["chiser"], seq["chiser"])
    return t


def index_value(t, w=WEIGHTS):
    return float(sum(w[k] * t[k] for k in w))


def index_geom(t, w=WEIGHTS, floor=1.0 / 120_000):
    return float(np.exp(sum(w[k] * np.log(max(t[k], floor)) for k in w)))


def _mid(v):
    r = v.argsort().argsort()
    return 1.0 - (r + 0.5) / len(r)


def null_index(null, seq, rng, w=WEIGHTS, form="additive"):
    """Block-product null: count scores from the joint multinomial draws, the
    two pair scores resampled jointly from the sequence null."""
    B = len(null["chi"])
    T = {}
    for k in COUNT:
        h = _mid(null[SIGNED.get(k, k)]); T[k] = np.minimum(1.0, 2.0 * np.minimum(h, 1.0 - h))
    pick = rng.integers(0, len(seq["r1"]), size=B)
    for k in SEQ:
        h = _mid(seq[k]); tt = h if k == "r1" else np.minimum(1.0, 2.0 * np.minimum(h, 1.0 - h))
        T[k] = tt[pick]
    if form == "geometric":
        return np.exp(sum(w[k] * np.log(np.maximum(T[k], 1.0 / B)) for k in w))
    return sum(w[k] * T[k] for k in w)


def tilt(x, delta, rng):
    p = (1.0 + delta * (DIGITS - 4.5) / 4.5) / 10.0
    p /= p.sum()
    return rng.choice(10, size=x.size, p=p).astype(np.int8)


def block_sort(x, m, fraction, rng):
    y = x.copy()
    n = (y.size // m) * m
    blocks = y[:n].reshape(-1, m)
    if fraction >= 1.0:
        blocks.sort(axis=1)
    else:
        idx = rng.choice(blocks.shape[0], size=int(blocks.shape[0] * fraction), replace=False)
        blocks[idx] = np.sort(blocks[idx], axis=1)
    y[:n] = blocks.reshape(-1)
    return y


def load_joint(files):
    """Load one or more joint_null.py outputs and concatenate them in the given order."""
    zs = [np.load(f) for f in files]
    keys = [k for k in zs[0].files if k != "n_done"]
    return {k: np.concatenate([z[k][:int(z["n_done"])] for z in zs]) for k in keys}
