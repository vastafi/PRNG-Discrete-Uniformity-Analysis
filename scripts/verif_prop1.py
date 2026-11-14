import itertools, numpy as np
from fractions import Fraction as F

# 1) exact enumeration: Cov(O_l, sqrt(N) W_ij) for small k, N, compared with the closed form
def exact(k, N):
    p = F(1, k); seqs = list(itertools.product(range(k), repeat=N)); M = len(seqs)
    worst = 0
    for l in range(k):
        for i in range(k):
            for j in range(k):
                EO = EW = EOW = F(0)
                for x in seqs:
                    O = [x.count(a) for a in range(k)]
                    T = sum(1 for t in range(N-1) if x[t] == i and x[t+1] == j)
                    W = T - (N-1)*p*p - p*(O[i]-N*p) - p*(O[j]-N*p)
                    EO += O[l]; EW += W; EOW += O[l]*W
                cov = EOW/M - (EO/M)*(EW/M)
                closed = 2*p**3 - ((l == i) + (l == j))*p**2
                assert cov == closed, (k, N, l, i, j, cov, closed)
    return True
for k, N in [(2, 6), (3, 5), (4, 4)]:
    print('exact enumeration k=%d N=%d: Cov(O_l, sqrt(N) W_ij) = 2p^3-(d_li+d_lj)p^2 for all l,i,j:' % (k, N), exact(k, N))

# 2) Monte Carlo at k=10: Cov(Z_l, W_ij) should be (2p^3-(d_li+d_lj)p^2)/N -> 0
rng = np.random.default_rng(1); k = 10; p = 1/k; N = 100_000; R = 4000
Z = np.empty((R, k)); Wd = np.empty((R, k*k)); r1 = np.empty(R); lin = np.empty(R)
mu, s2 = 4.5, 8.25; dev = np.arange(k) - mu
for r in range(R):
    x = rng.integers(0, k, N)
    O = np.bincount(x, minlength=k)
    T = np.bincount(x[:-1]*k + x[1:], minlength=k*k).reshape(k, k)
    W = (T - (N-1)*p*p - p*(O[:, None]-N*p) - p*(O[None, :]-N*p)) / np.sqrt(N)
    Z[r] = (O - N*p)/np.sqrt(N); Wd[r] = W.ravel()
    xc = x - x.mean(); r1[r] = (xc[:-1]*xc[1:]).sum()/(xc*xc).sum()
    lin[r] = (np.outer(dev, dev)*W).sum()/s2
C = np.corrcoef(np.hstack([Z, Wd]).T)[:k, k:]
print('MC k=10, N=1e5, R=4000: max |corr(Z_l, W_ij)| over 1000 pairs = %.3f (noise level ~ %.3f)' % (np.abs(C).max(), 3.3/np.sqrt(R)))
print('linearisation: corr(sqrt(N)(r1+1/N), sum(i-mu)(j-mu)W_ij/sigma^2) = %.5f' % np.corrcoef(np.sqrt(N)*(r1+1/N), lin)[0, 1])
G = 2*np.array([np.sum(np.where(z*np.sqrt(N)+N*p > 0, (z*np.sqrt(N)+N*p)*np.log((z*np.sqrt(N)+N*p)/(N*p)), 0)) for z in Z])
from scipy import stats
print('Spearman(G, |r1|) = %.3f, p = %.2f' % stats.spearmanr(G, np.abs(r1)))
