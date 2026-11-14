"""exp2_attack.py -- passive key-recovery attack from flooded decoding leakage.

Toy ring-LWE PKE, n=64, ternary secret/noise/randomness.  Each leakage query
gives the adversary an honest ciphertext (with its encryption randomness) and
the flooded exposure  y = eta + f,  eta = e*r + e2 - s*e1 (ring products).
The adversary solves least squares for x=(e,s) and rounds.  We record the
first number of queries N at which s is recovered exactly."""
import csv, math, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from floodlib import *

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
n = 64
V_ETA = n * (TERN * TERN * 2) + TERN
IDX = (np.arange(n)[:, None] - np.arange(n)[None, :]) % n
SGN = np.where(np.arange(n)[:, None] >= np.arange(n)[None, :], 1, -1)


def batch_M(r, e1):
    """r,e1: (K,n) -> M: (K,n,2n) with eta = M @ [e; s] + e2."""
    Rr = r[:, IDX] * SGN
    R1 = -(e1[:, IDX] * SGN)
    return np.concatenate([Rr, R1], axis=2).astype(np.float64)


def one_trial(sigma, rng, n_cap, grid):
    x = rng.integers(-1, 2, 2 * n)            # x = (e, s)
    xf = x.astype(np.float64)
    G = np.zeros((2 * n, 2 * n)); h = np.zeros(2 * n)
    N = 0
    for target in grid:
        if target > n_cap:
            break
        while N < target:
            K = min(target - N, 256)
            r = rng.integers(-1, 2, (K, n)); e1 = rng.integers(-1, 2, (K, n)); e2 = rng.integers(-1, 2, (K, n))
            M = batch_M(r, e1)
            eta = M @ xf + e2
            y = eta + (np.rint(rng.normal(0, sigma, eta.shape)) if sigma > 0 else 0)
            Mm = M.reshape(-1, 2 * n)
            G += Mm.T @ Mm
            h += Mm.T @ y.reshape(-1)
            N += K
        try:
            xh = np.linalg.solve(G + 1e-9 * np.eye(2 * n), h)
        except np.linalg.LinAlgError:
            continue
        if np.array_equal(np.rint(xh[n:]).astype(int), x[n:]):
            return N
    return None


def main():
    rng = np.random.default_rng(777)
    grid = sorted(set(int(round(1.2 ** i)) for i in range(2, 80)))
    rows = []
    cap = 150000
    for sigma, trials in [(0, 40), (8, 40), (16, 40), (32, 30), (64, 24), (128, 16), (256, 6)]:
        res = []
        t0 = time.time()
        for _ in range(trials):
            r = one_trial(sigma, rng, cap, grid)
            res.append(r if r is not None else cap)
        res = np.array(res)
        ln_X = n * math.log(3.0)          # success = recovering s (3^n candidates)
        row = dict(sigma=sigma, trials=trials, N_median=int(np.median(res)), N_min=int(res.min()),
                   N_max=int(res.max()),
                   N_floor=(round(N_floor(ln_X, n, V_ETA, sigma, 0.5), 2) if sigma > 0 else 0),
                   N_kl1=(round(2 * sigma ** 2 / B2(V_ETA, n, 2.0 ** -40) ** 2, 4) if sigma > 0 else 0))
        rows.append(row)
        print(row, "time %.1fs" % (time.time() - t0), flush=True)
    with open(os.path.join(DATA, "chart_attack.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("V_ETA", V_ETA, "B2", B2(V_ETA, n, 2.0 ** -40))


if __name__ == "__main__":
    main()
