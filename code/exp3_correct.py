"""exp3_correct.py -- end-to-end functional validation of flooded decryption.

Ring-LWE PKE over Z_q[x]/(x^n+1), n=256, q=3329 (ternary secrets/noise).
For each flooding std sigma we encrypt random 256-bit messages, decrypt with
flooded decoding and count decryptions with at least one wrong bit.  The
measured rate is compared with the Theorem-2 prediction."""
import csv, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from floodlib import *

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
n, q = 256, 3329
par = Params("toy-256-3329", n, 1, q, TERN, TERN, TERN, TERN, TERN)


def mulq(a, b):
    return negacyclic_mul_float(a, b) % q


def run(sigma, trials, rng):
    fails = 0
    bit_err = 0
    for _ in range(trials):
        s = rng.integers(-1, 2, n); e = rng.integers(-1, 2, n)
        a = rng.integers(0, q, n)
        b = (mulq(a, s) + e) % q
        r = rng.integers(-1, 2, n); e1 = rng.integers(-1, 2, n); e2 = rng.integers(-1, 2, n)
        m = rng.integers(0, 2, n)
        u = (mulq(a, r) + e1) % q
        v = (mulq(b, r) + e2 + (q // 2) * m) % q
        w = (v - mulq(u, s)) % q                       # = (q//2) m + eta
        w = (w + np.rint(rng.normal(0, sigma, n)).astype(np.int64)) % q   # flooding
        c = np.where(w > q // 2, w - q, w)             # centered
        dec = (np.abs(c) > q // 4).astype(int)         # 1 iff closer to q/2 than to 0
        wrong = int((dec != m).sum())
        bit_err += wrong
        fails += int(wrong > 0)
    return fails / trials, bit_err / (trials * n)


def main():
    rng = np.random.default_rng(31337)
    rows = []
    for sigma in [0, 120, 140, 160, 180, 200, 220, 240, 260]:
        trials = 6000 if sigma >= 160 else 3000
        f, be = run(sigma, trials, rng)
        row = dict(sigma=sigma, trials=trials, fail_measured=round(f, 5), bit_err_measured=round(be, 7),
                   fail_gauss=round(fail_prob_gauss(par, sigma), 5),
                   fail_union=round(fail_prob_union(par, sigma), 5))
        rows.append(row); print(row, flush=True)
    with open(os.path.join(DATA, "chart_correct.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print("v_eta", par.v_eta, "tau", par.tau)


if __name__ == "__main__":
    main()
