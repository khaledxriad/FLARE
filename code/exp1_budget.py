"""exp1_budget.py -- parameter tables, Monte-Carlo validation of the noise
model, leakage budgets and baseline comparison.  Writes CSV files to ../data."""
import csv, math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from floodlib import *

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
rng = np.random.default_rng(20261114)


def sample_cbd(eta, size):
    return (rng.integers(0, 2, size + (eta,)).sum(-1) - rng.integers(0, 2, size + (eta,)).sum(-1))


def sample_dist(var, size):
    """Sample a small-coefficient vector with the requested variance
    (CBD / ternary / rounded Gaussian, matching the parameter definitions)."""
    if abs(var - TERN) < 1e-9:
        return rng.integers(-1, 2, size)
    if abs(var - G32) < 1e-9:
        return np.rint(rng.normal(0, 3.2, size)).astype(np.int64)
    eta = int(round(var * 2))
    return sample_cbd(eta, size)


def mc_eta(par, trials):
    n, k = par.n, par.k
    out = np.empty((trials, n))
    for t in range(trials):
        s = sample_dist(par.vs, (k, n)); e = sample_dist(par.ve, (k, n))
        r = sample_dist(par.vr, (k, n)); e1 = sample_dist(par.v1, (k, n))
        e2 = sample_dist(par.v2, (n,))
        eta = e2.astype(np.int64).copy()
        for i in range(k):
            eta += negacyclic_mul_float(e[i], r[i])
            eta -= negacyclic_mul_float(s[i], e1[i])
        if par.du is not None:             # compression error terms
            hu = par.q / 2 ** (par.du + 1)
            du = np.rint(rng.uniform(-hu, hu, (k, n)))
            hv = par.q / 2 ** (par.dv + 1)
            dv = np.rint(rng.uniform(-hv, hv, n))
            eta = eta + dv
            for i in range(k):
                eta += negacyclic_mul_float(s[i], du[i].astype(np.int64))
        out[t] = eta
    return out


def main():
    N0 = 2 ** 20
    lam0, lam_t = 128, 112
    # ---- Table 1: parameters + Monte-Carlo validation ---------------------
    rows = []
    for par in PARAMS:
        trials = 400 if par.n <= 256 else (150 if par.n == 1024 else 40)
        eta = mc_eta(par, trials)
        var_mc = float((eta ** 2).mean())
        l2 = np.sqrt((eta ** 2).sum(1))
        bound = B2(par.v_eta, par.ell, 2.0 ** -40)
        rows.append(dict(name=par.name, n=par.n, k=par.k, log2q=round(par.log2q, 1), ell=par.ell,
                         var_analytic=round(par.v_eta, 2), var_mc=round(var_mc, 2),
                         std_analytic=round(math.sqrt(par.v_eta), 2),
                         l2_mean_mc=round(float(l2.mean()), 1), l2_max_mc=round(float(l2.max()), 1),
                         B2_bound=round(bound, 1), trials=trials,
                         exceed=int((l2 > bound).sum())))
        print(rows[-1], flush=True)
    write_csv("table_params.csv", rows)

    # ---- Table 2: flooding design at N=2^20 -------------------------------
    rows = []
    for par in PARAMS:
        sneed = sigma_needed(par, N0, lam0, lam_t)
        smax = sigma_max_correct(par, 2.0 ** -40)
        nmax = N_max(par, smax, lam0, lam_t) if smax > 0 else 0.0
        feasible = sneed <= smax
        rows.append(dict(name=par.name, sigma_needed=round(sneed, 1), log2_sigma_needed=round(math.log2(sneed), 2),
                         sigma_max=round(smax, 1), feasible=int(feasible),
                         log2_Nmax_at_sigmamax=(round(math.log2(nmax), 2) if nmax > 0 else "none"),
                         log2q_needed=round(math.log2(q_required(par, sneed, 2.0 ** -40)), 1),
                         log2q_have=round(par.log2q, 1)))
        print(rows[-1], flush=True)
    write_csv("table_design.csv", rows)

    # ---- Chart: N_max vs sigma (RLWE sets) --------------------------------
    rows = []
    for par in PARAMS[3:]:
        smax = sigma_max_correct(par, 2.0 ** -40)
        grid = [x for x in range(14, 41, 2) if x < math.log2(smax)] + ([math.log2(smax)] if math.log2(smax) <= 41 else [])
        for x in grid:
            nm = N_max(par, 2.0 ** x, lam0, lam_t)
            rows.append(dict(name=par.name, log2_sigma=round(x, 3), log2_Nmax=(round(math.log2(nm), 3) if nm > 0 else "")))
    write_csv("chart_budget.csv", rows)

    # ---- Chart: required sigma vs lambda, 4 methods -----------------------
    par = PARAMS[3]
    rows = []
    for lam in [20, 40, 60, 80, 100, 128]:
        B = B_for(par, N0, lam)
        ours = sigma_needed(par, N0, lam0=lam, lam_t=lam - 16)
        rows.append(dict(lam=lam, ours=round(math.log2(ours), 3),
                         dp_adv=round(math.log2(sigma_dp_advanced(B, N0, lam)), 3),
                         gauss_sd=round(math.log2(sigma_gauss_pinsker(B, N0, lam)), 3),
                         uniform=round(math.log2(sigma_uniform_smudging(B, par.ell, N0, lam)), 3)))
    write_csv("chart_sigma_vs_lambda.csv", rows)

    # ---- Chart: required modulus vs N (lambda = 128) ----------------------
    rows = []
    lam = 128
    for lg in range(0, 41, 4):
        N = 2 ** lg
        B = B_for(par, N, lam)
        ours = sigma_needed(par, N, lam0=lam, lam_t=lam - 16)
        d = dict(log2N=lg)
        for key, s in (("ours", ours), ("dp_adv", sigma_dp_advanced(B, N, lam)),
                       ("gauss_sd", sigma_gauss_pinsker(B, N, lam)),
                       ("uniform", sigma_uniform_smudging(B, par.ell, N, lam))):
            d[key] = round(math.log2(q_required(par, s, 2.0 ** -40)), 2)
        rows.append(d)
    write_csv("chart_modulus_vs_N.csv", rows)

    # ---- Chart: security bits vs alpha ------------------------------------
    par = PARAMS[4]
    sig = 2.0 ** 24
    rows = []
    for lgN in (20, 25, 28):
        N = 2 ** lgN
        B = B_for(par, N, lam_t)
        for a in [1.25, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512]:
            rows.append(dict(log2N=lgN, alpha=a, bits=round(bits_after_flooding(a, N, B, sig, lam0), 3)))
        print("alpha* for N=2^%d: %.2f  best bits %.2f" % (lgN, best_bits(N, B, sig, lam0)[1], best_bits(N, B, sig, lam0)[0]))
    write_csv("chart_alpha.csv", rows)
    print("sigma used for alpha chart:", sig)


def write_csv(name, rows):
    with open(os.path.join(DATA, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
