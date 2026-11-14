"""floodlib.py -- shared analytic tools for the flooding framework.

All formulas here are the ones proved in the paper (Lemmas 1-4, Theorems 1-3,
Proposition 1).  Nothing in this file talks to the network.
"""
import math
import numpy as np
from scipy.optimize import brentq

LN2 = math.log(2.0)


# --------------------------------------------------------------------------
# Parameter sets (module/ring-LWE style PKE:  eta = e^T r + e2 - s^T e1 + comp)
# --------------------------------------------------------------------------
class Params:
    def __init__(self, name, n, k, q, vs, ve, vr, v1, v2, du=None, dv=None):
        self.name, self.n, self.k, self.q = name, n, k, q
        self.vs, self.ve, self.vr, self.v1, self.v2 = vs, ve, vr, v1, v2
        self.du, self.dv = du, dv
        self.ell = n                      # coefficients exposed by one decoding

    @property
    def log2q(self):
        return math.log2(self.q)

    @property
    def tau(self):                        # decoding margin  floor(q/4) - 1
        return self.q // 4 - 1

    @property
    def v_eta(self):
        nk = self.n * self.k
        v = nk * (self.ve * self.vr + self.vs * self.v1) + self.v2
        if self.du is not None:           # Kyber-style compression error
            v += nk * self.vs * ((self.q / 2 ** self.du) ** 2) / 12.0
            v += ((self.q / 2 ** self.dv) ** 2) / 12.0
        return v


def cbd_var(eta):
    return eta / 2.0


TERN = 2.0 / 3.0         # variance of uniform {-1,0,1}
G32 = 3.2 ** 2           # variance of sigma=3.2 Gaussian noise

PARAMS = [
    Params("Kyber512-like", 256, 2, 3329, cbd_var(3), cbd_var(3), cbd_var(3), cbd_var(2), cbd_var(2), 10, 4),
    Params("Kyber768-like", 256, 3, 3329, cbd_var(2), cbd_var(2), cbd_var(2), cbd_var(2), cbd_var(2), 10, 4),
    Params("Kyber1024-like", 256, 4, 3329, cbd_var(2), cbd_var(2), cbd_var(2), cbd_var(2), cbd_var(2), 11, 5),
    Params("RLWE-1024-27", 1024, 1, 2 ** 27, TERN, G32, TERN, G32, G32),
    Params("RLWE-2048-54", 2048, 1, 2 ** 54, TERN, G32, TERN, G32, G32),
    Params("RLWE-4096-109", 4096, 1, 2 ** 109, TERN, G32, TERN, G32, G32),
]


# --------------------------------------------------------------------------
# Norm bounds (Laurent-Massart chi-square tail; Gaussian heuristic for eta_j)
# --------------------------------------------------------------------------
def B2(v_eta, ell, delta):
    """High-probability l2 bound on the exposure vector eta in R^ell."""
    t = math.log(1.0 / delta)
    return math.sqrt(v_eta * (ell + 2 * math.sqrt(ell * t) + 2 * t))


def Binf(v_eta, ell, delta):
    return math.sqrt(2 * v_eta * math.log(2 * ell / delta))


# --------------------------------------------------------------------------
# Renyi accounting (Lemma 1, Lemma 2, Lemma 3 / Theorem 1)
# --------------------------------------------------------------------------
def c_term(N, B, sigma):
    """c = N B^2 / (2 sigma^2 ln 2): total flooding 'cost' in bits (Theorem 1)."""
    return N * B * B / (2.0 * sigma * sigma * LN2)


def log2_RD_total(alpha, N, B, sigma):
    """log2 of R = E_Q[(p/q)^alpha] for an N-query transcript (Lemmas 1-3)."""
    return alpha * (alpha - 1.0) * c_term(N, B, sigma)


def bits_after_flooding(alpha, N, B, sigma, lam0):
    """Theorem 1:  -log2 Pr_real >= (1-1/alpha)*lam0 - (alpha-1)*c."""
    return (1.0 - 1.0 / alpha) * lam0 - (alpha - 1.0) * c_term(N, B, sigma)


def best_bits(N, B, sigma, lam0):
    """Optimise alpha in closed form: alpha* = sqrt(lam0/c), value (sqrt(lam0)-sqrt(c))^2."""
    c = c_term(N, B, sigma)
    if c >= lam0:
        return 0.0, 1.0
    return (math.sqrt(lam0) - math.sqrt(c)) ** 2, math.sqrt(lam0 / c)


def B_for(par, N, lam_t, extra_bits=1):
    """Norm bound with failure prob. delta_B = 2^-(lam_t + log2 N + extra),
    so that the additive term N*delta_B is <= 2^-(lam_t+1)."""
    d = 2.0 ** (-(lam_t + math.log2(max(N, 1)) + extra_bits))
    return B2(par.v_eta, par.ell, d)


def N_max(par, sigma, lam0=128, lam_t=112):
    """Largest N whose guaranteed security stays >= lam_t (internally lam_t+1,
    the other half of the budget is the additive N*delta_B term)."""
    def ok(N):
        B = B_for(par, N, lam_t)
        return best_bits(N, B, sigma, lam0)[0] >= lam_t + 1
    if not ok(1):
        return 0.0
    lo, hi = 1.0, 2.0
    while ok(hi) and hi < 1e80:
        lo, hi = hi, hi * 2
    for _ in range(100):
        mid = math.sqrt(lo * hi)
        if ok(mid):
            lo = mid
        else:
            hi = mid
    return lo


def sigma_needed(par, N, lam0=128, lam_t=112):
    """Smallest sigma such that N queries retain lam_t bits (closed form)."""
    B = B_for(par, N, lam_t)
    cmax = (math.sqrt(lam0) - math.sqrt(lam_t + 1)) ** 2
    return B * math.sqrt(N / (2.0 * LN2 * cmax))


# --------------------------------------------------------------------------
# Correctness (Theorem 2)
# --------------------------------------------------------------------------
def sigma_max_correct(par, delta_dec):
    """Largest flooding std keeping the decoding failure <= delta_dec
    (Gaussian-sum union bound of Theorem 2)."""
    x = par.tau ** 2 / (2.0 * math.log(2 * par.ell / delta_dec)) - par.v_eta
    return math.sqrt(x) if x > 0 else 0.0


def fail_prob_union(par, sigma):
    """Theorem 2 bound: ell * 2 * exp(-tau^2 / (2 (v_eta + sigma^2)))."""
    return min(1.0, par.ell * 2.0 * math.exp(-par.tau ** 2 / (2.0 * (par.v_eta + sigma ** 2))))


def fail_prob_gauss(par, sigma):
    """Exact Gaussian-tail prediction for any-bit failure of one decoding."""
    from scipy.special import erfc
    s = math.sqrt(par.v_eta + sigma ** 2)
    p = erfc(par.tau / (s * math.sqrt(2)))
    return 1.0 - (1.0 - p) ** par.ell


def q_required(par, sigma, delta_dec):
    """Modulus needed so that flooding with std sigma keeps failure <= delta_dec."""
    s = math.sqrt(par.v_eta + sigma ** 2)
    return 4.0 * math.sqrt(2 * math.log(2 * par.ell / delta_dec)) * s + 4.0


# --------------------------------------------------------------------------
# Baselines (Section "Baselines")
# --------------------------------------------------------------------------
def sigma_uniform_smudging(B2v, ell, N, lam):
    """Classical uniform smudging U[-S,S]^ell, SD <= N*||eta||_1/(2S) <= 2^-lam,
    ||eta||_1 <= sqrt(ell)*||eta||_2.  Returns the std S/sqrt(3)."""
    S = N * math.sqrt(ell) * B2v * 2.0 ** (lam - 1)
    return S / math.sqrt(3.0)


def sigma_gauss_pinsker(B2v, N, lam):
    """Gaussian flooding + Pinsker: SD <= B sqrt(N) / (2 sigma) <= 2^-lam."""
    return B2v * math.sqrt(N) * 2.0 ** (lam - 1)


def sigma_dp_advanced(B2v, N, lam, eps_total=1.0):
    """Gaussian mechanism + advanced composition (Dwork-Rothblum-Vadhan).
    delta budget: N*delta0 + delta' <= 2^-lam; total eps = eps_total."""
    dprime = 2.0 ** (-(lam + 1))
    delta0 = 2.0 ** (-(lam + 1)) / N
    L = math.log(1.0 / dprime)
    g = lambda e0: math.sqrt(2 * N * L) * e0 + N * e0 * (math.exp(e0) - 1.0) - eps_total
    e0 = brentq(g, 1e-12, 1.0)
    return B2v * math.sqrt(2 * math.log(1.25 / delta0)) / e0


# --------------------------------------------------------------------------
# Proposition 1: information-theoretic floor on the number of queries
# --------------------------------------------------------------------------
def N_floor(ln_X, ell, v_eta, sigma, p_succ=0.5):
    per_query = 0.5 * ell * math.log(1.0 + v_eta / sigma ** 2)       # nats
    need = p_succ * ln_X - math.log(2.0)
    return max(need, 0.0) / per_query


# --------------------------------------------------------------------------
# Ring helpers (negacyclic arithmetic) used by the experiments
# --------------------------------------------------------------------------
def negacyclic_matrix(a):
    """Matrix R with R @ x = a * x in Z[x]/(x^n+1)."""
    n = len(a)
    i = np.arange(n)[:, None]
    j = np.arange(n)[None, :]
    idx = (i - j) % n
    sign = np.where(i >= j, 1, -1)
    return a[idx] * sign


def negacyclic_mul_float(a, b):
    """Exact (for small coefficients) negacyclic product via FFT of length 2n."""
    n = len(a)
    fa = np.fft.rfft(a, 2 * n)
    fb = np.fft.rfft(b, 2 * n)
    c = np.fft.irfft(fa * fb, 2 * n)
    c = np.rint(c).astype(np.int64)
    return c[:n] - c[n:]
