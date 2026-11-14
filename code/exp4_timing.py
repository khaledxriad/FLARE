"""exp4_timing.py -- cost of flooded decoding in the Python/NumPy reference
implementation.  q = 2^27 for all sizes (word size, not security, is what the
timing depends on).  Median of repeated runs, single thread."""
import csv, os, sys, time, platform
import numpy as np
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from floodlib import negacyclic_mul_float

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
q = 2 ** 27
SIGMA = 2.0 ** 23            # flooding std used for the timing (inside sigma_max for n>=1024)


def bm_gauss_urandom(sigma, n):
    """Cryptographic-source Gaussian sampler: os.urandom + Box-Muller."""
    u = np.frombuffer(os.urandom(16 * (n // 2 + 1)), dtype=np.uint64).astype(np.float64) / 2.0 ** 64
    u1, u2 = u[: n // 2 + 1], u[n // 2 + 1: 2 * (n // 2 + 1)]
    u1 = np.maximum(u1, 1e-300)
    R = np.sqrt(-2.0 * np.log(u1))
    z = np.concatenate([R * np.cos(2 * np.pi * u2), R * np.sin(2 * np.pi * u2)])[:n]
    return sigma * z


def timeit(fn, reps):
    ts = []
    for _ in range(reps):
        t = time.perf_counter(); fn(); ts.append(time.perf_counter() - t)
    return float(np.median(ts)) * 1e6        # microseconds


def timeit_interleaved(fns, reps):
    """Interleave the variants so that drift affects all of them equally."""
    ts = [[] for _ in fns]
    for _ in range(reps):
        for i, fn in enumerate(fns):
            t = time.perf_counter(); fn(); ts[i].append(time.perf_counter() - t)
    return [float(np.median(x)) * 1e6 for x in ts]


def main():
    rows = []
    rngs = {"PCG64": np.random.Generator(np.random.PCG64(1)),
            "Philox": np.random.Generator(np.random.Philox(2))}
    for n in (1024, 2048, 4096):
        base = np.random.default_rng(5)
        s = base.integers(-1, 2, n); c1 = base.integers(0, q, n); c2 = base.integers(0, q, n)

        def dec_plain():
            w = (c2 - negacyclic_mul_float(c1, s)) % q
            c = np.where(w > q // 2, w - q, w)
            return (np.abs(c) > q // 4)

        fl = {}
        for name, g in rngs.items():
            def dec_flood(g=g):
                w = (c2 - negacyclic_mul_float(c1, s)
                     + np.rint(g.standard_normal(n) * SIGMA).astype(np.int64)) % q
                c = np.where(w > q // 2, w - q, w)
                return (np.abs(c) > q // 4)
            fl[name] = dec_flood

        def dec_flood_u():
            w = (c2 - negacyclic_mul_float(c1, s) + np.rint(bm_gauss_urandom(SIGMA, n)).astype(np.int64)) % q
            c = np.where(w > q // 2, w - q, w)
            return (np.abs(c) > q // 4)
        names = ["plain"] + list(fl.keys()) + ["urandom"]
        fns = [dec_plain] + list(fl.values()) + [dec_flood_u]
        for _ in range(100):           # warm-up
            for f in fns: f()
        med = timeit_interleaved(fns, 3000)
        row = dict(n=n)
        for nm, v in zip(names, med):
            row[nm + "_us"] = round(v, 1)
        rows.append(row); print(row, flush=True)
    with open(os.path.join(DATA, "chart_timing.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(platform.platform(), platform.python_version(), np.__version__)


if __name__ == "__main__":
    main()
