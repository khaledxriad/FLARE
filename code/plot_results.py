"""plot_results.py -- quick matplotlib plots of the CSV results (PNG, written to ../plots)."""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(HERE, "..", "plots")
os.makedirs(OUT, exist_ok=True)


def rows(name):
    with open(os.path.join(DATA, name)) as f:
        return list(csv.DictReader(f))


def save(fig, name):
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, name), dpi=160)
    plt.close(fig)


r = rows("chart_sigma_vs_lambda.csv")
fig, ax = plt.subplots(figsize=(5, 3.4))
for key, lab in [("ours", "Renyi (this work)"), ("dp_adv", "Gaussian + advanced composition"),
                 ("gauss_sd", "Gaussian + Pinsker"), ("uniform", "Uniform smudging")]:
    ax.plot([float(x["lam"]) for x in r], [float(x[key]) for x in r], marker="o", ms=3, label=lab)
ax.set_xlabel("target security (bits)"); ax.set_ylabel("log2 of required sigma"); ax.legend(fontsize=7)
save(fig, "sigma_vs_lambda.png")

r = rows("chart_modulus_vs_N.csv")
fig, ax = plt.subplots(figsize=(5, 3.4))
for key in ["ours", "dp_adv", "gauss_sd", "uniform"]:
    ax.plot([float(x["log2N"]) for x in r], [float(x[key]) for x in r], marker="o", ms=3, label=key)
ax.set_xlabel("log2 N (leaked decodings)"); ax.set_ylabel("modulus bits needed"); ax.legend(fontsize=7)
save(fig, "modulus_vs_N.png")

r = rows("chart_attack.csv")
fig, ax = plt.subplots(figsize=(5, 3.4))
s = [float(x["sigma"]) for x in r]
ax.loglog([v if v > 0 else 1 for v in s], [float(x["N_median"]) for x in r], marker="o", label="median queries")
ax.loglog([v if v > 0 else 1 for v in s], [max(float(x["N_floor"]), 0.5) for x in r], marker="s", label="information floor")
ax.set_xlabel("sigma"); ax.set_ylabel("queries N"); ax.legend(fontsize=7)
save(fig, "attack.png")

r = rows("chart_timing.csv")
fig, ax = plt.subplots(figsize=(5, 3.4))
n = [x["n"] for x in r]
for key in ["plain_us", "PCG64_us", "Philox_us", "urandom_us"]:
    ax.plot(n, [float(x[key]) for x in r], marker="o", label=key)
ax.set_xlabel("ring degree n"); ax.set_ylabel("time per decode (us)"); ax.legend(fontsize=7)
save(fig, "timing.png")
print("plots written to", OUT)
