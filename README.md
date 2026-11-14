# FLARE: code and data

Reference code and data for **FLARE: Flooding for Leakage-Resilient Lattice Encryption**
(FLooding with Renyi-Exponent accounting).

Author: Khaled Riad, Computer Science Department, College of Computer Sciences & Information
Technology, King Faisal University, Al-Ahsa, 31982, Saudi Arabia (kriad@kfu.edu.sa).

## Contents

| Path | What it is |
|---|---|
| `code/floodlib.py` | Closed-form tools: parameter sets, Renyi/KL formulas, square-root law, sigma and N_max, baselines |
| `code/exp1_budget.py` | Parameter tables, Monte Carlo exposure statistics, budget and baseline curves |
| `code/exp2_attack.py` | Passive least-squares key-recovery attack on a toy ring-LWE instance (n = 64) |
| `code/exp3_correct.py` | End-to-end flooded decryption and correctness check |
| `code/exp4_timing.py` | Timing of flooding with three random generators |
| `code/check_results.py` | Recomputes headline claims from the CSV files (PASS/FAIL) |
| `code/plot_results.py` | Matplotlib plots of the CSV data (written to `plots/`) |
| `data/*.csv` | Every number reported in the paper's tables and charts |

## Quick start

```
pip install -r requirements.txt
sh code/run_all.sh        # re-runs all experiments, checks the claims, writes plots
```

All experiments use fixed seeds, so `chart_*.csv` (except `chart_timing.csv`) and
`table_*.csv` are reproduced exactly. Timing depends on the machine.
Runtime is a few minutes on a laptop.

## Scope and limits

The experiments use toy sizes and plain Python/NumPy. They validate the noise model, the
correctness bound and the direction of the leakage bounds. They are not a constant-time or
production implementation, and they do not measure physical side-channel traces.

## Citation

See `CITATION.cff`. Archived on Zenodo: [Zendo url]. Repository: [GitHub Link].

## License

MIT, see `LICENSE`.
