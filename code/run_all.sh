#!/bin/sh
# Re-creates every CSV in ../data (timing is machine dependent) and the plots.
set -e
cd "$(dirname "$0")"
python3 exp1_budget.py
python3 exp2_attack.py
python3 exp3_correct.py
python3 exp4_timing.py
python3 check_results.py
python3 plot_results.py
