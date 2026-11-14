"""check_results.py -- recompute headline numbers from data/*.csv and print PASS/FAIL."""
import csv, os, sys
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
rd = lambda n: list(csv.DictReader(open(os.path.join(D, n))))
ok = True
def chk(name, cond):
    global ok
    ok &= bool(cond); print(("PASS " if cond else "FAIL ") + name)
att = {int(r["sigma"]): r for r in rd("chart_attack.csv")}
chk("attack: 3 queries without flooding", int(att[0]["N_median"]) == 3)
chk("attack: median above floor for all sigma>0", all(float(r["N_median"]) > float(r["N_floor"]) for r in att.values() if int(r["sigma"]) > 0))
des = {r["name"]: r for r in rd("table_design.csv")}
chk("design: Kyber-like sets infeasible at q=3329", all(des[k]["feasible"] == "0" for k in des if k.startswith("Kyber")))
chk("design: RLWE-2048 and RLWE-4096 feasible", all(des[k]["feasible"] == "1" for k in des if k.startswith("RLWE-2048") or k.startswith("RLWE-4096")))
par = rd("table_params.csv")
chk("params: no Monte Carlo norm exceeds the bound", all(int(r["exceed"]) == 0 for r in par))
sys.exit(0 if ok else 1)
