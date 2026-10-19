"""Compare eckart-ref with the certified grid (conformance/eckart_certified.csv)."""
import csv, os, statistics, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
rows = list(csv.DictReader(open(os.path.join(ROOT, "conformance", "eckart_certified.csv"))))
inp = os.path.join(HERE, "grid_in.csv")
with open(inp, "w") as g:
    for r in rows: g.write(f"{r['v1']},{r['v2']},{r['ustar']}\n")
out = subprocess.run([os.path.join(ROOT, "target", "release", "eckart-ref"), "grid", inp], capture_output=True, text=True).stdout.strip().split("\n")
res, errs = [], 0
for r, line in zip(rows, out):
    p = line.split(",")
    if p[3].startswith("ERR"): errs += 1; continue
    k = float(p[3]); mid = float(r["kappa_mid"])
    res.append((abs(k / mid - 1), r["domain"], r["alpha1"], r["alpha2"], r["ustar"][:8], r["cos_branch"]))
res.sort(reverse=True)
print(f"points {len(rows)}  errors {errs}")
print(f"max rel err {res[0][0]:.3e} at {res[0][1:]}")
print(f"median rel err {statistics.median(x[0] for x in res):.2e}")
for b in ("0", "1"):
    s = [x[0] for x in res if x[5] == b]; print(f"cos_branch={b}: {len(s)} points, max {max(s):.2e}")
