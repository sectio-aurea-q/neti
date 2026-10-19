"""Build conformance/eckart_certified.csv from the raw Arb reference grid in reference/.

* Rows of reference/eckart_reference.csv that failed or have relative radius > 1e-10 are
  replaced by the recomputed rows in reference/eckart_reference_repair.csv.
* refgrid.py passed alpha and u* to Arb as IEEE doubles, so each enclosure is for the exact
  double value. The `ustar` column therefore holds the exact decimal expansion of that double,
  and v1, v2 are computed from it (40 digits).
* The raw lo/hi were printed with 25 significant digits by round-to-nearest. They are widened
  here by a relative 1e-24 so that [kappa_lo, kappa_hi] remains a valid enclosure.
"""
import csv, os
from decimal import Decimal
import mpmath
mpmath.mp.dps = 50
HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "..", "reference")
rows = {}
for name in ("eckart_reference.csv", "eckart_reference_repair.csv"):
    for r in csv.DictReader(open(os.path.join(REF, name))):
        if r["kappa_mid"] in ("FAIL", "") or not r["rel_radius"]:
            continue
        if name == "eckart_reference.csv" and float(r["rel_radius"]) > 1e-10:
            continue
        rows[(r["alpha1"], r["alpha2"], r["ustar"])] = r
widen = mpmath.mpf("1e-24")
with open(os.path.join(HERE, "eckart_certified.csv"), "w", newline="") as g:
    w = csv.writer(g)
    w.writerow(["domain", "alpha1", "alpha2", "ustar", "v1", "v2", "cos_branch", "kappa_lo", "kappa_hi", "kappa_mid"])
    for key in sorted(rows, key=lambda k: (rows[k]["domain"], float(k[0]), float(k[1]), float(k[2]))):
        r = rows[key]
        a1 = mpmath.mpf(float(key[0])); a2 = mpmath.mpf(float(key[1]))
        u = mpmath.mpf(float(key[2]))
        v1 = a1 * u / (2 * mpmath.pi); v2 = a2 * u / (2 * mpmath.pi)
        lo = mpmath.mpf(r["kappa_lo"]) * (1 - widen); hi = mpmath.mpf(r["kappa_hi"]) * (1 + widen)
        w.writerow([r["domain"], key[0], key[1], str(Decimal(float(key[2]))), mpmath.nstr(v1, 40), mpmath.nstr(v2, 40),
                    r["cos_branch"], mpmath.nstr(lo, 28), mpmath.nstr(hi, 28), r["kappa_mid"]])
print(len(rows), "rows written")
