"""Conformance test for any Eckart implementation against certified enclosures.

Each row of eckart_certified.csv is a ball-arithmetic (Arb) enclosure [kappa_lo, kappa_hi]
of the exact Eckart tunnelling factor (Johnston & Heicklen 1962 / Brown 1981 definition,
energies from the higher asymptote, kappa = e^{v_min} * int_0^inf P(e) e^{-e} de).

Usage as a library:

    from conformance import run
    def my_kappa(v1, v2, u):     # barriers in k_B T, u = h nu / k_B T
        ...
    run(my_kappa, tol=1e-6)

Usage from the shell (tests eckart-ref itself):

    python3 conformance.py
"""
import csv, math, os

HERE = os.path.dirname(os.path.abspath(__file__))

def load():
    with open(os.path.join(HERE, "eckart_certified.csv")) as f:
        return list(csv.DictReader(f))

def run(fn, tol=1e-6, show=10):
    rows = load()
    bad, worst, failed = [], 0.0, 0
    for r in rows:
        v1, v2, u = float(r["v1"]), float(r["v2"]), float(r["ustar"])
        ref = float(r["kappa_mid"])
        try:
            k = float(fn(v1, v2, u))
        except Exception as e:
            failed += 1; bad.append((math.inf, r, f"exception: {e}")); continue
        rel = abs(k / ref - 1) if math.isfinite(k) else math.inf
        worst = max(worst, rel)
        if rel > tol:
            bad.append((rel, r, k))
    bad.sort(key=lambda x: -x[0])
    print(f"{len(rows)} certified points, tolerance {tol:g}")
    n_inf = sum(1 for b in bad if b[0] == math.inf)
    print(f"  pass: {len(rows) - len(bad)}   fail: {len(bad)}   (of which exception/NaN/inf: {n_inf})")
    print(f"  worst relative error: {worst:.3e}")
    for rel, r, k in bad[:show]:
        print(f"  alpha1={r['alpha1']:>5} alpha2={r['alpha2']:>5} u*={r['ustar']:>5} cos={r['cos_branch']}  "
              f"certified={float(r['kappa_mid']):.10g}  yours={k}  rel={rel:.2e}")
    return len(bad) == 0

if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.join(HERE, "..", "python"))
    import eckart_ref
    ok = run(eckart_ref.kappa, tol=1e-12)
    raise SystemExit(0 if ok else 1)
