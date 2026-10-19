"""
Certified reference grid for the Eckart tunnelling factor kappa.

Parametrisation as in Johnston & Heicklen (1962) and Brown (1981):
    alpha_i = 2 pi V_i / (h nu),   u* = h nu / (k_B T),   v_i = V_i/(k_B T) = alpha_i u* / (2 pi).
kappa is symmetric in (alpha1, alpha2); rows store alpha1 <= alpha2.
Convention: kappa = exp(v_min) * int_{E >= higher asymptote} P(E) exp(-E/kT) dE/kT.

Domains
  A  Brown's stated range: 0.5 <= alpha <= 20, 0 < u* <= 16,
     u* <= 12 if alpha1 >= 8, u* <= 10 if alpha1 >= 16.
  B  Extension: alpha up to 100, u* up to 30, u* down to 0.05, strong asymmetry.
  Flagged: cos branch (4 alpha1 alpha2 < pi^2, d imaginary).

Each row is an Arb ball enclosure. Note: u* is passed as the IEEE double nearest to the
decimal in the grid, so each enclosure is for that double value (see build_conformance.py).
Resumable: finished points are skipped.
"""
import csv, math, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from kappa_only import kappa

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "eckart_reference.csv")
FIELDS = ["domain", "alpha1", "alpha2", "ustar", "v1", "v2", "cos_branch",
          "kappa_mid", "kappa_lo", "kappa_hi", "rel_radius", "seconds"]

A_AL = [0.5, 1, 2, 4, 8, 12, 16, 20]
A_U = [0.05, 0.1, 0.2, 0.5, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16]
B_AL = [0.5, 1, 2, 4, 8, 16, 25, 40, 60, 80, 100]
B_U = [0.05, 0.2, 0.5, 1, 2, 4, 6, 8, 10, 12, 16, 20, 25, 30]


def in_A(a1, a2, u):
    if not (0.5 <= a1 <= a2 <= 20 and 0 < u <= 16):
        return False
    if a1 >= 16 and u > 10:
        return False
    if a1 >= 8 and u > 12:
        return False
    return True


def points():
    seen = set()
    for al, us in ((A_AL, A_U), (B_AL, B_U)):
        for i, a1 in enumerate(al):
            for a2 in al[i:]:
                for u in us:
                    key = (a1, a2, u)
                    if key in seen:
                        continue
                    seen.add(key)
                    yield ("A" if in_A(a1, a2, u) else "B"), a1, a2, u


def main():
    done = set()
    if os.path.exists(OUT):
        with open(OUT) as f:
            for r in csv.DictReader(f):
                done.add((float(r["alpha1"]), float(r["alpha2"]), float(r["ustar"])))
    new = not os.path.exists(OUT)
    f = open(OUT, "a", newline="")
    w = csv.DictWriter(f, FIELDS)
    if new:
        w.writeheader()
    pts = list(points())
    print(f"{len(pts)} Punkte, davon fertig {len(done)}", flush=True)
    for dom, a1, a2, u in pts:
        if (float(a1), float(a2), float(u)) in done:
            continue
        v1 = flint.arb(a1) * u / (2 * flint.arb.pi())
        v2 = flint.arb(a2) * u / (2 * flint.arb.pi())
        t = time.time()
        k = kappa(v1, v2, flint.arb(u))
        dt = time.time() - t
        row = dict(domain=dom, alpha1=a1, alpha2=a2, ustar=u,
                   v1=v1.str(12), v2=v2.str(12),
                   cos_branch=int(4 * a1 * a2 < math.pi ** 2), seconds=round(dt, 2))
        if k is None or not k.is_finite():
            row.update(kappa_mid="FAIL", kappa_lo="", kappa_hi="", rel_radius="")
        else:
            row.update(kappa_mid=k.mid().str(25, radius=False),
                       kappa_lo=k.lower().str(25, radius=False),
                       kappa_hi=k.upper().str(25, radius=False),
                       rel_radius=f"{float(k.rad() / abs(k.mid())):.2e}")
        w.writerow(row)
        f.flush()
    f.close()
    print("fertig", flush=True)


if __name__ == "__main__":
    main()
