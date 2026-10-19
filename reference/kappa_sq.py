"""Certified kappa as in kappa_only.py, but with the substitution t = s^2:

    kappa = e^{v_min} * int_0^inf P(s^2) e^{-s^2} 2 s ds.

P contains sqrt(t), a branch point at t = 0 where adaptive ball integration can exhaust its
evaluation limit for very small barriers. After t = s^2 the integrand is smooth on the real
axis. Same tail bound. Used to recompute grid points that failed or were not sharp enough
(repair.py).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flint
from flint import arb, acb
flint.ctx.prec = 333
from ptilde import Ptilde
NAN = acb("nan")


def kappa_sq(v1, v2, u, margin=50.0):
    v1, v2, u = arb(v1), arb(v2), arb(u)
    vmin = v1 if float(v1) <= float(v2) else v2
    vmax = v2 if float(v1) <= float(v2) else v1
    dlt = vmax - vmin
    T0 = float(vmin.upper()) + margin
    S0 = T0 ** 0.5

    def f(z, analytic):
        if analytic and not (z.real.lower() > 0):
            return NAN
        try:
            t = z * z
            out = Ptilde(t, vmin, dlt, acb(u)) * (-t).exp() * 2 * z
            return out if out.is_finite() else NAN
        except Exception:
            return NAN
    cuts = [0.0] + [c for c in (1e-3, 1e-2, 0.03, 0.1, 0.3, 1.0, 2.0, 4.0) if c < S0] + [S0]
    I = acb(0)
    for a, b in zip(cuts[:-1], cuts[1:]):
        J = acb.integral(f, acb(a), acb(b), eval_limit=2000000, depth_limit=400)
        if not J.real.is_finite():
            return None
        I += J
    tail = arb(T0).neg().exp()
    val = arb(I.real.mid(), I.real.rad()) + arb(0, float(tail)) + tail / 2
    return vmin.exp() * val


if __name__ == "__main__":
    import math, time
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from kappa_only import kappa
    for a1, a2, u in [(4, 8, 6), (0.5, 0.5, 1), (20, 20, 10), (1, 4, 0.05), (4, 4, 0.1)]:
        v1 = arb(a1) * u / (2 * arb.pi()); v2 = arb(a2) * u / (2 * arb.pi())
        t = time.time(); k2 = kappa_sq(v1, v2, arb(u)); t2 = time.time() - t
        t = time.time(); k1 = kappa(v1, v2, arb(u)); t1 = time.time() - t
        print(a1, a2, u, "sq:", k2.str(15) if k2 is not None else None, f"{t2:.1f}s",
              " alt:", k1.str(15) if k1 is not None else None, f"{t1:.1f}s",
              " ueberlappen:", (k1.overlaps(k2) if (k1 is not None and k2 is not None) else "-"), flush=True)
