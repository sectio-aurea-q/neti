"""Certified kappa (Arb ball arithmetic) by direct integration in t = energy/k_B T.

Tail: P <= 1 for real t, so int_T0^inf P e^{-t} dt lies in [0, e^{-T0}]; with
T0 = v_min + 50 the absolute error in kappa is at most e^{-50}.
Used for the first pass of the reference grid (refgrid.py).
"""
import os, sys, flint
from flint import arb, acb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
flint.ctx.prec = 333
from ptilde import Ptilde
NAN = acb("nan")

def kappa(v1, v2, u, margin=50.0):
    """kappa = e^{v_min} int_0^inf P e^{-t} dt; T0 = v_min + margin, so the tail error is at most e^{-margin}."""
    v1, v2, u = arb(v1), arb(v2), arb(u)
    vmin = v1 if float(v1) <= float(v2) else v2
    vmax = v2 if float(v1) <= float(v2) else v1
    dlt = vmax - vmin
    T0 = float(vmin.upper()) + margin
    def f(z, analytic):
        if analytic and not (z.real.lower() > 0):
            return NAN
        try:
            out = Ptilde(z, vmin, dlt, acb(u)) * (-z).exp()
            return out if out.is_finite() else NAN
        except Exception:
            return NAN
    I = acb.integral(f, acb(0), acb(T0), eval_limit=400000, depth_limit=160)
    if not I.real.is_finite():
        return None
    tail = arb(T0).neg().exp()
    val = arb(I.real.mid(), I.real.rad()) + arb(0, float(tail)) + tail/2
    return vmin.exp() * val
