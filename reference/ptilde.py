"""Eckart transmission probability in ball arithmetic (Arb via python-flint).

Ptilde(t, vmin, dlt, w): t = energy above the higher asymptote in k_B T, vmin = smaller
barrier, dlt = difference of the barriers, w = u = h nu / k_B T. The imaginary-d case
(small barriers) is handled automatically: sqrt of a negative ball gives an imaginary
ball and cosh(i x) = cos(x).
"""
import sys
import flint
from flint import arb, acb, acb_series

flint.ctx.prec = 333
NAN = acb("nan")


# ----------------------------------------------------------------- Integrand

def _cosh(z):
    """cosh for acb and acb_series alike (acb_series has no cosh)."""
    return (z.exp() + (-z).exp()) / 2


def _c(x):
    """Convert to acb, but leave power series untouched (so Ptilde works for series in u or v_min)."""
    return x if isinstance(x, acb_series) else acb(x)


def _G(vmin, vmax):
    return 1 / _c(vmin).sqrt() + 1 / _c(vmax).sqrt()


def Ptilde(t, vmin, dlt, w):
    """Transmission probability at energy t (in k_B T) above the higher asymptote.

    t            acb, >= 0 (may contain 0)
    vmin, dlt    smaller barrier and barrier difference (acb/arb or acb_series)
    w            u = h nu / k_B T (acb/arb or acb_series)
    """
    pi = acb.pi()
    vmax = vmin + dlt
    G = _G(vmin, vmax)
    inv = 1 / w                       # keine Wurzel aus w -- entscheidend
    beta = 4 * pi * inv / G
    rt0 = t.sqrt()
    rt1 = (t + dlt).sqrt()
    sP = beta * (rt0 + rt1)
    dm = beta * (rt0 - rt1)
    tpd = pi * inv * (16 * _c(vmin) * _c(vmax) - w * w).sqrt()
    return (_cosh(sP) - _cosh(dm)) / (_cosh(sP) + _cosh(tpd))

