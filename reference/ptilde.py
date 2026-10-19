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
    """cosh, einheitlich fuer acb und acb_series (acb_series hat kein cosh)."""
    return (z.exp() + (-z).exp()) / 2


def _c(x):
    """nach acb wandeln, Potenzreihen aber unangetastet lassen.

    Damit ist Ptilde typneutral: Reihenvariable darf u SEIN (Taylor in der
    Frequenz) ODER v_min (Taylor in der Rueckbarriere).  Vorher erzwang eine
    acb()-Wandlung den ersten Fall.
    """
    return x if isinstance(x, acb_series) else acb(x)


def _G(vmin, vmax):
    return 1 / _c(vmin).sqrt() + 1 / _c(vmax).sqrt()


def Ptilde(t, vmin, dlt, w):
    """Durchlasswahrscheinlichkeit in t-Koordinate.

    t            acb, >= 0 (darf 0 enthalten)
    vmin, dlt    kleinere Barriere und Abstand, acb/arb oder acb_series
    w            die Frequenzvariable u, acb/arb oder acb_series

    Genau eine der Groessen vmin/dlt bzw. w traegt die Reihenvariable; die
    t-Abhaengigkeit steht ausserhalb jeder Reihenoperation.
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

