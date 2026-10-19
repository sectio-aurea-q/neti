r"""
Zertifiziertes kappa allein -- ohne Ableitungen, also ohne Potenzreihen.

Gebraucht fuer den Pruefsatz gegen Arkane: dort ist u = h nu/(k_B T) oft weit unter 1
(hohe Temperatur, niedrige Sattelfrequenz), und in diesem Bereich bricht die
Reihenrechnung fuer die Ableitungen zusammen. kappa selbst braucht sie nicht.

Schwanz: fuer t >= T0 gilt |1 - Ptilde| <= 1 (Ptilde ist eine Wahrscheinlichkeit,
0 <= Ptilde <= 1 fuer reelles t), also

    int_T0^inf Ptilde e^{-t} dt  in  [0, e^{-T0}]

Das ist grob, aber bei T0 = 60 ist e^{-60} = 8,8e-27 und damit weit unter jeder
hier interessierenden Genauigkeit. Das Lemma aus cert2.py liefert es schaerfer,
braucht aber eine Scheibe in u; hier genuegt die elementare Schranke.
"""
import os, sys, flint
from flint import arb, acb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
flint.ctx.prec = 333
from ptilde import Ptilde
NAN = acb("nan")

def kappa(v1, v2, u, margin=50.0):
    """kappa = e^{v_min} int_0^inf Ptilde e^{-t} dt.

    T0 MUSS mit v_min skalieren. Der Schwanz ist int_T0^inf Ptilde e^{-t} dt in
    [0, e^{-T0}], also der absolute Fehler in kappa hoechstens e^{v_min - T0}.
    Ein festes T0 = 60 ergibt bei v_min = 82,9 (Gasphasenbarriere bei niedriger
    Temperatur) einen Fehler von e^{22,9} ~ 1e10 -- und die Masse des Integranden
    sitzt dort ohnehin erst bei t ~ 56. Mit T0 = v_min + 50 ist der Fehler
    e^{-50} = 2e-22.
    """
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
