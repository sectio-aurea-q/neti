"""What the conformance test reports for the common textbook arrangement.

P = 1 - (cosh(A-B) + cosh(D)) / (cosh(A+B) + cosh(D)),  D = 2*sqrt(|alpha1*alpha2 - pi^2/4|)
integrated with the same accurate quadrature (scipy), so only the arithmetic of P differs.
"""
import math, os, sys
from scipy.integrate import quad
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "conformance"))
from conformance import run

def naive_kappa(v1, v2, u):
    lo, hi = min(v1, v2), max(v1, v2)
    k = 2 * math.pi / u
    a_lo, a_hi = k * lo, k * hi
    s = 1 / math.sqrt(a_lo) + 1 / math.sqrt(a_hi)
    D = 2 * math.sqrt(abs(a_lo * a_hi - math.pi ** 2 / 4))
    def P(e):
        A = 2 * math.sqrt(k * e) / s
        B = 2 * math.sqrt(a_hi - a_lo + k * e) / s
        try:
            return 1 - (math.cosh(A - B) + math.cosh(D)) / (math.cosh(A + B) + math.cosh(D))
        except OverflowError:
            return float("nan")
    f = lambda t: 2 * t * P(t * t) * math.exp(lo - t * t)
    tmax = math.sqrt(lo + 80)
    pts = sorted({math.sqrt(x) for x in (0.25*lo, 0.5*lo, lo, lo + 2, lo + 8) if 0 < x < lo + 80})
    return quad(f, 0, tmax, points=pts, limit=400, epsabs=0, epsrel=1e-12)[0]

if __name__ == "__main__":
    run(naive_kappa, tol=1e-6)
