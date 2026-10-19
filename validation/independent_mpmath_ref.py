"""Independent Eckart kappa from the Eckart potential itself (Landau-Lifshitz form),
hbar = m = k_B T = 1. No J&H alpha parametrisation is used.
V(x) = A y/(1+y) + B y/(1+y)^2, y = exp(x/L); A = V1 - V2, B = (sqrt V1 + sqrt V2)^2.
(Closed form verified against direct numerical Schroedinger integration in ode_check.py.)
P(E) = [cosh(2pi L (k1+k2)) - cosh(2pi L (k1-k2))] / [cosh(2pi L (k1+k2)) + cosh(pi sqrt(8 B L^2 - 1))]
k1 = sqrt(2E), k2 = sqrt(2(E-A)), E measured from the left (V1-side) asymptote.
L fixed by |V''(x_max)| = omega^2 = u^2 (hbar omega = h nu = u k_B T)."""
import mpmath as mp

def setup(v1, v2, u):
    v1, v2, u = mp.mpf(v1), mp.mpf(v2), mp.mpf(u)
    A = v1 - v2
    B = (mp.sqrt(v1) + mp.sqrt(v2))**2
    # curvature at max for L=1, numerically (independent of closed forms)
    V = lambda x: A*mp.exp(x)/(1+mp.exp(x)) + B*mp.exp(x)/(1+mp.exp(x))**2
    ystar = (B + A)/(B - A)
    xs = mp.log(ystar)
    assert abs(V(xs) - v1) < mp.mpf(10)**(-mp.mp.dps+10)*max(1, v1), (V(xs), v1)
    c2 = -mp.diff(V, xs, 2)
    L = mp.sqrt(c2)/u
    tpd = mp.pi*mp.sqrt(mp.mpc(8*B*L**2 - 1))
    D = mp.re(mp.cosh(tpd))
    return A, B, L, D

def P_of_E(E, A, L, D):
    k1 = mp.sqrt(2*E); k2 = mp.sqrt(2*(E - A))
    # numerator via exact identity cosh(x+y)-cosh(x-y) = 2 sinh x sinh y
    num = 2*mp.sinh(2*mp.pi*L*k1)*mp.sinh(2*mp.pi*L*k2)
    return num/(mp.cosh(2*mp.pi*L*(k1+k2)) + D)

def kappa(v1, v2, u, dps=50, tailmargin=150):
    with mp.workdps(dps + 20):
        v1m, v2m, um = mp.mpf(v1), mp.mpf(v2), mp.mpf(u)
        A, B, L, D = setup(v1m, v2m, um)
        E0 = max(mp.mpf(0), A)          # higher asymptote
        vmin = min(v1m, v2m)
        f = lambda t: 2*t*P_of_E(E0 + t*t, A, L, D)*mp.exp(vmin - t*t)
        w = um/(2*mp.pi)                # sigmoid width at barrier top
        eps = set()
        for j in (-40, -20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20, 40):
            eps.add(vmin + j*w)
        for c in (0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 1.1):
            eps.add(c*vmin)
        for c in (1e-6, 1e-4, 1e-2, 0.1, 0.5, 1, 2, 5, 12, 30):
            eps.add(mp.mpf(c))
        for c in (2, 5, 15, 40):
            eps.add(vmin + c)
        emax = vmin + tailmargin
        pts = sorted(set([mp.mpf(0)] + [e for e in eps if 0 < e < emax] + [emax]))
        ts = [mp.sqrt(e) for e in pts]
        val, err = mp.quad(f, ts, error=True, maxdegree=10)
        return +val, err

if __name__ == "__main__":
    import sys
    for a in [(4, 4, 2), (1e-3, 1e-3, 1.0), (0.05, 200, 8), (40, 45, 60), (0.5, 0.7, 0.01), (1, 1, 3),
              (0.9999999, 0.9999999, 4.0), (1.0000001, 1.0000001, 4.0)]:
        print(a, mp.nstr(kappa(*a)[0], 20), mp.nstr(kappa(*a)[1], 3))
