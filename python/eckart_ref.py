"""Python interface to eckart-ref (ctypes; no build tools needed beyond `cargo build --release`).

    import eckart_ref as er
    er.kappa(v1, v2, u)                 # reduced units: barriers in k_B T, u = h nu / k_B T
    er.kappa_si(V1, V2, nu_cm, T)       # J/mol, J/mol, cm^-1, K (CODATA 2018)
    er.transmission(v1, v2, u, eps)     # P(eps), eps above the higher asymptote in k_B T
"""
import ctypes, os, sys

def _load():
    here = os.path.dirname(os.path.abspath(__file__))
    name = {"darwin": "libeckart_ref.dylib", "win32": "eckart_ref.dll"}.get(sys.platform, "libeckart_ref.so")
    for d in (os.environ.get("ECKART_REF_LIB", ""), os.path.join(here, "..", "target", "release"), here):
        p = os.path.join(d, name)
        if os.path.exists(p):
            return ctypes.CDLL(p)
    raise OSError(f"{name} not found; run `cargo build --release` or set ECKART_REF_LIB")

_lib = _load()
_d = ctypes.c_double
_lib.eckart_kappa.argtypes = [_d, _d, _d, ctypes.POINTER(_d), ctypes.POINTER(_d)]
_lib.eckart_kappa.restype = ctypes.c_int
_lib.eckart_kappa_si.argtypes = [_d, _d, _d, _d, ctypes.POINTER(_d), ctypes.POINTER(_d)]
_lib.eckart_kappa_si.restype = ctypes.c_int
_lib.eckart_transmission.argtypes = [_d, _d, _d, _d]
_lib.eckart_transmission.restype = _d

def _call(fn, *args):
    v, e = _d(), _d()
    rc = fn(*args, ctypes.byref(v), ctypes.byref(e))
    if rc == 1: raise ValueError("invalid input (all arguments must be finite and > 0)")
    if rc == 2: raise RuntimeError("quadrature did not converge")
    return v.value

def kappa(v1, v2, u):
    return _call(_lib.eckart_kappa, v1, v2, u)

def kappa_si(v1_j_mol, v2_j_mol, nu_cm, t_k):
    return _call(_lib.eckart_kappa_si, v1_j_mol, v2_j_mol, nu_cm, t_k)

def transmission(v1, v2, u, eps):
    return _lib.eckart_transmission(v1, v2, u, eps)
