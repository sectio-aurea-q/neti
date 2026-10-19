"""Python interface to eckart-ref (ctypes; build the library first with `cargo build --release`).

    import eckart_ref as er
    er.kappa(v1, v2, u)                      # reduced units: barriers in k_B T, u = h nu / k_B T
    er.kappa_si(V1, V2, nu_cm, T)            # J/mol, J/mol, cm^-1 (sign ignored), K
    er.kappa_full(v1, v2, u)                 # dict with value, ln_value, quad_err, tail_bound, ...
    er.log_transmission(v1, v2, u, eps)      # ln P(eps), eps above the higher asymptote in k_B T
    er.transmission_si(E, V1, V2, nu_cm)     # P(E), E from the higher asymptote, all in J/mol
"""
import ctypes, os, sys

def _load():
    here = os.path.dirname(os.path.abspath(__file__))
    name = {"darwin": "libeckart_ref.dylib", "win32": "eckart_ref.dll"}.get(sys.platform, "libeckart_ref.so")
    for d in (os.environ.get("ECKART_REF_LIB", ""), os.path.join(here, "..", "target", "release"), here):
        p = os.path.join(d, name)
        if d and os.path.exists(p):
            return ctypes.CDLL(p)
    raise OSError(f"{name} not found; run `cargo build --release` or set ECKART_REF_LIB")

class _Result(ctypes.Structure):
    _fields_ = [("value", ctypes.c_double), ("ln_value", ctypes.c_double), ("quad_err", ctypes.c_double),
                ("tail_bound", ctypes.c_double), ("cos_branch", ctypes.c_int), ("evals", ctypes.c_uint64)]

_lib = _load()
_d = ctypes.c_double
_lib.eckart_kappa.argtypes = [_d, _d, _d, ctypes.POINTER(_Result)]
_lib.eckart_kappa.restype = ctypes.c_int
_lib.eckart_kappa_si.argtypes = [_d, _d, _d, _d, ctypes.POINTER(_Result)]
_lib.eckart_kappa_si.restype = ctypes.c_int
_lib.eckart_log_transmission.argtypes = [_d, _d, _d, _d]
_lib.eckart_log_transmission.restype = _d
_lib.eckart_transmission_si.argtypes = [_d, _d, _d, _d]
_lib.eckart_transmission_si.restype = _d

def _call(fn, *args):
    r = _Result()
    rc = fn(*args, ctypes.byref(r))
    if rc == 1: raise ValueError("invalid input (arguments must be finite and > 0)")
    if rc == 2: raise RuntimeError("quadrature did not converge")
    return r

def kappa_full(v1, v2, u):
    r = _call(_lib.eckart_kappa, v1, v2, u)
    return dict(value=r.value, ln_value=r.ln_value, quad_err=r.quad_err, tail_bound=r.tail_bound,
                cos_branch=bool(r.cos_branch), evals=r.evals)

def kappa(v1, v2, u):
    return _call(_lib.eckart_kappa, v1, v2, u).value

def kappa_si(v1_j_mol, v2_j_mol, nu_cm, t_k):
    return _call(_lib.eckart_kappa_si, v1_j_mol, v2_j_mol, nu_cm, t_k).value

def log_transmission(v1, v2, u, eps):
    return _lib.eckart_log_transmission(v1, v2, u, eps)

def transmission_si(e_j_mol, v1_j_mol, v2_j_mol, nu_cm):
    return _lib.eckart_transmission_si(e_j_mol, v1_j_mol, v2_j_mol, nu_cm)
