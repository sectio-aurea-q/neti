# eckart-ref

Accurate Eckart tunnelling transmission coefficients κ(T), checked against certified
ball-arithmetic enclosures.

The Eckart barrier is the standard one-dimensional tunnelling correction in transition-state
theory and is implemented in most kinetics codes. Its usual arrangement,
`P = 1 − (cosh(A−B) + cosh D)/(cosh(A+B) + cosh D)`, has two numerical defects:

1. **Cancellation.** `1 − (…)` returns exactly 0 once the true value falls below ~1e-16.
   At low temperature these deep-tunnelling energies can carry most of the thermal
   integral, so κ comes out too small.
2. **Small barriers.** For `4·sqrt(V1·V2) < hν`, 2πd is imaginary and `cosh(2πd)` must become
   `cos(2π|d|)` (Johnston & Heicklen 1962; Brown 1981). Using `cosh` with `|…|` makes κ too
   small, in some cases below 1.

eckart-ref evaluates `P = 2 sinh A sinh B / (cosh(A+B) + D)` in log space, with the correct
branch for D. It integrates the thermal average adaptively and reports a quadrature error
estimate together with a rigorous bound on the truncated tail.

## Accuracy

The test set is 1197 points of the certified grid in `conformance/eckart_certified.csv`. It
covers α₁, α₂ from 0.5 to 100, u* = hν/k_BT from 0.05 to 30, strong asymmetry and the cos
branch. Each value is an Arb enclosure whose relative radius is at most 1e-10, and in most
cases much smaller.

| | eckart-ref (float64) |
|---|---|
| points passing at 1e-12 | **1197 / 1197** |
| worst relative error | 4.0e-14 (α = 100, u* = 25, κ ≈ 1.3e98) |
| median relative error | 7e-16 |
| cos branch (102 points), worst | 1.0e-15 |

The largest errors appear where κ is very large. There they are limited by float64 rounding of
the input in exp(v_min), not by the quadrature. The textbook arrangement above, without
overflow guards but integrated with the same accurate quadrature, fails 256 of the 1197 points
at 1e-6 (`examples/naive_form.py`).

Reproduce:

```bash
cargo build --release
python3 conformance/conformance.py          # tests eckart-ref
python3 examples/naive_form.py              # tests the textbook arrangement
```

## Use

Rust:

```rust
let k = eckart_ref::kappa(v1, v2, u)?;          // barriers in k_B T, u = h nu / k_B T
let k = eckart_ref::kappa_si(v1_j_mol, v2_j_mol, nu_cm, t_k)?;   // CODATA 2018
println!("{} (quad err {:.1e})", k.value, k.quad_err);
```

Python (ctypes, no extra build tools):

```python
import eckart_ref as er          # python/eckart_ref.py, after cargo build --release
er.kappa(4.0, 4.0, 2.0)
er.kappa_si(160e3, 160e3, 1200.0, 200.0)
```

C / C++ / Fortran: `eckart_kappa`, `eckart_kappa_si`, `eckart_transmission` in
`libeckart_ref` (see `src/lib.rs`).

Command line:

```bash
eckart-ref reduced 4 4 2
eckart-ref si 160 160 1200 200          # kJ/mol, kJ/mol, cm^-1, K
eckart-ref grid points.csv              # rows v1,v2,u
```

## Testing your own code

```python
from conformance import run
run(my_kappa, tol=1e-6)        # my_kappa(v1, v2, u) with barriers in k_B T, u = h nu / k_B T
```

Definition used, the same as in Johnston & Heicklen (1962) and Brown (1981): α_i = 2πV_i/(hν).
Energies are measured from the higher asymptote, and κ = e^{v_min} ∫₀^∞ P(ε) e^{−ε} dε with
v_min = min(V1, V2)/k_BT. κ is symmetric in V1 and V2. Codes that use a coarse fixed energy grid
add a separate quadrature error, which this test also reports.

## Limits

- This is the one-dimensional Eckart model. It is not multidimensional tunnelling,
  instanton theory or SCT.
- eckart-ref is validated against certified values, but it is not a certified computation
  in every call. A ball-arithmetic mode for per-call enclosures is planned.
- The reference grid was computed with python-flint (Arb). The scripts that produce it are
  in `reference/`.

## Author

Sabrina Hackfort (sectio-aurea-q, MEGALODON Research), ORCID
[0009-0007-7458-2487](https://orcid.org/0009-0007-7458-2487).

Tools: Rust, Arb / python-flint, mpmath, Claude.

## License

AGPL-3.0-or-later.
