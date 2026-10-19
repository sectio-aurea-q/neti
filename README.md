# eckart-ref

eckart-ref computes accurate Eckart tunnelling transmission coefficients κ(T) and checks
them against certified ball-arithmetic enclosures.

The Eckart barrier is the standard one-dimensional tunnelling correction in transition-state
theory, and most kinetics codes implement it. The usual arrangement,
`P = 1 − (cosh(A−B) + cosh D)/(cosh(A+B) + cosh D)`, has two numerical defects:

1. **Cancellation.** `1 − (…)` returns exactly 0 once the true value is below about 1e-16.
   At low temperature these deep-tunnelling energies can carry most of the thermal integral,
   so κ comes out too small.
2. **Small barriers.** For `4·sqrt(V1·V2) < hν`, 2πd is imaginary, and `cosh(2πd)` must become
   `cos(2π|d|)` (Johnston & Heicklen 1962; Brown 1981). Using `cosh` with `|…|` makes κ too
   small, sometimes below 1.

eckart-ref fixes both, and avoids two further losses of accuracy:

- **No cancellation in P.** It evaluates `P = 2 sinh A sinh B / (cosh(A+B) + D)` in log
  space, with the correct branch for D.
- **No cancellation in the cos branch.** It writes the denominator there as
  `2 sinh²((A+B)/2) + 2 cos²(π|d|)`, a sum of two non-negative terms, so it stays accurate
  even when cos(2π|d|) ≈ −1 (very small barriers).
- **Accurate 1 − P.** Where P > 1/2 it computes 1 − P directly, so the region above the
  barrier top is also accurate.
- **Adaptive integration.** It integrates the thermal average with global adaptive
  Gauss–Kronrod quadrature in t = √ε. The breakpoints follow the barrier top, where P rises
  over a width of about hν/(2π k_BT). Each call has an evaluation limit.
- **Reported uncertainty.** It returns ln κ, which stays finite when κ exceeds float64 range,
  a quadrature error estimate, and a rigorous bound on the truncated tail.

## Accuracy

| Test | Points | Worst relative error | Median |
|---|---|---|---|
| Certified Arb grid, α 0.5–100, u* 0.05–30 (`conformance/eckart_certified.csv`) | 1197 / 1197 within 1e-12 | 4.5e-14 | 4.4e-16 |
| — of which cos branch | 102 | 7.8e-16 | |
| Independent 50-digit reference, random, v up to 500, u 0.01–80 | 303 / 303 | 7.7e-14 | 2.2e-16 |
| Independent 50-digit reference, random, v 1e-4 – 5e4, u 0.0035–300 (`validation/sweep2.py`) | 197 / 197 (185 within 1e-13) | 2.2e-12 | < 1e-16 |

The independent reference (`validation/independent_mpmath_ref.py`) does not use the
α-parametrisation. It starts from the Eckart potential itself, with the closed-form P(E)
checked against a direct numerical solution of the Schrödinger equation.

The largest errors occur where κ is very large (κ ≈ 1e98 at α = 100, u* = 25). They come from
float64 rounding in the large exponent e^{v_min}: about half from rounding the input, about
half from the library's own arithmetic.

For comparison, the textbook arrangement above fails 256 of the 1197 certified points at
1e-6 (102 in the cos branch, 121 from cancellation, 33 from overflow). It was integrated with
the same accurate quadrature but without overflow guards (`examples/naive_form.py`).

A typical call needs about 500–800 integrand evaluations.

Reproduce:

```bash
cargo build --release && cargo test --release
python3 conformance/conformance.py          # eckart-ref against the certified grid
python3 validation/validate.py              # same, with statistics
python3 examples/naive_form.py              # the textbook arrangement
```

## Use

Rust:

```rust
let k = eckart_ref::kappa(v1, v2, u)?;          // barriers in k_B T, u = h nu / k_B T
let k = eckart_ref::kappa_si(v1_j_mol, v2_j_mol, nu_cm, t_k)?;   // J/mol, J/mol, cm^-1, K
println!("{} (ln {}, quad err {:.1e})", k.value, k.ln_value, k.quad_err);
let p = eckart_ref::transmission_si(e_j_mol, v1_j_mol, v2_j_mol, nu_cm)?;  // P(E)
```

Python (ctypes; build the library first):

```python
import eckart_ref as er          # python/eckart_ref.py
er.kappa(4.0, 4.0, 2.0)
er.kappa_si(160e3, 160e3, -1200.0, 200.0)      # negative (imaginary) frequency is accepted
er.kappa_full(800, 800, 400)                   # dict with value (inf here), ln_value, quad_err, ...
```

C / C++ / Fortran: see `include/eckart_ref.h` and `examples/c_example.c`
(`gcc -Iinclude examples/c_example.c -Ltarget/release -leckart_ref`).

Command line:

```bash
eckart-ref reduced 4 4 2
eckart-ref si 160 160 1200 200          # kJ/mol, kJ/mol, cm^-1, K  (J/mol in the API)
eckart-ref grid points.csv              # rows v1,v2,u; one output row per input row
```

## Conventions

- **Barriers.** V1 and V2 are the barrier heights seen from the reactant and the product
  side. κ is symmetric in them. They are used exactly as given: pass zero-point-corrected
  barriers if you want the usual ZPE-corrected Eckart model.
- **Energy zero.** Energies are measured from the higher of the two asymptotes, and
  κ = e^{v_min} ∫₀^∞ P(ε) e^{−ε} dε with v_min = min(V1, V2)/k_BT. Codes that measure E from
  the barrier top or from the reactant shift by a constant.
- **Frequency.** ν is the magnitude of the imaginary frequency; its sign is ignored.
- **Constants.** These are the exact SI values (identical in CODATA 2018 and 2022). κ depends
  on the constants through about v_min·δ, so compare with codes that use older constants in
  reduced units (v1, v2, u).
- **Definition.** α_i = 2πV_i/(hν), as in Johnston & Heicklen (1962) and Brown (1981).

## Testing your own code

```python
from conformance import run      # conformance/conformance.py
run(my_kappa, tol=1e-6)          # my_kappa(v1, v2, u), barriers in k_B T, u = h nu / k_B T
```

How the certified grid was made:
- It was computed with python-flint (Arb), using the scripts in `reference/`.
- `conformance/build_conformance.py` assembles it.
- Each row is an enclosure for the IEEE double nearest to the listed grid value, and the
  `ustar` column gives that double exactly.
- The bounds are widened by a relative 1e-24 to account for decimal printing.
- Codes that use a coarse fixed energy grid add their own quadrature error, which the test
  also shows.

## Limits

- This is the one-dimensional Eckart model. It is not multidimensional tunnelling,
  instanton theory or small-curvature tunnelling.
- eckart-ref is validated against certified and independent high-precision values. It is
  not a certified computation in every call. A ball-arithmetic mode for per-call enclosures
  is planned.
- The certified grid covers α ≤ 100. Beyond that the evidence is the independent 50-digit
  reference; the error grows slowly with α (2e-12 at α ≈ 1e6).

## Author

Sabrina Hackfort (sectio-aurea-q, MEGALODON Research), ORCID
[0009-0007-7458-2487](https://orcid.org/0009-0007-7458-2487).

Tools: Rust, Arb / python-flint, mpmath, Claude.

## License

LICENSE_LINE
