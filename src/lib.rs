//! # eckart-ref
//!
//! Eckart tunnelling transmission coefficient κ(T) computed without the two numerical
//! defects found in most implementations (see the audit in this repository):
//!
//! 1. **No cancellation.** The transmission probability is evaluated as
//!    `P = 2 sinh(A) sinh(B) / (cosh(A+B) + D)` in log space, instead of
//!    `1 - (cosh(A-B) + D)/(cosh(A+B) + D)`, which returns 0 below ~1e-16.
//! 2. **Small-barrier branch.** `D = cosh(2πd)` if `α₁α₂ ≥ π²/4`, otherwise
//!    `D = cos(2π|d|)` (Johnston & Heicklen 1962; Brown 1981). Many codes use `cosh`
//!    with `|…|` in both cases.
//!
//! The thermal average is integrated with adaptive Gauss–Kronrod (7/15) quadrature in the
//! variable `t = sqrt(ε)`, which removes the square-root behaviour at the threshold. The
//! integral is truncated at `ε_max = v_min + 80`; since `P ≤ 1`, the omitted tail is
//! rigorously bounded by `exp(v_min − ε_max) = exp(−80)` (in units of κ's integrand) and
//! is reported.
//!
//! ## Conventions (Johnston & Heicklen 1962, Brown 1981)
//! * `v1`, `v2`: barrier heights from the reactant and product sides, in units of k_B T (> 0).
//! * `u = h ν / (k_B T)` with ν the magnitude of the imaginary frequency (> 0).
//! * `α_i = 2π v_i / u`.
//! * Energies are measured from the higher of the two asymptotes, so
//!   `κ = exp(v_min) ∫₀^∞ P(ε) exp(−ε) dε`, with `v_min = min(v1, v2)`.
//!
//! κ is symmetric in (v1, v2).

use std::f64::consts::PI;

/// CODATA 2018 exact constants (SI).
pub mod si {
    pub const H: f64 = 6.626_070_15e-34;
    pub const C: f64 = 299_792_458.0;
    pub const KB: f64 = 1.380_649e-23;
    pub const NA: f64 = 6.022_140_76e23;
    pub const R: f64 = KB * NA;
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kappa {
    /// κ(T)
    pub value: f64,
    /// Quadrature error estimate (absolute), sum of |Kronrod − Gauss| over panels.
    /// Conservative: the Kronrod result is normally far more accurate than this.
    pub quad_err: f64,
    /// Rigorous upper bound on the truncated tail (absolute).
    pub tail_bound: f64,
    /// True if the small-barrier branch D = cos(2π|d|) applies.
    pub cos_branch: bool,
    /// Number of integrand evaluations.
    pub evals: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum EckartError {
    NonPositive(&'static str),
    NotFinite(&'static str),
    NoConvergence,
}

impl std::fmt::Display for EckartError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EckartError::NonPositive(s) => write!(f, "{s} must be > 0"),
            EckartError::NotFinite(s) => write!(f, "{s} must be finite"),
            EckartError::NoConvergence => write!(f, "quadrature did not converge"),
        }
    }
}
impl std::error::Error for EckartError {}

#[inline]
fn log_sinh(x: f64) -> f64 {
    // log(sinh x), x > 0, without overflow
    if x > 1.0 { x - std::f64::consts::LN_2 + (-(-2.0 * x).exp()).ln_1p() } else { x.sinh().ln() }
}

#[inline]
fn log_cosh(x: f64) -> f64 {
    // log(cosh x), x >= 0, without overflow
    if x > 1.0 { x - std::f64::consts::LN_2 + (-2.0 * x).exp().ln_1p() } else { x.cosh().ln() }
}

#[inline]
fn log_add_exp(x: f64, y: f64) -> f64 {
    if x > y { x + (y - x).exp().ln_1p() } else { y + (x - y).exp().ln_1p() }
}

/// Precomputed barrier geometry.
#[derive(Debug, Clone, Copy)]
pub struct Barrier {
    v_min: f64,
    k: f64,        // 2π/u
    d_alpha: f64,  // α_hi − α_lo, computed without cancellation
    inv_s: f64,    // 1 / (α_lo^{-1/2} + α_hi^{-1/2})
    q: f64,        // α_lo α_hi − π²/4
    log_d: f64,    // log D for the cosh branch (q ≥ 0)
    log_cos2: f64, // log cos²(π|d|) for the cos branch (q < 0)
}

impl Barrier {
    pub fn new(v1: f64, v2: f64, u: f64) -> Result<Self, EckartError> {
        for (x, n) in [(v1, "v1"), (v2, "v2"), (u, "u")] {
            if !x.is_finite() { return Err(EckartError::NotFinite(n)); }
            if !(x > 0.0) { return Err(EckartError::NonPositive(n)); }
        }
        let (lo, hi) = if v1 <= v2 { (v1, v2) } else { (v2, v1) };
        let k = 2.0 * PI / u;
        let a_lo = k * lo;
        let a_hi = k * hi;
        let q = a_lo * a_hi - PI * PI / 4.0;
        let (log_d, log_cos2) = if q >= 0.0 {
            (log_cosh(2.0 * q.sqrt()), f64::NAN)
        } else {
            // cos(π|d|) with π|d| = sqrt(π²/4 − α_lo α_hi); written as sin(π/2 − π|d|) and
            // π/2 − sqrt(π²/4 − p) = p / (π/2 + sqrt(π²/4 − p)), free of cancellation.
            let p = a_lo * a_hi;
            let w = p / (PI / 2.0 + (PI * PI / 4.0 - p).sqrt());
            (f64::NAN, 2.0 * w.sin().ln())
        };
        Ok(Barrier {
            v_min: lo,
            k,
            d_alpha: k * (hi - lo),
            inv_s: 1.0 / (1.0 / a_lo.sqrt() + 1.0 / a_hi.sqrt()),
            q,
            log_d,
            log_cos2,
        })
    }

    pub fn cos_branch(&self) -> bool { self.q < 0.0 }

    /// log P(ε), ε = energy above the higher asymptote in units of k_B T.
    pub fn log_transmission(&self, eps: f64) -> f64 {
        if !(eps > 0.0) { return f64::NEG_INFINITY; }
        let a = 2.0 * (self.k * eps).sqrt() * self.inv_s;
        let b = 2.0 * (self.d_alpha + self.k * eps).sqrt() * self.inv_s;
        let log_num = std::f64::consts::LN_2 + log_sinh(a) + log_sinh(b);
        let log_den = if self.q >= 0.0 {
            log_add_exp(log_cosh(a + b), self.log_d)
        } else {
            // cosh(S) + cos(2π|d|) = 2 sinh²(S/2) + 2 cos²(π|d|): two non-negative terms,
            // no cancellation even when cos(2π|d|) is close to −1.
            std::f64::consts::LN_2 + log_add_exp(2.0 * log_sinh(0.5 * (a + b)), self.log_cos2)
        };
        (log_num - log_den).min(0.0)
    }

    /// P(ε)
    pub fn transmission(&self, eps: f64) -> f64 { self.log_transmission(eps).exp() }

    /// Integrand in t = sqrt(ε): 2 t P(t²) exp(v_min − t²)
    #[inline]
    fn integrand_t(&self, t: f64) -> f64 {
        if t <= 0.0 { return 0.0; }
        let e = t * t;
        2.0 * t * (self.log_transmission(e) + self.v_min - e).exp()
    }
}

// Gauss–Kronrod 7/15 nodes and weights (QUADPACK qk15).
const XGK: [f64; 8] = [
    0.991_455_371_120_812_639_206_854_697_526_329,
    0.949_107_912_342_758_524_526_189_684_047_851,
    0.864_864_423_359_769_072_789_712_788_640_926,
    0.741_531_185_599_394_439_863_864_773_280_788,
    0.586_087_235_467_691_130_294_144_845_693_013,
    0.405_845_151_377_397_166_906_606_412_076_961,
    0.207_784_955_007_898_467_600_689_403_773_245,
    0.0,
];
const WGK: [f64; 8] = [
    0.022_935_322_010_529_224_963_732_008_058_970,
    0.063_092_092_629_978_553_290_700_663_189_204,
    0.104_790_010_322_250_183_839_876_322_541_518,
    0.140_653_259_715_525_918_745_189_590_510_238,
    0.169_004_726_639_267_902_826_583_426_598_550,
    0.190_350_578_064_785_409_913_256_402_421_014,
    0.204_432_940_075_298_892_414_161_999_234_649,
    0.209_482_141_084_727_828_012_999_174_891_714,
];
const WG: [f64; 4] = [
    0.129_484_966_168_869_693_270_611_432_679_082,
    0.279_705_391_489_276_667_901_467_771_423_780,
    0.381_830_050_505_118_944_950_369_775_488_975,
    0.417_959_183_673_469_387_755_102_040_816_327,
];

/// One G7/K15 step on [a, b]: (Kronrod estimate, |Kronrod − Gauss|)
pub fn gk15<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64) -> (f64, f64) {
    let c = 0.5 * (a + b);
    let h = 0.5 * (b - a);
    let fc = f(c);
    let mut k = fc * WGK[7];
    let mut g = fc * WG[3];
    for j in 0..7 {
        let x = h * XGK[j];
        let s = f(c - x) + f(c + x);
        k += WGK[j] * s;
        if j % 2 == 1 { g += WG[j / 2] * s; }
    }
    (k * h, ((k - g) * h).abs())
}

fn adapt<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64, tol: f64, depth: u32, evals: &mut usize) -> Result<(f64, f64), EckartError> {
    let (k, e) = gk15(f, a, b);
    *evals += 15;
    if e <= tol || depth == 0 || (b - a) <= f64::EPSILON * a.abs().max(1.0) * 16.0 {
        if depth == 0 && e > tol * 1e3 { return Err(EckartError::NoConvergence); }
        return Ok((k, e));
    }
    let m = 0.5 * (a + b);
    let (k1, e1) = adapt(f, a, m, 0.5 * tol, depth - 1, evals)?;
    let (k2, e2) = adapt(f, m, b, 0.5 * tol, depth - 1, evals)?;
    Ok((k1 + k2, e1 + e2))
}

/// Relative tolerance used by [`kappa`].
pub const REL_TOL: f64 = 1e-13;
/// Truncation distance above v_min (in k_B T).
pub const TAIL: f64 = 80.0;

/// κ in reduced units: barrier heights v1, v2 (in k_B T) and u = hν/(k_B T).
pub fn kappa(v1: f64, v2: f64, u: f64) -> Result<Kappa, EckartError> {
    let bar = Barrier::new(v1, v2, u)?;
    let f = |t: f64| bar.integrand_t(t);
    let vm = bar.v_min;
    let eps_max = vm + TAIL;
    // breakpoints in ε: threshold region, around the barrier top, tail
    // (scaled with the barrier, plus fixed points in units of k_B T for small barriers)
    let mut cand: Vec<f64> = vec![0.25 * vm, 0.5 * vm, 0.75 * vm, 0.9 * vm, vm, 1.1 * vm,
                                  vm + 2.0, vm + 8.0, vm + 25.0, 0.5, 2.0, 5.0, 12.0, 30.0];
    cand.retain(|&x| x > 0.0 && x < eps_max);
    cand.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut br: Vec<f64> = vec![0.0];
    for x in cand {
        if x > *br.last().unwrap() * (1.0 + 1e-9) { br.push(x); }
    }
    br.push(eps_max);
    let ts: Vec<f64> = br.iter().map(|e| e.sqrt()).collect();
    // first pass to get the scale for the tolerance
    let mut evals = 0usize;
    let mut rough = 0.0;
    for w in ts.windows(2) { rough += gk15(&f, w[0], w[1]).0; evals += 15; }
    let tol_total = REL_TOL * rough.abs().max(f64::MIN_POSITIVE);
    let mut val = 0.0;
    let mut err = 0.0;
    let n = (ts.len() - 1) as f64;
    for w in ts.windows(2) {
        let (k, e) = adapt(&f, w[0], w[1], tol_total / n, 40, &mut evals)?;
        val += k;
        err += e;
    }
    Ok(Kappa { value: val, quad_err: err, tail_bound: (vm - eps_max).exp(), cos_branch: bar.cos_branch(), evals })
}

/// κ in SI-style units: barrier heights in J/mol from the reactant (`v1_j_mol`) and product
/// (`v2_j_mol`) sides, |imaginary frequency| in cm⁻¹, temperature in K. CODATA 2018 constants.
pub fn kappa_si(v1_j_mol: f64, v2_j_mol: f64, nu_cm: f64, t_k: f64) -> Result<Kappa, EckartError> {
    if !(t_k > 0.0) { return Err(EckartError::NonPositive("T")); }
    let rt = si::R * t_k;
    let u = si::H * si::C * 100.0 * nu_cm / (si::KB * t_k);
    kappa(v1_j_mol / rt, v2_j_mol / rt, u)
}

// ─── C ABI (for C, C++, Fortran, Python ctypes) ─────────────────────────────

/// Reduced-unit κ. Returns 0 on success, 1 on invalid input, 2 on non-convergence.
/// `value` and `quad_err` may be null.
#[no_mangle]
pub extern "C" fn eckart_kappa(v1: f64, v2: f64, u: f64, value: *mut f64, quad_err: *mut f64) -> i32 {
    match kappa(v1, v2, u) {
        Ok(k) => {
            unsafe {
                if !value.is_null() { *value = k.value; }
                if !quad_err.is_null() { *quad_err = k.quad_err; }
            }
            0
        }
        Err(EckartError::NoConvergence) => 2,
        Err(_) => 1,
    }
}

/// SI-unit κ: barriers in J/mol, |imaginary frequency| in cm⁻¹, T in K (CODATA 2018).
#[no_mangle]
pub extern "C" fn eckart_kappa_si(v1_j_mol: f64, v2_j_mol: f64, nu_cm: f64, t_k: f64, value: *mut f64, quad_err: *mut f64) -> i32 {
    match kappa_si(v1_j_mol, v2_j_mol, nu_cm, t_k) {
        Ok(k) => {
            unsafe {
                if !value.is_null() { *value = k.value; }
                if !quad_err.is_null() { *quad_err = k.quad_err; }
            }
            0
        }
        Err(EckartError::NoConvergence) => 2,
        Err(_) => 1,
    }
}

/// Transmission probability P(ε), ε above the higher asymptote in units of k_B T.
/// Returns NaN on invalid barrier parameters.
#[no_mangle]
pub extern "C" fn eckart_transmission(v1: f64, v2: f64, u: f64, eps: f64) -> f64 {
    match Barrier::new(v1, v2, u) { Ok(b) => b.transmission(eps), Err(_) => f64::NAN }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gk15_integrates_polynomials_exactly() {
        // Kronrod-15 is exact up to degree 22 (Gauss-7 up to 13).
        for deg in 0..=22u32 {
            let f = |x: f64| x.powi(deg as i32);
            let (k, _) = gk15(&f, 0.0, 1.0);
            let exact = 1.0 / (deg as f64 + 1.0);
            assert!((k - exact).abs() < 1e-14, "degree {deg}: {k} vs {exact}");
        }
    }

    #[test]
    fn symmetric_in_barriers() {
        let a = kappa(5.0, 9.0, 4.0).unwrap().value;
        let b = kappa(9.0, 5.0, 4.0).unwrap().value;
        assert!((a / b - 1.0).abs() < 1e-14);
    }

    #[test]
    fn deep_tunnelling_probability_not_zero() {
        // A=5, B=5, 2πd=80 -> P = 3.975088773825091e-31 (60-digit mpmath)
        let lp = std::f64::consts::LN_2 + log_sinh(5.0) + log_sinh(5.0) - log_add_exp(log_cosh(10.0), log_cosh(80.0));
        assert!((lp.exp() / 3.975088773825091e-31 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_input() {
        assert!(kappa(0.0, 1.0, 1.0).is_err());
        assert!(kappa(1.0, f64::NAN, 1.0).is_err());
        assert!(kappa(1.0, 1.0, -1.0).is_err());
    }
    /// Values certified with Arb (ball arithmetic, 333 bits) by reference/kappa_sq.py.
    #[test]
    fn matches_certified_edge_cases() {
        let cases: [(f64, f64, f64, f64, f64); 6] = [
            (0.999_999_9, 0.999_999_9, 4.0, 1.813_079_015_081_841_1, 1e-14), // q just below 0 (cos)
            (1.000_000_1, 1.000_000_1, 4.0, 1.813_079_101_074_999_1, 1e-14), // q just above 0
            (0.05, 200.0, 8.0, 0.365_082_203_731_967_56, 1e-14),               // strong asymmetry, κ < 1
            (40.0, 45.0, 60.0, 4.406_813_543_707_574_8e14, 1e-13),             // deep tunnelling
            (0.5, 0.7, 0.01, 1.000_011_109_665_044_7, 1e-13),                  // high temperature
            (1.0, 1.0, 3.0, 1.580_149_972_097_182_7, 1e-14),
        ];
        for (v1, v2, u, exact, tol) in cases {
            let k = kappa(v1, v2, u).unwrap().value;
            assert!((k / exact - 1.0).abs() < tol, "({v1},{v2},{u}): {k} vs {exact}");
        }
    }

    #[test]
    fn tiny_barriers_converge_quickly() {
        // cos branch with cos(2π|d|) ≈ −1: needs the cancellation-free denominator
        let k = kappa(1e-3, 1e-3, 1.0).unwrap();
        assert!((k.value / 1.001_000_425_104_319_4 - 1.0).abs() < 1e-14); // 50-digit mpmath
        assert!(k.evals < 5_000, "evals {}", k.evals);
        assert!(kappa(1e-8, 1e-8, 0.5).is_ok());
    }
}
