//! # NETI
//!
//! Eckart tunnelling transmission coefficient κ(T) computed without the two numerical
//! defects found in many implementations:
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
//! integral is truncated at `ε_max = v_min + 80`; since `P ≤ 1`, the omitted tail of κ is
//! rigorously bounded by `exp(v_min − ε_max) = exp(−80)` (absolute) and is reported.
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
    /// ln κ (finite even when κ overflows f64).
    pub ln_value: f64,
    /// QUADPACK-style quadrature error estimate (absolute). An estimate, not a bound;
    /// it includes a floor of 50·ε_mach·|κ| for rounding in the integrand.
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
            EckartError::NoConvergence => write!(f, "quadrature did not converge within the evaluation limit"),
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

fn check(x: f64, name: &'static str) -> Result<(), EckartError> {
    if !x.is_finite() { return Err(EckartError::NotFinite(name)); }
    if !(x > 0.0) { return Err(EckartError::NonPositive(name)); }
    Ok(())
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
        check(v1, "v1")?; check(v2, "v2")?; check(u, "u")?;
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

    /// ln of the two denominators' common form: ln(cosh(x) + D) for x ≥ 0.
    #[inline]
    fn log_cosh_plus_d(&self, x: f64) -> f64 {
        if self.q >= 0.0 {
            log_add_exp(log_cosh(x), self.log_d)
        } else {
            // cosh(x) + cos(2π|d|) = 2 sinh²(x/2) + 2 cos²(π|d|): two non-negative terms,
            // no cancellation even when cos(2π|d|) is close to −1.
            std::f64::consts::LN_2 + log_add_exp(2.0 * log_sinh_or_neginf(0.5 * x), self.log_cos2)
        }
    }

    /// (ln P, ln(1 − P)) at ε (energy above the higher asymptote, units of k_B T).
    ///
    /// P = 2 sinh A sinh B / (cosh(A+B) + D) and 1 − P = (cosh(A−B) + D)/(cosh(A+B) + D).
    /// Whichever of P, 1 − P is smaller is computed directly; the other via log1p/expm1,
    /// so both are accurate to a few ulp over the whole range.
    pub fn log_p_and_q(&self, eps: f64) -> (f64, f64) {
        if eps.is_nan() { return (f64::NAN, f64::NAN); }
        if eps <= 0.0 { return (f64::NEG_INFINITY, 0.0); }
        let a = 2.0 * (self.k * eps).sqrt() * self.inv_s;
        let b = 2.0 * (self.d_alpha + self.k * eps).sqrt() * self.inv_s;
        let log_den = self.log_cosh_plus_d(a + b);
        let log_p = std::f64::consts::LN_2 + log_sinh(a) + log_sinh(b) - log_den;
        if log_p < -std::f64::consts::LN_2 {
            (log_p, (-log_p.exp()).ln_1p())
        } else {
            let log_q = self.log_cosh_plus_d((a - b).abs()) - log_den;
            let log_q = if log_q > 0.0 { 0.0 } else { log_q };
            (( -log_q.exp_m1()).ln(), log_q)
        }
    }

    /// ln P(ε). NaN for NaN input, −∞ for ε ≤ 0.
    pub fn log_transmission(&self, eps: f64) -> f64 { self.log_p_and_q(eps).0 }

    /// P(ε), ε above the higher asymptote in units of k_B T.
    pub fn transmission(&self, eps: f64) -> f64 { self.log_transmission(eps).exp() }

    /// ln of the integrand in t = sqrt(ε): ln(2t) + ln P(t²) + v_min − t²
    #[inline]
    fn log_integrand_t(&self, t: f64) -> f64 {
        if t <= 0.0 { return f64::NEG_INFINITY; }
        let e = t * t;
        (2.0 * t).ln() + self.log_transmission(e) + self.v_min - e
    }
}

#[inline]
fn log_sinh_or_neginf(x: f64) -> f64 { if x > 0.0 { log_sinh(x) } else { f64::NEG_INFINITY } }

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

/// One G7/K15 step on [a, b] with the QUADPACK (qk15) error estimate.
/// Returns (result, abs error estimate, resabs).
pub fn gk15<F: Fn(f64) -> f64>(f: &F, a: f64, b: f64) -> (f64, f64, f64) {
    let c = 0.5 * (a + b);
    let h = 0.5 * (b - a);
    let mut fv1 = [0.0; 7];
    let mut fv2 = [0.0; 7];
    let fc = f(c);
    let mut resk = fc * WGK[7];
    let mut resg = fc * WG[3];
    let mut resabs = fc.abs() * WGK[7];
    for j in 0..7 {
        let x = h * XGK[j];
        let (f1, f2) = (f(c - x), f(c + x));
        fv1[j] = f1; fv2[j] = f2;
        resk += WGK[j] * (f1 + f2);
        resabs += WGK[j] * (f1.abs() + f2.abs());
        if j % 2 == 1 { resg += WG[j / 2] * (f1 + f2); }
    }
    let reskh = resk * 0.5;
    let mut resasc = WGK[7] * (fc - reskh).abs();
    for j in 0..7 { resasc += WGK[j] * ((fv1[j] - reskh).abs() + (fv2[j] - reskh).abs()); }
    let (result, resabs, resasc) = (resk * h, resabs * h.abs(), resasc * h.abs());
    let mut err = ((resk - resg) * h).abs();
    if resasc != 0.0 && err != 0.0 { err = resasc * (1.0f64).min((200.0 * err / resasc).powf(1.5)); }
    let floor = 50.0 * f64::EPSILON * resabs;
    if resabs > f64::MIN_POSITIVE / (50.0 * f64::EPSILON) { err = err.max(floor); }
    (result, err, resabs)
}

/// Relative tolerance requested from the quadrature (a floor of 50 ε_mach applies).
pub const REL_TOL: f64 = 1e-13;
/// Truncation distance above v_min (in k_B T).
pub const TAIL: f64 = 80.0;
/// Maximum integrand evaluations per κ.
pub const MAX_EVALS: usize = 200_000;

struct Panel { a: f64, b: f64, val: f64, err: f64 }

/// Global adaptive Gauss–Kronrod (QAG-like): always bisect the panel with the largest error.
fn qag<F: Fn(f64) -> f64>(f: &F, breaks: &[f64], rel_tol: f64) -> Result<(f64, f64, usize), EckartError> {
    let mut panels: Vec<Panel> = Vec::new();
    let mut evals = 0usize;
    for w in breaks.windows(2) {
        let (v, e, _) = gk15(f, w[0], w[1]);
        evals += 15;
        panels.push(Panel { a: w[0], b: w[1], val: v, err: e });
    }
    loop {
        let total: f64 = panels.iter().map(|p| p.val).sum();
        let err: f64 = panels.iter().map(|p| p.err).sum();
        if !total.is_finite() || !err.is_finite() { return Err(EckartError::NoConvergence); }
        let floor = 64.0 * f64::EPSILON * total.abs();
        if err <= (rel_tol * total.abs()).max(floor) { return Ok((total, err, evals)); }
        if evals >= MAX_EVALS {
            // accept if the remaining error is at the rounding level of the integrand
            if err <= 1e3 * floor { return Ok((total, err, evals)); }
            return Err(EckartError::NoConvergence);
        }
        let (i, _) = panels.iter().enumerate().fold((0, -1.0), |acc, (i, p)| if p.err > acc.1 { (i, p.err) } else { acc });
        let p = panels.swap_remove(i);
        let m = 0.5 * (p.a + p.b);
        if !(m > p.a && m < p.b) {
            // cannot split further: keep it as is but stop refining it
            panels.push(Panel { err: 0.0, ..p });
            continue;
        }
        let (v1, e1, _) = gk15(f, p.a, m);
        let (v2, e2, _) = gk15(f, m, p.b);
        evals += 30;
        panels.push(Panel { a: p.a, b: m, val: v1, err: e1 });
        panels.push(Panel { a: m, b: p.b, val: v2, err: e2 });
    }
}

/// κ in reduced units: barrier heights v1, v2 (in k_B T) and u = hν/(k_B T).
pub fn kappa(v1: f64, v2: f64, u: f64) -> Result<Kappa, EckartError> {
    let bar = Barrier::new(v1, v2, u)?;
    let vm = bar.v_min;
    let eps_max = vm + TAIL;
    // breakpoints in ε (scaled with the barrier, plus fixed points in k_B T for small barriers)
    let mut cand: Vec<f64> = vec![0.25 * vm, 0.5 * vm, 0.75 * vm, 0.9 * vm, vm, 1.1 * vm,
                                  vm + 2.0, vm + 8.0, vm + 25.0, 0.5, 2.0, 5.0, 12.0, 30.0];
    // around the barrier top P rises from 0 to 1 over a width of about u/(2π) k_B T
    // (parabolic limit); for wide barriers at high T this is far narrower than the panels
    let w = u / (2.0 * PI);
    for n in [0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0] { cand.push(vm - n * w); cand.push(vm + n * w); }
    cand.retain(|&x| x > 0.0 && x < eps_max);
    cand.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut br: Vec<f64> = vec![0.0];
    for x in cand { if x > *br.last().unwrap() * (1.0 + 1e-9) { br.push(x); } }
    br.push(eps_max);
    let ts: Vec<f64> = br.iter().map(|e| e.sqrt()).collect();
    // scale: integrate exp(log f − c) with c ≈ max log f, so neither overflow nor underflow occurs
    let mut c = f64::NEG_INFINITY;
    for w in ts.windows(2) {
        for j in 1..8 { c = c.max(bar.log_integrand_t(w[0] + (w[1] - w[0]) * j as f64 / 8.0)); }
    }
    if !c.is_finite() { return Err(EckartError::NoConvergence); }
    let f = |t: f64| (bar.log_integrand_t(t) - c).exp();
    let (val, err, evals) = qag(&f, &ts, REL_TOL)?;
    let ln_value = c + val.ln();
    let scale = c.exp();
    Ok(Kappa {
        value: ln_value.exp(),
        ln_value,
        quad_err: if scale.is_finite() { err * scale } else { f64::INFINITY },
        tail_bound: (vm - eps_max).exp(),
        cos_branch: bar.cos_branch(),
        evals,
    })
}

/// κ in SI-style units: barrier heights in J/mol from the reactant (`v1_j_mol`) and product
/// (`v2_j_mol`) sides, imaginary frequency in cm⁻¹ (sign ignored), temperature in K.
/// Exact SI constants (CODATA 2018/2022). Pass zero-point-corrected barriers if you want the
/// usual ZPE-corrected Eckart model; the library does not add or remove ZPE.
pub fn kappa_si(v1_j_mol: f64, v2_j_mol: f64, nu_cm: f64, t_k: f64) -> Result<Kappa, EckartError> {
    check(t_k, "T")?;
    check(v1_j_mol, "V1")?;
    check(v2_j_mol, "V2")?;
    if !nu_cm.is_finite() { return Err(EckartError::NotFinite("nu")); }
    let nu = nu_cm.abs();
    if !(nu > 0.0) { return Err(EckartError::NonPositive("|nu|")); }
    let rt = si::R * t_k;
    let u = si::H * si::C * 100.0 * nu / (si::KB * t_k);
    kappa(v1_j_mol / rt, v2_j_mol / rt, u)
}

/// Microcanonical P(E) in SI units: E in J/mol measured from the higher asymptote, barriers in
/// J/mol, imaginary frequency in cm⁻¹ (sign ignored). Independent of temperature.
pub fn transmission_si(e_j_mol: f64, v1_j_mol: f64, v2_j_mol: f64, nu_cm: f64) -> Result<f64, EckartError> {
    check(v1_j_mol, "V1")?;
    check(v2_j_mol, "V2")?;
    let nu = nu_cm.abs();
    check(nu, "|nu|")?;
    // reduced units with an arbitrary reference energy (here: hν), P is scale-invariant
    let hnu = si::H * si::C * 100.0 * nu * si::NA; // J/mol
    let bar = Barrier::new(v1_j_mol / hnu, v2_j_mol / hnu, 1.0)?;
    Ok(bar.transmission(e_j_mol / hnu))
}

// ─── C ABI (for C, C++, Fortran, Python ctypes) ─────────────────────────────
// Return codes: 0 ok, 1 invalid input, 2 no convergence / overflow. On error *value = NaN.

#[repr(C)]
pub struct NetiResult { pub value: f64, pub ln_value: f64, pub quad_err: f64, pub tail_bound: f64, pub cos_branch: i32, pub evals: u64 }

fn to_c(r: Result<Kappa, EckartError>, out: *mut NetiResult) -> i32 {
    let (code, k) = match r {
        Ok(k) => (0, k),
        Err(e) => (if matches!(e, EckartError::NoConvergence) { 2 } else { 1 },
                   Kappa { value: f64::NAN, ln_value: f64::NAN, quad_err: f64::NAN, tail_bound: f64::NAN, cos_branch: false, evals: 0 }),
    };
    if !out.is_null() {
        unsafe { *out = NetiResult { value: k.value, ln_value: k.ln_value, quad_err: k.quad_err, tail_bound: k.tail_bound, cos_branch: k.cos_branch as i32, evals: k.evals as u64 }; }
    }
    code
}

/// Reduced-unit κ (barriers in k_B T, u = hν/k_B T).
#[no_mangle]
pub extern "C" fn neti_kappa(v1: f64, v2: f64, u: f64, out: *mut NetiResult) -> i32 { to_c(kappa(v1, v2, u), out) }

/// SI-unit κ: barriers in J/mol, imaginary frequency in cm⁻¹ (sign ignored), T in K.
#[no_mangle]
pub extern "C" fn neti_kappa_si(v1_j_mol: f64, v2_j_mol: f64, nu_cm: f64, t_k: f64, out: *mut NetiResult) -> i32 {
    to_c(kappa_si(v1_j_mol, v2_j_mol, nu_cm, t_k), out)
}

/// ln P(ε) in reduced units (ε above the higher asymptote in k_B T). NaN on invalid input.
#[no_mangle]
pub extern "C" fn neti_log_transmission(v1: f64, v2: f64, u: f64, eps: f64) -> f64 {
    match Barrier::new(v1, v2, u) { Ok(b) => b.log_transmission(eps), Err(_) => f64::NAN }
}

/// P(E) in SI units (E, V1, V2 in J/mol, E from the higher asymptote; ν in cm⁻¹). NaN on invalid input.
#[no_mangle]
pub extern "C" fn neti_transmission_si(e_j_mol: f64, v1_j_mol: f64, v2_j_mol: f64, nu_cm: f64) -> f64 {
    transmission_si(e_j_mol, v1_j_mol, v2_j_mol, nu_cm).unwrap_or(f64::NAN)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gk15_integrates_polynomials_exactly() {
        // Kronrod-15 is exact up to degree 22 (Gauss-7 up to 13).
        for deg in 0..=22u32 {
            let f = |x: f64| x.powi(deg as i32);
            let (k, _, _) = gk15(&f, 0.0, 1.0);
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
    /// Wide barriers at high temperature: P rises from 0 to 1 within ~u/(2π) of the barrier
    /// top. Values from an independent 50-digit reference (validation/independent_mpmath_ref.py).
    #[test]
    fn wide_barriers_high_temperature() {
        for (v1, v2, u, exact) in [(500.0, 500.0, 0.01, 1.000_004_175_012_131_9), (250.0, 250.0, 0.0123, 1.000_006_328_992_720_4)] {
            let k = kappa(v1, v2, u).unwrap();
            assert!((k.value / exact - 1.0).abs() < 1e-13, "({v1},{v2},{u}): {} vs {exact}", k.value);
            assert!(k.evals < MAX_EVALS);
        }
    }

    #[test]
    fn huge_kappa_reports_ln_instead_of_hanging() {
        let k = kappa(800.0, 800.0, 400.0).unwrap();
        assert!(k.value.is_infinite() && k.ln_value.is_finite() && k.ln_value > 700.0);
    }

    #[test]
    fn nan_energy_propagates() {
        let b = Barrier::new(4.0, 4.0, 2.0).unwrap();
        assert!(b.log_transmission(f64::NAN).is_nan());
        assert_eq!(b.log_transmission(0.0), f64::NEG_INFINITY);
    }

    #[test]
    fn p_close_to_one_is_accurate() {
        // far above the top, 1 − P is computed directly; P + (1 − P) must be 1
        let b = Barrier::new(20.0, 20.0, 0.01).unwrap();
        let (lp, lq) = b.log_p_and_q(25.0);
        assert!((lp.exp() + lq.exp() - 1.0).abs() < 1e-15);
    }
}
