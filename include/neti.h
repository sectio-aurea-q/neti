/* NETI C interface. Link with -lneti (target/release/libneti.{so,dylib}).
 * Return codes: 0 ok, 1 invalid input (non-finite or non-positive), 2 no convergence/overflow.
 * On error all fields of *out are NaN. out may be NULL. */
#ifndef NETI_H
#define NETI_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double value;      /* kappa (may be +inf if it exceeds double range; see ln_value) */
    double ln_value;   /* ln kappa */
    double quad_err;   /* quadrature error estimate (absolute); an estimate, not a bound */
    double tail_bound; /* rigorous bound on the truncated tail (absolute) */
    int32_t cos_branch;/* 1 if 4*sqrt(V1*V2) < h*nu (D = cos(2*pi*|d|)) */
    uint64_t evals;    /* integrand evaluations */
} NetiResult;

/* Reduced units: v1, v2 = barrier heights / k_B T (reactant and product side), u = h*nu / k_B T. */
int32_t neti_kappa(double v1, double v2, double u, NetiResult *out);

/* SI units: barriers in J/mol, |imaginary frequency| in cm^-1 (sign ignored), T in K.
 * Exact SI constants. Barriers are used as given (pass ZPE-corrected values if wanted). */
int32_t neti_kappa_si(double v1_j_mol, double v2_j_mol, double nu_cm, double t_k, NetiResult *out);

/* ln P(eps), eps = energy above the higher asymptote in units of k_B T. NaN on invalid input. */
double neti_log_transmission(double v1, double v2, double u, double eps);

/* P(E), E in J/mol above the higher asymptote, barriers in J/mol, nu in cm^-1. NaN on invalid input. */
double neti_transmission_si(double e_j_mol, double v1_j_mol, double v2_j_mol, double nu_cm);

#ifdef __cplusplus
}
#endif
#endif
