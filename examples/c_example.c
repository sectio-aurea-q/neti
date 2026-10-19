#include <stdio.h>
#include <math.h>
#include "eckart_ref.h"
int main(void) {
    EckartResult r;
    int rc = eckart_kappa_si(160e3, 160e3, -1200.0, 200.0, &r);   /* RMG-style negative frequency is fine */
    printf("rc=%d kappa=%.10g ln=%.10g cos=%d evals=%llu\n", rc, r.value, r.ln_value, r.cos_branch, (unsigned long long)r.evals);
    rc = eckart_kappa(0.0, 1.0, 1.0, &r);
    printf("invalid input: rc=%d value=%g\n", rc, r.value);
    printf("P(E=100 kJ/mol) = %.6e\n", eckart_transmission_si(100e3, 160e3, 160e3, 1200.0));
    return 0;
}
