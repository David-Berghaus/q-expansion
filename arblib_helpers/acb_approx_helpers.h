#include "acb.h"
#include "acb_mat.h"

#include "acb.h"
// #include <complex.h>

void
acb_approx_mul(acb_t res, const acb_t x, const acb_t y, int prec);

void
acb_approx_add(acb_t res, const acb_t x, const acb_t y, int prec);

void
acb_approx_sub(acb_t res, const acb_t x, const acb_t y, int prec);

void
acb_approx_set(acb_t res, const acb_t x);

void
acb_approx_set_arb(acb_t res, const arb_t x);

void
acb_approx_swap(acb_t res, acb_t x);

void
acb_approx_div_arb(acb_t res, const acb_t x, const arb_t y, int prec);

void
acb_approx_inv(acb_t z, const acb_t x, int prec);

void
acb_approx_div(acb_t z, const acb_t x, const acb_t y, int prec);

void
acb_approx_mul_arb(acb_t res, const acb_t x, const arb_t y, int prec);

void
arb_approx_hypot(arb_t z, const arb_t x, const arb_t y, int prec);

void acb_approx_abs(arb_t r, const acb_t z, int prec);

// double complex acb_to_dc(acb_t z);