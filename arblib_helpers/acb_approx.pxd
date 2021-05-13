#Some functions that are required to do approximate arithmetic with acb (matrices).
#All wrapped functions are custom implementations which distinguishes it from "acb_mat_approx"

from sage.libs.arb.arb cimport *
from sage.libs.arb.types cimport *
from sage.libs.arb.acb cimport *

cdef extern from "arb_wrap.h":
    #Functions provided by arb that are missing in sage.libs.arb.acb_mat
    void acb_mat_approx_solve_triu(acb_mat_t X, const acb_mat_t U, const acb_mat_t B, int unit, long prec)
    void acb_mat_approx_solve_tril(acb_mat_t X, const acb_mat_t L, const acb_mat_t B, int unit, long prec)
    int acb_mat_approx_lu(long *P, acb_mat_t LU, const acb_mat_t A, long prec)
    void acb_mat_approx_solve_lu_precomp(acb_mat_t X, const long *perm, const acb_mat_t A, const acb_mat_t B, long prec)
    int acb_mat_approx_solve(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, long prec)
    int acb_mat_approx_inv(acb_mat_t X, const acb_mat_t A, long prec)
    void acb_mat_approx_mul(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)
    void acb_mat_approx_mul_classical(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)
    void acb_mat_neg(acb_mat_t dest, const acb_mat_t src)

    void acb_mat_window_init(acb_mat_t window, const acb_mat_t mat, long r1, long c1, long r2, long c2)
    void acb_mat_window_clear(acb_mat_t window)

    void arf_printd(const arf_t x, long d)
    void arf_set_d(arf_t res, double x)
    void acb_set_d_d(acb_t z, double x, double y)

cdef extern from "acb_approx_helpers.c":
    # Helper functions that we have implemented ourselves
    void acb_approx_mul(acb_t res, const acb_t x, const acb_t y, int prec)
    void acb_approx_add(acb_t res, const acb_t x, const acb_t y, int prec)
    void acb_approx_sub(acb_t res, const acb_t x, const acb_t y, int prec)
    void acb_approx_set(acb_t res, const acb_t x)
    void acb_approx_set_arb(acb_t res, const arb_t x)
    void acb_approx_div_arb(acb_t res, const acb_t x, const arb_t y, int prec)
    void acb_approx_mul_arb(acb_t res, const acb_t x, const arb_t y, int prec)
    void acb_approx_abs(arb_t r, const acb_t z, int prec)
    void arb_approx_hypot(arb_t z, const arb_t x, const arb_t y, int prec)
    void acb_approx_inv(acb_t z, const acb_t x, int prec)
    double complex acb_to_dc(acb_t z)

cdef extern from "acb_mat_approx_helpers.c":
    void acb_mat_approx_add(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, int prec)
    void acb_mat_approx_sub(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, int prec)
    void acb_mat_approx_norm(arb_t res, acb_mat_t x, int prec)
    void acb_mat_approx_scalar_mul(acb_mat_t res, const acb_mat_t A, const acb_t c, int prec)
    void acb_mat_approx_scalar_mul_arb(acb_mat_t res, const acb_mat_t A, const arb_t c, int prec)
    void acb_mat_approx_scalar_addmul(acb_mat_t res, acb_mat_t B, const acb_mat_t A, const acb_t c, int prec)
    void acb_mat_approx_dotc(acb_t res, acb_mat_t x, acb_mat_t y, int prec)
    void acb_mat_approx_dot(acb_t res, acb_mat_t x, acb_mat_t y, int prec)
    void acb_mat_change_prec(acb_mat_t res, acb_mat_t A, int prec)
    void acb_mat_approx_left_mul_diag(acb_mat_t res, const acb_mat_t D, const acb_mat_t A, int prec)
    void acb_mat_approx_right_mul_diag(acb_mat_t res, const acb_mat_t A, const acb_mat_t D, int prec)
    void acb_approx_complex_sign(acb_t res, acb_t z, arb_t z_abs, int prec)
    void lartg(acb_t c, acb_t s, acb_t r, acb_t f, acb_t g, int prec)

cdef extern from "acb_dft_helpers.c":
    void acb_compute_dft_matrix(acb_mat_t A, int N, int bit_prec)
    void acb_test_fft(int N, int bit_prec)

cdef extern from "acb_poly_approx_helpers.c":
    void acb_poly_approx_evaluate_horner(acb_t res, const acb_poly_t f, const acb_t a, int prec)