#Some functions that are missing in sage.libs.arb.acb_mat

# distutils: depends = acb_mat.h

from sage.libs.arb.types cimport acb_t, acb_ptr, acb_srcptr, acb_mat_t, acb_poly_t, mag_t, arf_t

# acb_mat.h
cdef extern from "arb_wrap.h":
    void acb_mat_approx_solve_triu(acb_mat_t X, const acb_mat_t U, const acb_mat_t B, int unit, long prec)
    void acb_mat_approx_solve_tril(acb_mat_t X, const acb_mat_t L, const acb_mat_t B, int unit, long prec)
    int acb_mat_approx_lu(long *P, acb_mat_t LU, const acb_mat_t A, long prec)
    void acb_mat_approx_solve_lu_precomp(acb_mat_t X, const long *perm, const acb_mat_t A, const acb_mat_t B, long prec)
    int acb_mat_approx_solve(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, long prec)
    int acb_mat_approx_inv(acb_mat_t X, const acb_mat_t A, long prec)
    void acb_mat_approx_mul(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)
    void acb_mat_approx_mul_classical(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)

    void acb_mat_window_init(acb_mat_t window, const acb_mat_t mat, long r1, long c1, long r2, long c2)
    void acb_mat_window_clear(acb_mat_t window)

    void arf_printd(const arf_t x, long d)