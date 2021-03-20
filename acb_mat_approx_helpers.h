//We declare this as a header file because we have *.c in gitignore...

#include "acb.h"
#include "acb_mat.h"

static void
acb_approx_mag(mag_t res, const acb_t x)
{
    mag_t t;
    mag_init(t);
    arf_get_mag(res, arb_midref(acb_realref(x)));
    arf_get_mag(t, arb_midref(acb_imagref(x)));
    mag_hypot(res, res, t);
    mag_clear(t);
}

static void
acb_approx_mul(acb_t res, const acb_t x, const acb_t y, long prec)
{
    arf_complex_mul(arb_midref(acb_realref(res)), arb_midref(acb_imagref(res)),
        arb_midref(acb_realref(x)), arb_midref(acb_imagref(x)), 
        arb_midref(acb_realref(y)), arb_midref(acb_imagref(y)), prec, ARF_RND_DOWN);
}

static void
acb_approx_add(acb_t res, const acb_t x, const acb_t y, long prec)
{
    arf_add(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)), arb_midref(acb_realref(y)), prec, ARF_RND_DOWN);
    arf_add(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)), arb_midref(acb_imagref(y)), prec, ARF_RND_DOWN);
}

static void
acb_approx_sub(acb_t res, const acb_t x, const acb_t y, long prec)
{
    arf_sub(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)), arb_midref(acb_realref(y)), prec, ARF_RND_DOWN);
    arf_sub(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)), arb_midref(acb_imagref(y)), prec, ARF_RND_DOWN);
}

static void
acb_approx_set(acb_t res, const acb_t x)
{
    arf_set(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)));
    arf_set(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)));
}

static void
acb_approx_div_arb(acb_t res, const acb_t x, const arb_t y, long prec)
{
    arf_div(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)), arb_midref(y), prec, ARF_RND_DOWN);
    arf_div(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)), arb_midref(y), prec, ARF_RND_DOWN);
}

static void
acb_approx_inv(acb_t z, const acb_t x, long prec)
{
    arf_set(arb_midref(acb_realref(z)), arb_midref(acb_realref(x)));
    arf_set(arb_midref(acb_imagref(z)), arb_midref(acb_imagref(x)));

    mag_zero(arb_radref(acb_realref(z)));
    mag_zero(arb_radref(acb_imagref(z)));

    acb_inv(z, z, prec);

    mag_zero(arb_radref(acb_realref(z)));
    mag_zero(arb_radref(acb_imagref(z)));
}

static void
acb_approx_div(acb_t z, const acb_t x, const acb_t y, long prec)
{
    acb_t t;
    acb_init(t);
    acb_approx_inv(t, y, prec);
    acb_approx_mul(z, x, t, prec);
    acb_clear(t);
}

static void
acb_approx_mul_arb(acb_t res, const acb_t x, const arb_t y, long prec)
{
    arf_mul(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)), arb_midref(y), prec, ARF_RND_DOWN);
    arf_mul(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)), arb_midref(y), prec, ARF_RND_DOWN);
}

void
acb_mat_approx_add(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)
{
    long i, j, r, c;
    r = acb_mat_nrows(mat1);
    c = acb_mat_ncols(mat1);

    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            acb_approx_add(acb_mat_entry(res, i, j), acb_mat_entry(mat1, i, j), acb_mat_entry(mat2, i, j), prec);
}

void //Computes A*c+B and stores it in res
acb_mat_approx_scalar_addmul(acb_mat_t res, acb_mat_t B, const acb_mat_t A, const acb_t c, long prec)
{
    long i, j, rows, cols;
    rows = acb_mat_nrows(A);
    cols = acb_mat_ncols(A);
    acb_t t;
    acb_init(t);

    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            acb_approx_mul(t,acb_mat_entry(A, i, j),c,prec);
            acb_approx_add(acb_mat_entry(res, i, j),t,acb_mat_entry(B, i, j),prec);

    acb_clear(t);
}

void //Computes c*A where c is a real scalar and A is a matrix
acb_mat_approx_scalar_mul_arb(acb_mat_t res, const acb_mat_t A, const arb_t c, long prec)
{
    long i, j, rows, cols;
    rows = acb_mat_nrows(A);
    cols = acb_mat_ncols(A);

    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            acb_approx_mul_arb(acb_mat_entry(res, i, j),acb_mat_entry(A, i, j),c,prec);
}

void //Computes x dot y where x&y are len*1 matrices.
acb_mat_approx_dot(acb_t res, acb_mat_t x, acb_mat_t y, long prec)
{
    int len = acb_mat_nrows(x);
    acb_ptr vec_x = _acb_vec_init(len);
    acb_ptr vec_y = _acb_vec_init(len);

    for (int i = 0; i < len; i++) //Add column entries to vector
    {
        acb_approx_set(vec_x + i,acb_mat_entry(x,i,0));
        acb_approx_set(vec_y + i,acb_mat_entry(y,i,0));
    }

    acb_approx_dot(res, 0, 0, vec_x, 1, vec_y, 1, len, prec);

    _acb_vec_clear(vec_x, len);
    _acb_vec_clear(vec_y, len);
}

void //Computes conjugate(x) dot y where x&y are len*1 matrices.
acb_mat_approx_dotc(acb_t res, acb_mat_t x, acb_mat_t y, long prec)
{
    int len = acb_mat_nrows(x);
    acb_ptr vec_x = _acb_vec_init(len);
    acb_ptr vec_y = _acb_vec_init(len);

    for (int i = 0; i < len; i++) //Add column entries to vector & conjugate x
    {
        acb_conj(vec_x + i,acb_mat_entry(x,i,0));
        acb_approx_set(vec_y + i,acb_mat_entry(y,i,0));
    }

    acb_approx_dot(res, 0, 0, vec_x, 1, vec_y, 1, len, prec);

    _acb_vec_clear(vec_x, len);
    _acb_vec_clear(vec_y, len);
}

void acb_mat_approx_norm(arb_t res, acb_mat_t x, long prec){
    acb_t tmp;
    acb_init(tmp);
    acb_mat_approx_dotc(tmp, x, x, prec); //There might be a faster way to compute this
    arb_set(res,acb_realref(tmp));
    arb_sqrt(res,res,prec);
    acb_clear(tmp);
}

void test_acb_mat_approx_dot()
    {
        acb_mat_t x, y;
        acb_mat_init(x, 3, 1);
        acb_mat_init(y, 3, 1);
        acb_mat_ones(x);
        acb_mat_ones(y);
        acb_onei(acb_mat_entry(x,0,0));
        acb_t res;
        acb_init(res);
        acb_mat_approx_dot(res, x, y, 53);
        acb_printd(res,10);
    }