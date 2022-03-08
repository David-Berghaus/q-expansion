#include "acb.h"
#include "acb_mat.h"
#include "acb_approx_helpers.h"

void
acb_mat_approx_add(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, int prec)
{
    int i, j, r, c;
    r = acb_mat_nrows(mat1);
    c = acb_mat_ncols(mat1);

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            acb_approx_add(acb_mat_entry(res, i, j), acb_mat_entry(mat1, i, j), acb_mat_entry(mat2, i, j), prec);
        }
    }
}

void
acb_mat_approx_sub(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, int prec)
{
    int i, j, r, c;
    r = acb_mat_nrows(mat1);
    c = acb_mat_ncols(mat1);

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            acb_approx_sub(acb_mat_entry(res, i, j), acb_mat_entry(mat1, i, j), acb_mat_entry(mat2, i, j), prec);
        }
    }
}

void //Computes A*c+B and stores it in res
acb_mat_approx_scalar_addmul(acb_mat_t res, acb_mat_t B, const acb_mat_t A, const acb_t c, int prec)
{
    int i, j, rows, cols;
    rows = acb_mat_nrows(A);
    cols = acb_mat_ncols(A);
    acb_t t;
    acb_init(t);

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            acb_approx_mul(t, acb_mat_entry(A,i,j), c, prec);
            acb_approx_add(acb_mat_entry(res,i,j), t, acb_mat_entry(B, i, j), prec);
        }
    }

    acb_clear(t);
}

void //Computes c*A where c is a complex scalar and A is a matrix
acb_mat_approx_scalar_mul(acb_mat_t res, const acb_mat_t A, const acb_t c, int prec)
{
    int i, j, rows, cols;
    rows = acb_mat_nrows(A);
    cols = acb_mat_ncols(A);

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            acb_approx_mul(acb_mat_entry(res, i, j),acb_mat_entry(A, i, j),c,prec);
        }
    }
}

void //Computes c*A where c is a real scalar and A is a matrix
acb_mat_approx_scalar_mul_arb(acb_mat_t res, const acb_mat_t A, const arb_t c, int prec)
{
    int i, j, rows, cols;
    rows = acb_mat_nrows(A);
    cols = acb_mat_ncols(A);

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            acb_approx_mul_arb(acb_mat_entry(res, i, j),acb_mat_entry(A, i, j),c,prec);
        }
    }
}

void //Computes conjugate(x) dot y where x&y are len*1 matrices.
acb_mat_approx_dotc(acb_t res, acb_mat_t x, acb_mat_t y, int prec)
{
    int len = acb_mat_nrows(x);
    arb_ptr tmp1, tmp2;
    TMP_INIT;

    TMP_START;
    tmp1 = TMP_ALLOC(sizeof(arb_struct) * len);
    tmp2 = TMP_ALLOC(sizeof(arb_struct) * len);

    int i;
    for (i = 0; i < len; i++)
    {
        tmp1[i] = *acb_realref(acb_mat_entry(x, i, 0));
        tmp2[i] = *acb_realref(acb_mat_entry(y, i, 0));
    }
    arb_approx_dot(acb_realref(res), NULL, 0, tmp1, 1, tmp2, 1, len, prec);

    for (i = 0; i < len; i++)
        tmp2[i] = *acb_imagref(acb_mat_entry(y, i, 0));
    arb_approx_dot(acb_imagref(res), NULL, 0, tmp1, 1, tmp2, 1, len, prec);

    for (i = 0; i < len; i++)
        tmp1[i] = *acb_imagref(acb_mat_entry(x, i, 0));
    arb_approx_dot(acb_realref(res), acb_realref(res), 0, tmp1, 1, tmp2, 1, len, prec);

    for (i = 0; i < len; i++)
        tmp2[i] = *acb_realref(acb_mat_entry(y, i, 0));
    arb_approx_dot(acb_imagref(res), acb_imagref(res), 1, tmp1, 1, tmp2, 1, len, prec);

    TMP_END;
}

void //Computes x dot y where x&y are len*1 matrices.
acb_mat_approx_dot(acb_t res, acb_mat_t x, acb_mat_t y, int prec)
{
    int len = acb_mat_nrows(x);
    acb_ptr tmp1, tmp2;
    TMP_INIT;

    TMP_START;
    tmp1 = TMP_ALLOC(sizeof(acb_struct) * len);
    tmp2 = TMP_ALLOC(sizeof(acb_struct) * len);

    int i;
    for (i = 0; i < len; i++)
    {
        tmp1[i] = *acb_mat_entry(x, i, 0);
        tmp2[i] = *acb_mat_entry(y, i, 0);
    }
    acb_approx_dot(res, NULL, 0, tmp1, 1, tmp2, 1, len, prec);

    TMP_END;
}

void acb_mat_approx_norm(arb_t res, acb_mat_t x, int prec){ //returns sqrt(real(dotc(x, x)))
    arb_ptr tmp;
    TMP_INIT;

    TMP_START;
    int i;
    int len = acb_mat_nrows(x);
    tmp = TMP_ALLOC(sizeof(arb_struct) * len);

    for (i = 0; i < len; i++)
        tmp[i] = *acb_realref(acb_mat_entry(x, i, 0));
    arb_approx_dot(res, NULL, 0, tmp, 1, tmp, 1, len, prec);

    for (i = 0; i < len; i++)
        tmp[i] = *acb_imagref(acb_mat_entry(x, i, 0));
    arb_approx_dot(res, res, 0, tmp, 1, tmp, 1, len, prec);

    arf_sqrt(arb_midref(res), arb_midref(res), prec, ARF_RND_DOWN);    

    TMP_END;
}

void acb_approx_complex_sign(acb_t res, acb_t z, arb_t z_abs, int prec)
{
    if (acb_is_zero(z) != 0)
        acb_one(res);
    else
        acb_approx_div_arb(res,z,z_abs,prec);
}

void //See Algorithm 1 of https://www.netlib.org/lapack/lawnspdf/lawn148.pdf
lartg(acb_t c, acb_t s, acb_t r, acb_t f, acb_t g, int prec)
{
    arb_t f_abs, g_abs, fg_norm, arb_tmp;
    acb_t f_sign;
    if (acb_is_zero(g) != 0)
    {
        // c = 1
        // s = 0
        // r = f
        acb_one(c);
        acb_zero(s);
        acb_approx_set(r,f);
    }
    else if (acb_is_zero(f) != 0)
    {
        // g_abs = abs(g)
        // c = 0
        // s = complex_sign(np.conjugate(g),g_abs)
        // r = g_abs
        arb_init(g_abs);
        acb_approx_abs(g_abs,g,prec);
        acb_zero(c);
        acb_conj(s,g);
        acb_approx_complex_sign(s,s,g_abs,prec);
        acb_approx_set_arb(r,g_abs);
        arb_clear(g_abs);
    }
    else
    {
        // f_abs, g_abs = abs(f), abs(g)
        // fg_norm = math.hypot(f_abs,g_abs)
        // f_sign = complex_sign(f,f_abs)
        // c = f_abs/fg_norm
        // s = f_sign*np.conjugate(g)/fg_norm
        // r = f_sign*fg_norm
        arb_init(f_abs);
        arb_init(g_abs);
        arb_init(fg_norm);
        arb_init(arb_tmp);
        acb_init(f_sign);
        acb_approx_abs(f_abs,f,prec);
        acb_approx_abs(g_abs,g,prec);
        arb_approx_hypot(fg_norm,f_abs,g_abs,prec);
        acb_approx_complex_sign(f_sign,f,f_abs,prec);
        arf_div(arb_midref(arb_tmp),arb_midref(f_abs),arb_midref(fg_norm),prec,ARF_RND_DOWN);
        acb_approx_set_arb(c,arb_tmp);
        acb_conj(s,g);
        acb_approx_mul(s,s,f_sign,prec);
        acb_approx_div_arb(s,s,fg_norm,prec);
        acb_approx_mul_arb(r,f_sign,fg_norm,prec);
        arb_clear(f_abs);
        arb_clear(g_abs);
        arb_clear(fg_norm);
        arb_clear(arb_tmp);
        acb_clear(f_sign);
    }
}

void acb_mat_change_prec(acb_mat_t res, acb_mat_t A, int prec) //Sets entries of A to prec (increases and decreases of precision are allowed)
{
    arf_t arf_tmp;
    arf_init(arf_tmp);
    int nrows, ncols, i, j;
    nrows = acb_mat_nrows(A);
    ncols = acb_mat_ncols(A);

    for (i = 0; i < nrows; i++)
    {
        for (j = 0; j < ncols; j++)
        {
            arf_set_round(arf_tmp, arb_midref(acb_realref(acb_mat_entry(A,i,j))), prec, ARF_RND_DOWN);
            arf_swap(arb_midref(acb_realref(acb_mat_entry(res,i,j))), arf_tmp);;
            arf_set_round(arf_tmp, arb_midref(acb_imagref(acb_mat_entry(A,i,j))), prec, ARF_RND_DOWN);
            arf_swap(arb_midref(acb_imagref(acb_mat_entry(res,i,j))), arf_tmp);
        }
    }

    arf_clear(arf_tmp);
}

void //Computes D*A where D is a diagonal matrix which is stored as a Nx1 vector
acb_mat_approx_left_mul_diag(acb_mat_t res, const acb_mat_t D, const acb_mat_t A, int prec)
{
    int i, j, rows, cols;
    rows = acb_mat_nrows(A);
    cols = acb_mat_ncols(A);

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            acb_approx_mul(acb_mat_entry(res, i, j), acb_mat_entry(A, i, j), acb_mat_entry(D, i, 0), prec);
        }
    }
}

void //Computes A*D where D is a diagonal matrix which is stored as a Nx1 vector
acb_mat_approx_right_mul_diag(acb_mat_t res, const acb_mat_t A, const acb_mat_t D, int prec)
{
    int i, j, rows, cols;
    rows = acb_mat_nrows(A);
    cols = acb_mat_ncols(A);

    for (i = 0; i < rows; i++)
    {
        for (j = 0; j < cols; j++)
        {
            acb_approx_mul(acb_mat_entry(res, i, j), acb_mat_entry(A, i, j), acb_mat_entry(D, j, 0), prec);
        }
    }
}

//This function is the analogue of https://github.com/fredrik-johansson/arb/blob/637582ad487b364493ac39d45cf46fa09d5dfefb/acb/vec_set_powers.c
//for approximate arithmetic and using matrices instead of arrays
void acb_mat_set_powers_approx(acb_mat_t xs, const acb_t x, int prec)
{
    int i;
    int len = acb_mat_ncols(xs);

    for (i = 0; i < len; i++)
    {
        if (i == 0)
            acb_one(acb_mat_entry(xs,0,i));
        else if (i == 1)
            acb_set_round(acb_mat_entry(xs,0,i), x, prec);
        else if (i % 2 == 0)
            acb_approx_mul(acb_mat_entry(xs,0,i), acb_mat_entry(xs,0,i/2), acb_mat_entry(xs,0,i/2), prec);
        else
            acb_approx_mul(acb_mat_entry(xs,0,i), acb_mat_entry(xs,0,i-1), x, prec);
    }
}

void evaluate_modular_splitting_polynomial(acb_t res, acb_mat_t coeffs, acb_mat_t xs, acb_mat_t ys, int j, int k, int Ms, int Mf, int bit_prec)
{
    int coeff_index, l, m;
    acb_ptr coeffs_vec, ys_vec, xs_vec, partial_results;
    coeffs_vec = _acb_vec_init(k);
    ys_vec = _acb_vec_init(k);
    xs_vec = _acb_vec_init(j);
    partial_results = _acb_vec_init(j);

    for (m = 0; m < k; m++){
        acb_swap(ys_vec+m,acb_mat_entry(ys,0,m));
    }
    for (l = 0; l < j; l++){
        acb_swap(xs_vec+l,acb_mat_entry(xs,0,l));
    }
    for (l = 0; l < j; l++){
        for (m = 0; m < k; m++){
            coeff_index = j*m+l;
            if (coeff_index >= Ms && coeff_index <= Mf){ //Otherwise coeffs are left to be zero
                acb_swap(coeffs_vec+m,acb_mat_entry(coeffs,coeff_index-Ms,0)); //Recall that the coeffs are stored in a Nx1 vector...
            }
        }
        acb_approx_dot(partial_results+l, NULL, 0, coeffs_vec, 1, ys_vec, 1, k, bit_prec);
        //Swap back
        for (m = 0; m < k; m++){
            coeff_index = j*m+l;
            if (coeff_index >= Ms && coeff_index <= Mf){ //Otherwise coeffs are left to be zero
                acb_swap(coeffs_vec+m,acb_mat_entry(coeffs,coeff_index-Ms,0)); //Recall that the coeffs are stored in a Nx1 vector...
            }
        }
    }
    acb_approx_dot(res, NULL, 0, partial_results, 1, xs_vec, 1, j, bit_prec);
    //Swap back
    for (m = 0; m < k; m++){
        acb_swap(ys_vec+m,acb_mat_entry(ys,0,m));
    }
    //Swap back
    for (l = 0; l < j; l++){
        acb_swap(xs_vec+l,acb_mat_entry(xs,0,l));
    }
    _acb_vec_clear(coeffs_vec, k);
    _acb_vec_clear(ys_vec, k);
    _acb_vec_clear(xs_vec, j);
    _acb_vec_clear(partial_results, j);
}