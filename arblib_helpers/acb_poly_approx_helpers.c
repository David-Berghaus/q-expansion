#include "acb_poly.h"

void acb_poly_get_mid(acb_poly_t f, const acb_poly_t g)
{
    int i;

    for (i = 0; i < acb_poly_length(g); i++)
    {
        acb_get_mid(acb_poly_get_coeff_ptr(f,i), acb_poly_get_coeff_ptr(g,i));
    }
}

void //Return coefficients of f as acb_mat_t column-vector
acb_mat_set_poly(acb_mat_t res, const acb_poly_t f)
{
    int i;

    for (i = 0; i < acb_poly_length(f); i++)
    {
        acb_set(acb_mat_entry(res,i,0), acb_poly_get_coeff_ptr(f,i));
    }
}

void //Construct polynomial from acb_mat_t column-vector
acb_poly_set_mat(acb_poly_t f, const acb_mat_t c)
{
    int i;

    for (i = 0; i < acb_mat_nrows(c); i++)
    {
        acb_set(acb_poly_get_coeff_ptr(f,i), acb_mat_entry(c,i,0));
    }
}