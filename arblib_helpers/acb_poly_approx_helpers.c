/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson
    Copyright (C) 2021 David Berghaus
    This file is providing functionality for approximate polynomials, based on arb.
    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_approx_helpers.h"

void
_acb_poly_approx_evaluate_horner(acb_t y, acb_srcptr f, int len,
                           const acb_t x, int prec)
{
    if (len == 0)
    {
        acb_zero(y);
    }
    else if (len == 1 || acb_is_zero(x))
    {
        acb_set_round(y, f, prec);
    }
    else if (len == 2)
    {
        acb_approx_mul(y, x, f + 1, prec);
        acb_approx_add(y, y, f + 0, prec);
    }
    else
    {
        int i = len - 1;
        acb_t t, u;

        acb_init(t);
        acb_init(u);
        acb_set(u, f + i);

        for (i = len - 2; i >= 0; i--)
        {
            acb_approx_mul(t, u, x, prec);
            acb_approx_add(u, f + i, t, prec);
        }

        acb_swap(y, u);

        acb_clear(t);
        acb_clear(u);
    }
}

void
acb_poly_approx_evaluate_horner(acb_t res, const acb_poly_t f, const acb_t a, int prec)
{
    _acb_poly_approx_evaluate_horner(res, f->coeffs, f->length, a, prec);
}