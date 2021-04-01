//We declare this as a header file because we have *.c in gitignore...

#include "acb.h"

void
acb_approx_mul(acb_t res, const acb_t x, const acb_t y, slong prec)
{
    arf_complex_mul(arb_midref(acb_realref(res)), arb_midref(acb_imagref(res)),
        arb_midref(acb_realref(x)), arb_midref(acb_imagref(x)), 
        arb_midref(acb_realref(y)), arb_midref(acb_imagref(y)), prec, ARF_RND_DOWN);
}

void
acb_approx_add(acb_t res, const acb_t x, const acb_t y, slong prec)
{
    arf_add(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)), arb_midref(acb_realref(y)), prec, ARF_RND_DOWN);
    arf_add(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)), arb_midref(acb_imagref(y)), prec, ARF_RND_DOWN);
}

void
acb_approx_sub(acb_t res, const acb_t x, const acb_t y, slong prec)
{
    arf_sub(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)), arb_midref(acb_realref(y)), prec, ARF_RND_DOWN);
    arf_sub(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)), arb_midref(acb_imagref(y)), prec, ARF_RND_DOWN);
}

void
acb_approx_set(acb_t res, const acb_t x)
{
    arf_set(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)));
    arf_set(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)));
}

void
acb_approx_set_arb(acb_t res, const arb_t x)
{
    arf_set(arb_midref(acb_realref(res)), arb_midref(x));
    arf_zero(arb_midref(acb_imagref(res)));
}

void
acb_approx_div_arb(acb_t res, const acb_t x, const arb_t y, slong prec)
{
    arf_div(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)), arb_midref(y), prec, ARF_RND_DOWN);
    arf_div(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)), arb_midref(y), prec, ARF_RND_DOWN);
}

void
acb_approx_inv(acb_t z, const acb_t x, slong prec)
{
    arf_set(arb_midref(acb_realref(z)), arb_midref(acb_realref(x)));
    arf_set(arb_midref(acb_imagref(z)), arb_midref(acb_imagref(x)));

    mag_zero(arb_radref(acb_realref(z)));
    mag_zero(arb_radref(acb_imagref(z)));

    acb_inv(z, z, prec);

    mag_zero(arb_radref(acb_realref(z)));
    mag_zero(arb_radref(acb_imagref(z)));
}

void
acb_approx_div(acb_t z, const acb_t x, const acb_t y, slong prec)
{
    acb_t t;
    acb_init(t);
    acb_approx_inv(t, y, prec);
    acb_approx_mul(z, x, t, prec);
    acb_clear(t);
}

void
acb_approx_mul_arb(acb_t res, const acb_t x, const arb_t y, slong prec)
{
    arf_mul(arb_midref(acb_realref(res)), arb_midref(acb_realref(x)), arb_midref(y), prec, ARF_RND_DOWN);
    arf_mul(arb_midref(acb_imagref(res)), arb_midref(acb_imagref(x)), arb_midref(y), prec, ARF_RND_DOWN);
}

void
arb_approx_hypot(arb_t z, const arb_t x, const arb_t y, slong prec)
{
    if (arb_is_zero(y))
    {
        arb_abs(z, x);
    }
    else if (arb_is_zero(x))
    {
        arb_abs(z, y);
    }
    else
    {
        arb_t t;
        arb_init(t);
        /* TODO: use arb_fmma */
        arf_mul(arb_midref(t), arb_midref(x), arb_midref(x), prec + 4, ARF_RND_DOWN);
        arf_mul(arb_midref(z), arb_midref(y), arb_midref(y), prec + 4, ARF_RND_DOWN);
        arf_add(arb_midref(t), arb_midref(t), arb_midref(z), prec + 4, ARF_RND_DOWN);
        arf_sqrt(arb_midref(z), arb_midref(t), prec, ARF_RND_DOWN);
        arb_clear(t);
    }
}

void acb_approx_abs(arb_t r, const acb_t z, slong prec)
{
    arb_approx_hypot(r, acb_realref(z), acb_imagref(z), prec);
}