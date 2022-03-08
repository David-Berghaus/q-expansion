cimport numpy as np

from sage.rings.complex_arb cimport *

from arblib_helpers.acb_approx cimport *
from classes.acb_mat_class cimport Acb_Mat

cdef class Modular_splitting_polynomial():
    cdef int j, k, Ms, Mf
    cdef Acb_Mat xs, ys #Python cannot handle pointers so we store the arrays as matrices: https://stackoverflow.com/a/47008303
    cdef RealBall y_fund_fact

    cdef evaluate(self, acb_t res, Acb_Mat coeffs, int bit_prec)