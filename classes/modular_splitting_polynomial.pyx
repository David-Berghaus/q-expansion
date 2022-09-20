# ****************************************************************************
#       Copyright (C) 2022 David Berghaus <berghaus@th.physik.uni-bonn.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off

from math import sqrt

from sage.rings.complex_arb cimport *
from sage.libs.arb.acb_mat cimport *
from sage.rings.real_arb import RealBallField
from sage.rings.complex_arb import ComplexBallField

from arblib_helpers.acb_approx cimport *

cdef class Modular_splitting_polynomial():
    """
    Class for storing and evaluating rows of W using modular splitting.
    """
    def __init__(self, int Ms, int Mf, ComplexBall q, RealBall y_fund_fact, int bit_prec):
        cdef CBF = ComplexBallField(bit_prec)
        cdef ComplexBall tmp = CBF(0)
        self.y_fund_fact = y_fund_fact
        coeff_len = Mf+1
        n_sqrt = int(sqrt(coeff_len))
        k = n_sqrt
        if k**2 != coeff_len: #coeff_len is not a square number
            k += 1
        j = n_sqrt
        self.xs, self.ys = Acb_Mat(1,j), Acb_Mat(1,k)
        acb_mat_set_powers_approx(self.xs.value,q.value,bit_prec)
        acb_approx_mul(tmp.value,acb_mat_entry(self.xs.value,0,j-1),q.value,bit_prec) #This is y = x^j
        acb_mat_set_powers_approx(self.ys.value,tmp.value,bit_prec)
        self.j, self.k, self.Ms, self.Mf = j, k, Ms, Mf

    cdef evaluate(self, acb_t res, Acb_Mat coeffs, int bit_prec):
        sig_on()
        evaluate_modular_splitting_polynomial(res,coeffs.value,self.xs.value,self.ys.value,self.j,self.k,self.Ms,self.Mf,bit_prec)
        acb_approx_mul_arb(res,res,self.y_fund_fact.value,bit_prec)
        sig_off()