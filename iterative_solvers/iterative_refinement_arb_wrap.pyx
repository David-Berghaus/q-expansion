# ****************************************************************************
#       Copyright (C) 2022 David Berghaus <berghaus@th.physik.uni-bonn.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_off
import math

from sage.misc.banner import version_dict
vers_dict = version_dict()
from sage.libs.arb.arb cimport *
from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_mat cimport *
from sage.rings.complex_arb cimport *
from sage.matrix.matrix_complex_ball_dense cimport *
from sage.rings.real_arb import RealBallField
from sage.rings.complex_arb import ComplexBallField
if vers_dict['major'] == 9 and vers_dict['minor'] == 2: #We still need to support sage 9.2 for now
    from sage.rings.complex_field import ComplexField
else:
    from sage.rings.complex_mpfr import ComplexField

from arblib_helpers.acb_approx cimport *
from classes.acb_mat_class cimport Acb_Mat, Acb_Mat_Win
from classes.plu_class cimport PLU_Mat
from classes.block_factored_mat_class cimport Block_Factored_Mat
from iterative_solvers.gmres_arb_wrap cimport mat_vec_mul
from belyi.number_fields import get_decimal_digit_prec
from point_matching.point_matching_arb_wrap import digits_to_bits

cpdef iterative_refinement_arb_wrap(Block_Factored_Mat A, Acb_Mat b, int prec, RealBall tol, PLU_Mat PLU, x0=None, maxiter=None, mix_prec=True, starting_prec=0, is_scaled=True, imposed_zeros=[]):
    """Uses classical iterative refinement 
    to solve Ax = b
    Parameters
    ----------
    A : Block_Factored_Mat, linear system to solve
    b : Acb_Mat, right hand side
    x0 : array, matrix
        initial guess (or result from restarted computation)
    tol : float
        relative convergence tolerance
    maxiter : maximum amount of iterations
    PLU : matrix, approximate PLU-decomposition of matrix A.
    mix_prec : If True, gradually increases working precision, based on precision of residual
    starting_prec: int, specify (bit) precision of first iteration, only use this for restarting
    scaled : If True, we solve a scaled version of the linear system (make sure to adapt your preconditioner
        to your choice of is_scaled)
    """
    cdef int i
    cdef int dimen = A.nrows()
    cdef int low_prec = 64
    cdef int low_digit_prec = 16
    #We use these variables to cast python objects to C-classes on which we can call native acb-functions
    #This conversion is a bit tedious but should not be inefficient
    cdef Acb_Mat acb_mat_cast
    cdef arb_t normr
    arb_init(normr)
    RBF = RealBallField(53)
    CC = ComplexField(53)
    cdef RealBall rbf_normr, rbf_normr_old
    rbf_normr, rbf_normr_old = RBF(0), RBF(0) #we use these expressions to print values properly in python

    cdef Acb_Mat r = Acb_Mat(dimen,1)
    cdef Acb_Mat d = Acb_Mat(dimen,1)
    cdef Acb_Mat x = Acb_Mat(dimen,1)
    if x0 != None:
        acb_mat_cast = x0
        acb_mat_set(x.value, acb_mat_cast.value)
    else:
        PLU.solve(x,b,low_prec,imposed_zeros=imposed_zeros)
    if maxiter == None:
        maxiter = 2147483647

    for i in range(maxiter):
        if mix_prec == True:
            iter_prec = digits_to_bits(low_digit_prec*(i+2) + 5 + starting_prec) #Gradually increase working precision during each iteration
        else:
            iter_prec = prec
        # r = b - A*x
        mat_vec_mul(r, A, x, iter_prec, is_scaled, imposed_zeros=imposed_zeros)
        sig_on()
        acb_mat_approx_sub(r.value, b.value, r.value, iter_prec)
        sig_off()
        # d = A\r
        PLU.solve(d,r,low_prec,imposed_zeros=imposed_zeros)
        # x = x + d
        acb_mat_approx_add(x.value, x.value, d.value, iter_prec)
        #normr = norm(r)
        sig_on()
        acb_mat_approx_norm(normr, r.value, low_prec)
        sig_off()
        if i != 0:
            arb_set(rbf_normr_old.value,rbf_normr.value)
        arb_set(rbf_normr.value,normr)
        # print(str(i) + ".iteration: " + str(CC(rbf_normr)))
        print(str(i) + ", " + str(CC(rbf_normr)))
        # if normr < tol:
        if arb_lt(normr, tol.value) == 1:
            break
        if i == maxiter-1:
            raise ArithmeticError("Maximum amount of iterations reached without sufficient convergence!")
        if i != 0 and (rbf_normr_old-rbf_normr)/(rbf_normr_old+rbf_normr) < RBF(1e-3):
            print("Convergence is too slow, this probably happens because we cannot impose a Victor Miller normalization on the considered form!")
            raise ArithmeticError("Convergence is too slow, this probably happens because we cannot impose a Victor Miller normalization on the considered form!")
    
    arb_clear(normr)
    return x

    