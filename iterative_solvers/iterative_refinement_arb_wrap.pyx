from cysignals.signals cimport sig_on, sig_off
import math

from sage.libs.arb.arb cimport *
from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_mat cimport *
from sage.rings.complex_arb cimport *
from sage.matrix.matrix_complex_ball_dense cimport *
from sage.rings.real_arb import RealBallField
from sage.rings.complex_arb import ComplexBallField

from arblib_helpers.acb_approx cimport *
from classes.acb_mat_class cimport Acb_Mat, Acb_Mat_Win
from classes.plu_class cimport PLU_Mat
from classes.block_factored_mat_class cimport Block_Factored_Mat
from iterative_solvers.gmres_arb_wrap cimport mat_vec_mul
from belyi.number_fields import get_decimal_digit_prec
from point_matching.point_matching_arb_wrap import digits_to_bits

cpdef iterative_refinement_arb_wrap(Block_Factored_Mat A, Acb_Mat b, int prec, RealBall tol, PLU_Mat PLU, x0=None, maxiter=None, mix_prec=True, starting_prec=0, is_scaled=True):
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

    cdef Acb_Mat r = Acb_Mat(dimen,1)
    cdef Acb_Mat d = Acb_Mat(dimen,1)
    cdef Acb_Mat x = Acb_Mat(dimen,1)
    if x0 != None:
        acb_mat_cast = x0
        acb_mat_set(x.value, acb_mat_cast.value)
    else:
        PLU.solve(x,b,low_prec)
    if maxiter == None:
        maxiter = get_decimal_digit_prec(tol)//15 + 10 #We expect to get around 15 digits per iteration, otherwise something likely went wrong

    for i in range(maxiter):
        if mix_prec == True:
            iter_prec = digits_to_bits(low_digit_prec*(i+2) + 5 + starting_prec) #Gradually increase working precision during each iteration
        else:
            iter_prec = prec
        # r = b - A*x
        mat_vec_mul(r, A, x, iter_prec, is_scaled)
        sig_on()
        acb_mat_approx_sub(r.value, b.value, r.value, iter_prec)
        sig_off()
        # d = A\r
        PLU.solve(d,r,low_prec)
        # x = x + d
        acb_mat_approx_add(x.value, x.value, d.value, iter_prec)
        #normr = norm(r)
        sig_on()
        acb_mat_approx_norm(normr, r.value, low_prec)
        sig_off()
        arb_printd(normr, 10)
        print(' (' + str(i) + '. iteration)')
        # if normr < tol:
        if arb_lt(normr, tol.value) == 1:
            break
        if i == maxiter-1:
            raise ArithmeticError("Maximum amount of iterations reached without sufficient convergence!")
    
    arb_clear(normr)
    return x

    