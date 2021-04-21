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
from point_matching.point_matching_arb_wrap import digits_to_bits

def test_ir(S,int digit_prec,Y=0,int M=0):
    from psage.modform.maass.automorphic_forms_alg import get_M_for_holom   
    from point_matching.point_matching_arb_wrap import digits_to_bits, get_V_tilde_matrix_factored_b_arb_wrap
    import time
    bit_prec = digits_to_bits(digit_prec)
    RBF = RealBallField(bit_prec)
    CBF = ComplexBallField(bit_prec)
    if float(Y) == 0: #This comparison does not seem to be defined for arb-types...
        Y = RBF(S.group().minimal_height()*0.8)
    if M == 0:
        weight = S.weight()
        M = math.ceil(get_M_for_holom(Y,weight,digit_prec))
    print("Y = ", Y)
    print("M = ", M)
    print("dimen = ", S.group().ncusps()*M)
    cdef Block_Factored_Mat V
    cdef Acb_Mat b, V_inv, res
    cdef PLU_Mat plu
    V, b = get_V_tilde_matrix_factored_b_arb_wrap(S,M,Y,bit_prec)
    tol = RBF(10.0)**(-digit_prec)

    V_dp = V.construct_sc_np()
    plu = PLU_Mat(V_dp,prec=53)

    st = time.time()
    res = iterative_refinement_arb_wrap(V, b, bit_prec, tol, plu)
    print(time.time()-st)

    V.diag_inv_scale_vec(res, res, bit_prec)
    print("test result for Gamma0(1): ")
    acb_add_ui(acb_mat_entry(res.value,0,0), acb_mat_entry(res.value,0,0), 24, bit_prec)
    acb_printd(acb_mat_entry(res.value,0,0), digit_prec)
    print('')
    # res.str(10)


cpdef iterative_refinement_arb_wrap(Block_Factored_Mat A, Acb_Mat b, int prec, RealBall tol, PLU_Mat PLU, x0=None, maxiter=10**4, mix_prec=True, is_scaled=True):
    """Uses classical iterative refinement 
    to solve Ax = b
    Parameters
    ----------
    A : Block_Factored_Mat, linear system to solve
    b : Acb_Mat, right hand side
    x0 : array, matrix
        initial guess
    tol : float
        relative convergence tolerance
    maxiter : maximum amount of iterations
    PLU : matrix, approximate PLU-decomposition of matrix A.
    mix_prec : If True, gradually increases working precision, based on precision of residual
    scaled : If True, we solve a scaled version of the linear system (make sure to adapt your preconditioner
        to your choice of is_scaled)
    """
    cdef int i
    cdef int dimen = A.nrows()
    cdef int low_prec = 64
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
        acb_mat_approx_solve_lu_precomp(x.value, PLU.P, PLU.value, b.value, low_prec)

    for i in range(maxiter):
        if mix_prec == True:
            iter_prec = digits_to_bits(16*(i+2) + 5) #Gradually increase working precision during each iteration
        else:
            iter_prec = prec
        # r = b - A*x
        mat_vec_mul(r, A, x, iter_prec, is_scaled)
        sig_on()
        acb_mat_approx_sub(r.value, b.value, r.value, iter_prec)
        sig_off()
        # d = A\r
        acb_mat_approx_solve_lu_precomp(d.value, PLU.P, PLU.value, r.value, low_prec)
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
    
    arb_clear(normr)
    return x

    