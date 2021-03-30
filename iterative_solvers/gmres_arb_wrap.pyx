#This is an implementation of GMRES MGS for acb-matrices. 
#It is based on: https://github.com/pyamg/pyamg/blob/main/pyamg/krylov/_gmres_mgs.py

import math
from cysignals.signals cimport sig_on, sig_off

from sage.libs.arb.arb cimport *
from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_mat cimport *
from sage.rings.complex_arb cimport *
from sage.matrix.matrix_complex_ball_dense cimport *
from sage.rings.real_arb import RealBallField
from sage.rings.complex_arb import ComplexBallField
from sage.matrix.matrix_space import MatrixSpace

from arblib_helpers.acb_mat_approx cimport *
from arblib_helpers.acb_mat_approx_helpers cimport *
from classes.acb_mat_class cimport Acb_Mat, Acb_Mat_Win
from classes.plu_class cimport PLU_Mat

def test_gmres(S,int digit_prec,Y=0,int M=0):
    from psage.modform.maass.automorphic_forms_alg import get_M_for_holom   
    from point_matching.point_matching_arb_wrap import digits_to_bits, get_V_tilde_matrix_b_arb_wrap
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
    cdef Matrix_complex_ball_dense V, b
    V, b = get_V_tilde_matrix_b_arb_wrap(S,M,Y,bit_prec)
    dimen = acb_mat_nrows(V.value)
    cdef Acb_Mat V2, V2_inv, b2, res, res2
    V2 = Acb_Mat(dimen, dimen)
    V2_inv = Acb_Mat(dimen, dimen)
    V2._set_mcbd(V)
    b2, res, res2 = Acb_Mat(dimen, 1), Acb_Mat(dimen, 1), Acb_Mat(dimen, 1)
    b2._set_mcbd(b)

    tol = RBF(10.0)**(-digit_prec)
    low_prec = 64
    acb_mat_approx_inv(V2_inv.value, V2.value, low_prec)
    # cdef ComplexBall pert = CBF(1e-16+1e-16j)
    # for i in range(M):
    #     for j in range(M):
    #         acb_approx_sub(acb_mat_entry(V2_inv.value,i,j), acb_mat_entry(V2_inv.value,i,j), pert.value, inv_prec)
    # acb_mat_change_prec(V2_inv.value, V2_inv.value, bit_prec)
    x_gmres_arb_wrap = gmres_mgs_arb_wrap(V2, b2, res2, bit_prec, tol, M=V2_inv)
    res2 = x_gmres_arb_wrap[0]
    print("test result for Gamma0(1): ")
    acb_add_ui(acb_mat_entry(res2.value,0,0), acb_mat_entry(res2.value,0,0), 24, bit_prec)
    acb_printd(acb_mat_entry(res2.value,0,0), digit_prec)
    # print('')
    # x_gmres_arb_wrap[0].str(10)

cdef apply_givens(list Q, Acb_Mat_Win v, int k, int prec):
    """Apply the first k Givens rotations in Q to v.
    Parameters
    ----------
    Q : list
        list of consecutive 2x2 Givens rotations
    v : array
        vector to apply the rotations to
    k : int
        number of rotations to apply.
    Returns
    -------
    v is changed in place
    Notes
    -----
    This routine is specialized for GMRES.  It assumes that the first Givens
    rotation is for dofs 0 and 1, the second Givens rotation is for
    dofs 1 and 2, and so on.
    """
    cdef acb_t v1, v2, t
    acb_init(v1)
    acb_init(v2)
    acb_init(t)
    cdef int j
    cdef Acb_Mat Qloc
    for j in range(k):
        Qloc = Q[j]
        #We do the matrix multiplications "by hand" because it should be faster than calling acb_mat_mul_approx
        #for 2x2 matrices because we avoid some allocations
        # v1 = Qloc[0,0]*v[j]+Qloc[0,1]*v[j+1]
        acb_approx_mul(t,acb_mat_entry(Qloc.value,0,0),acb_mat_entry(v.value,j,0),prec)
        acb_approx_mul(v1,acb_mat_entry(Qloc.value,0,1),acb_mat_entry(v.value,j+1,0),prec)
        acb_approx_add(v1,v1,t,prec)
        # v2 = Qloc[1,0]*v[j]+Qloc[1,1]*v[j+1]
        acb_approx_mul(t,acb_mat_entry(Qloc.value,1,0),acb_mat_entry(v.value,j,0),prec)
        acb_approx_mul(v2,acb_mat_entry(Qloc.value,1,1),acb_mat_entry(v.value,j+1,0),prec)
        acb_approx_add(v2,v2,t,prec)
        # v[j] = v1
        # v[j+1] = v2
        acb_swap(acb_mat_entry(v.value,j,0),v1)
        acb_swap(acb_mat_entry(v.value,j+1,0),v2)
    acb_clear(v1)
    acb_clear(v2)
    acb_clear(t)

cpdef gmres_mgs_arb_wrap(Acb_Mat A, Acb_Mat b, Acb_Mat x0, int prec, RealBall tol, restrt=None, maxiter=None, M=None, PLU=None):
    """Generalized Minimum Residual Method (GMRES) based on MGS.
    GMRES iteratively refines the initial solution guess to the system
    Ax = b
    Modified Gram-Schmidt version
    Parameters
    ----------
    A : array, matrix, sparse matrix, LinearOperator
        n x n, linear system to solve
    b : array, matrix
        right hand side, shape is (n,) or (n,1)
    x0 : array, matrix
        initial guess
    tol : float
        relative convergence tolerance, i.e. tol is scaled by the norm
        of the initial preconditioned residual
    restrt : None, int
        - if int, restrt is max number of inner iterations
          and maxiter is the max number of outer iterations
        - if None, do not restart GMRES, and max number of inner iterations
          is maxiter
    maxiter : None, int
        - if restrt is None, maxiter is the max number of inner iterations
          and GMRES does not restart
        - if restrt is int, maxiter is the max number of outer iterations,
          and restrt is the max number of inner iterations
    M : matrix, inverted preconditioner, i.e. solve M A x = M b.
    PLU : matrix, approximate PLU-decomposition of matrix A. (as it gets computed by arb)
    Returns
    -------
    (xNew, info)
    xNew : an updated guess to the solution of Ax = b
    info : halting status of gmres
            ==  =============================================
            0   successful exit
            >0  convergence to tolerance not achieved,
                return iteration count instead.  This value
                is precisely the order of the Krylov space.
            <0  numerical breakdown, or illegal input
            ==  =============================================
    Notes
    -----
        - For robustness, modified Gram-Schmidt is used to orthogonalize the
          Krylov Space Givens Rotations are used to provide the residual norm
          each iteration
    References
    ----------
    .. [1] Yousef Saad, "Iterative Methods for Sparse Linear Systems,
       Second Edition", SIAM, pp. 151-172, pp. 272-275, 2003
       http://www-users.cs.umn.edu/~saad/books.html
    .. [2] C. T. Kelley, http://www4.ncsu.edu/~ctk/matlab_roots.html
    """
    cdef int outer, k, niter, max_outer, max_inner
    cdef int inner = 0 #To avoid compiler warning
    cdef int dimen = acb_mat_nrows(A.value)

    cdef Acb_Mat x = Acb_Mat(dimen,1)
    acb_mat_set(x.value, x0.value)

    cdef arb_t normr, normv, arb_tmp
    arb_init(normr)
    arb_init(normv)
    arb_init(arb_tmp)
    cdef acb_t c, s, acb_tmp, acb_tmp2, acb_tmp3
    acb_init(c)
    acb_init(s)
    acb_init(acb_tmp)
    acb_init(acb_tmp2)
    acb_init(acb_tmp3)

    # Set number of outer and inner iterations
    if restrt:
        if maxiter:
            max_outer = maxiter
        else:
            max_outer = 1
        if restrt > dimen:
            print('Setting number of inner iterations (restrt) to maximum\
                  allowed, which is A.shape[0] ')
            restrt = dimen
        max_inner = restrt
    else:
        max_outer = 1
        if maxiter > dimen:
            print('Setting number of inner iterations (maxiter) to maximum\
                  allowed, which is A.shape[0] ')
            maxiter = dimen
        elif maxiter is None:
            maxiter = min(dimen, 40)
        max_inner = maxiter

    # Prep for method
    # r = b - np.ravel(A*x)
    cdef Acb_Mat r = Acb_Mat(dimen,1)
    sig_on()
    acb_mat_approx_mul(r.value,A.value,x.value,prec)
    sig_off()
    sig_on()
    acb_mat_approx_sub(r.value,b.value,r.value,prec)
    sig_off()

    # Apply preconditioner
    #r = np.ravel(M*r)
    apply_preconditioner(r, M, PLU, prec)

    #normr = norm(r)
    sig_on()
    acb_mat_approx_norm(normr,r.value,prec)
    sig_off()

    # Scale tol by ||r_0||_2, we use the preconditioned residual
    # because this is left preconditioned GMRES.
    # if normr != 0.0:
    #     tol = tol*normr
    # if arb_is_zero(normr) == 0:
    #     arb_mul(tol.value, tol.value, normr, prec)

    # Use separate variable to track iterations.  If convergence fails, we
    # cannot simply report niter = (outer-1)*max_outer + inner.  Numerical
    # error could cause the inner loop to halt while the actual ||r|| > tol.
    niter = 0

    cdef Acb_Mat H, V, g, Qblock
    cdef Acb_Mat_Win v

    #We use these variables to cast python objects to C-classes on which we can call native acb-functions
    #This conversion is a bit tedious but should not be inefficient
    cdef Acb_Mat acb_mat_cast
    cdef Acb_Mat_Win acb_mat_win_cast, acb_mat_win_cast2

    # Begin GMRES
    for outer in range(max_outer):

        # Preallocate for Givens Rotations, Hessenberg matrix and Krylov Space
        # Space required is O(dimen*max_inner).
        Q = []  # Givens Rotations -> one might want to use libcpp.vector but the performance impact shouldn't be noticeable
        # Upper Hessenberg matrix, which is then
        #   converted to upper tri with Givens Rots
        H = Acb_Mat(max_inner+1, max_inner+1)
        V = Acb_Mat(dimen, max_inner+1)  # Krylov Space
        # vs stores the windows to each column of V.
        vs = []
        # v = r/normr
        # V[:, 0] = scal(1.0/normr, r)
        # vs.append(V[:, 0])
        vs.append(V.get_window(0,0,dimen,1))
        arb_inv(arb_tmp,normr,prec)
        v = vs[0]
        sig_on()
        acb_mat_approx_scalar_mul_arb(v.value, r.value, arb_tmp, prec)
        sig_off()

        # This is the RHS vector for the problem in the Krylov Space
        # g = np.zeros((dimen,), dtype=complex)
        # g[0] = normr
        g = Acb_Mat(dimen, 1)
        acb_approx_set_arb(acb_mat_entry(g.value,0,0), normr)

        for inner in range(max_inner):

            # New Search Direction
            # v = V[:, inner+1]
            # v[:] = np.ravel(M*(A*vs[-1]))
            # vs.append(v)
            # normv_old = norm(v)
            vs.append(V.get_window(0,inner+1,dimen,inner+2))
            v = vs[-1]
            acb_mat_win_cast = vs[-2]
            sig_on()
            acb_mat_approx_mul(v.value, A.value, acb_mat_win_cast.value, prec)
            sig_off()
            apply_preconditioner(v, M, PLU, prec)

            #  Modified Gram Schmidt
            for k in range(inner+1):
                # vk = vs[k]
                # alpha = dotc(vk, v)
                # H[k, inner] = alpha
                # v[:] = axpy(vk, v, dimen, -alpha)
                acb_mat_win_cast = vs[k]
                acb_mat_approx_dotc(acb_tmp, acb_mat_win_cast.value, v.value, prec)
                acb_approx_set(acb_mat_entry(H.value, k, inner), acb_tmp)
                acb_neg(acb_tmp, acb_tmp)
                acb_mat_approx_scalar_addmul(v.value, v.value, acb_mat_win_cast.value, acb_tmp, prec)

            # normv = norm(v)
            sig_on()
            acb_mat_approx_norm(normv,v.value,prec)
            sig_off()
            # H[inner+1, inner] = normv
            acb_approx_set_arb(acb_mat_entry(H.value, inner+1, inner), normv)

            # Check for breakdown
            # if H[inner+1, inner] != 0.0:
            #     v[:] = scal(1.0/H[inner+1, inner], v)
            if acb_is_zero(acb_mat_entry(H.value, inner+1, inner)) == 0:
                acb_inv(acb_tmp, acb_mat_entry(H.value, inner+1, inner), prec)
                sig_on()
                acb_mat_approx_scalar_mul(v.value, v.value, acb_tmp, prec)
                sig_off()

            # Apply previous Givens rotations to H
            if inner > 0:
                # apply_givens(Q, H[:, inner], inner)
                apply_givens(Q, H.get_window(0,inner,dimen,inner+1), inner, prec)

            # Calculate and apply next complex-valued Givens Rotation
            # ==> Note that if max_inner = dimen, then this is unnecessary
            # for the last inner
            #     iteration, when inner = dimen-1.
            if inner != dimen-1:
                # if H[inner+1, inner] != 0:
                if acb_is_zero(acb_mat_entry(H.value,inner+1,inner)) == 0:
                    # [c, s, r] = lartg(H[inner, inner], H[inner+1, inner])
                    # Qblock = np.array([[c, s], [-np.conjugate(s), c]],
                    #                   dtype=complex)
                    # Q.append(Qblock)
                    sig_on()
                    lartg(c, s, acb_tmp, acb_mat_entry(H.value,inner,inner), acb_mat_entry(H.value,inner+1,inner), prec)
                    sig_off()
                    Q.append(Acb_Mat(2,2))
                    Qblock = Q[-1]
                    acb_approx_set(acb_mat_entry(Qblock.value,0,0), c)
                    acb_approx_set(acb_mat_entry(Qblock.value,0,1), s)
                    acb_conj(acb_mat_entry(Qblock.value,1,0), s)
                    acb_neg(acb_mat_entry(Qblock.value,1,0), acb_mat_entry(Qblock.value,1,0))
                    acb_approx_set(acb_mat_entry(Qblock.value,1,1), c)

                    # Apply Givens Rotation to g,
                    #   the RHS for the linear system in the Krylov Subspace.
                    # g[inner:inner+2] = np.dot(Qblock, g[inner:inner+2])
                    acb_approx_mul(acb_tmp,acb_mat_entry(Qblock.value,0,0),acb_mat_entry(g.value,inner,0),prec)
                    acb_approx_mul(acb_tmp2,acb_mat_entry(Qblock.value,0,1),acb_mat_entry(g.value,inner+1,0),prec)
                    acb_approx_add(acb_tmp3,acb_tmp,acb_tmp2,prec)
                    acb_approx_mul(acb_tmp,acb_mat_entry(Qblock.value,1,0),acb_mat_entry(g.value,inner,0),prec)
                    acb_approx_mul(acb_tmp2,acb_mat_entry(Qblock.value,1,1),acb_mat_entry(g.value,inner+1,0),prec)
                    acb_approx_add(acb_mat_entry(g.value,inner+1,0),acb_tmp,acb_tmp2,prec)
                    acb_swap(acb_mat_entry(g.value,inner,0),acb_tmp3)

                    # Apply effect of Givens Rotation to H
                    # H[inner, inner] = dotu(Qblock[0, :],
                    #                        H[inner:inner+2, inner])
                    # H[inner+1, inner] = 0.0
                    acb_approx_mul(acb_tmp, acb_mat_entry(Qblock.value,0,0), acb_mat_entry(H.value,inner,inner), prec)
                    acb_approx_mul(acb_tmp2, acb_mat_entry(Qblock.value,0,1), acb_mat_entry(H.value,inner+1,inner), prec)
                    acb_approx_add(acb_mat_entry(H.value,inner,inner), acb_tmp, acb_tmp2, prec)
                    acb_zero(acb_mat_entry(H.value,inner+1,inner))

            niter += 1

            # Don't update normr if last inner iteration, because
            # normr is calculated directly after this loop ends.
            if inner < max_inner-1:
                # normr = np.abs(g[inner+1])
                acb_approx_abs(normr, acb_mat_entry(g.value,inner+1,0), prec)
                arb_printd(normr, 10)
                print('')
                # if normr < tol:
                if arb_lt(normr, tol.value) == 1:
                    break

        # end inner loop, back to outer loop

        # Find best update to x in Krylov Space V.  Solve inner x inner system.
        # y = sp.linalg.solve(H[0:inner+1, 0:inner+1], g[0:inner+1]) #H[0:inner+1, 0:inner+1] is upper triangular
        acb_mat_cast = Acb_Mat(inner+1,1) #This is y
        acb_mat_win_cast = H.get_window(0, 0, inner+1, inner+1)
        acb_mat_win_cast2 = g.get_window(0, 0, inner+1, 1)
        sig_on()
        acb_mat_approx_solve_triu(acb_mat_cast.value, acb_mat_win_cast.value, acb_mat_win_cast2.value, 0, prec)
        sig_off()
        # update = np.ravel(V[:, :inner+1].dot(y.reshape(-1, 1)))
        acb_mat_win_cast = V.get_window(0, 0, dimen, inner+1)
        sig_on()
        acb_mat_approx_mul(acb_mat_cast.value, acb_mat_win_cast.value, acb_mat_cast.value, prec) #acb_mat_cast is now update
        sig_off()
        # x = x + update
        sig_on()
        acb_mat_approx_add(x.value, x.value, acb_mat_cast.value, prec)
        sig_off()
        # r = b - np.ravel(A*x)
        sig_on()
        acb_mat_approx_mul(r.value, A.value, x.value, prec)
        sig_off()
        sig_on()
        acb_mat_approx_sub(r.value, b.value, r.value, prec)
        sig_off()

        # Apply preconditioner
        # r = np.ravel(M*r)
        apply_preconditioner(r, M, PLU, prec)
        # normr = norm(r)
        sig_on()
        acb_mat_approx_norm(normr,r.value,prec)
        sig_off()

        # test for convergence
        if arb_lt(normr, tol.value) == 1:
            return (x, 0)

    # end outer loop
    arb_clear(normr)
    arb_clear(normv)
    arb_clear(arb_tmp)
    acb_clear(c)
    acb_clear(s)
    acb_clear(acb_tmp)
    acb_clear(acb_tmp2)
    acb_clear(acb_tmp3)

    return (x, niter)

cdef apply_preconditioner(x, M, PLU, int prec): #Apply preconditioner (if available) on x
    """
    Parameters
    ----------
    x : vector / matrix (view)
    M : preconditioner matrix, i.e. sets x=M*x
    PLU : approximate LU-decomposition of matrix A, i.e. sets x=(P*L*U)^(-1)*x
    """

    cdef Acb_Mat M_cast, x_cast
    cdef Acb_Mat_Win x_win_cast #Use this in case x is a window
    cdef acb_mat_t aliasing_avoider_matrix
    cdef PLU_Mat PLU_cast

    if M != None and PLU != None:
        raise NameError("Please select only one preconditioner!")
        
    if M != None:
        M_cast = M
        if type(x) is Acb_Mat:
            x_cast = x
            sig_on()
            acb_mat_approx_mul(x_cast.value, M_cast.value, x_cast.value, prec)
            sig_off()
        elif type(x) is Acb_Mat_Win:
            x_win_cast = x
            dimen = acb_mat_nrows(x_win_cast.value)
            acb_mat_init(aliasing_avoider_matrix, dimen, 1) #Arb currently does not work for aliasing window-matrices
            acb_mat_set(aliasing_avoider_matrix, x_win_cast.value)
            sig_on()
            acb_mat_approx_mul(x_win_cast.value, M_cast.value, aliasing_avoider_matrix, prec)
            sig_off()
            acb_mat_clear(aliasing_avoider_matrix)
    
    if PLU != None:
        PLU_cast = PLU
        if type(x) is Acb_Mat:
            x_cast = x
            sig_on()
            acb_mat_approx_solve_lu_precomp(x_cast.value, PLU_cast.P, PLU_cast.value, x_cast.value, prec)
            sig_off()
        elif type(x) is Acb_Mat_Win:
            x_win_cast = x
            sig_on()
            acb_mat_approx_solve_lu_precomp(x_win_cast.value, PLU_cast.P, PLU_cast.value, x_win_cast.value, prec)
            sig_off()