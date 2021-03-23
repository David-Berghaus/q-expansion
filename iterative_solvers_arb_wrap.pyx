#Based on: https://github.com/pyamg/pyamg/blob/main/pyamg/krylov/_gmres_mgs.py

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

from acb_mat_approx cimport *
from acb_mat_class cimport Acb_Mat, Acb_Mat_Win
from acb_mat_class import Acb_Mat, Acb_Mat_Win

cdef extern from "acb_mat_approx_helpers.h":
    void acb_approx_mul(acb_t res, const acb_t x, const acb_t y, long prec)
    void acb_approx_add(acb_t res, const acb_t x, const acb_t y, long prec)
    void acb_approx_sub(acb_t res, const acb_t x, const acb_t y, long prec)
    void acb_approx_set(acb_t res, const acb_t x)
    void acb_approx_set_arb(acb_t res, const arb_t x)
    void acb_approx_div_arb(acb_t res, const acb_t x, const arb_t y, long prec)
    void acb_approx_mul_arb(acb_t res, const acb_t x, const arb_t y, long prec)
    void acb_approx_abs(arb_t r, const acb_t z, long prec)
    void arb_approx_hypot(arb_t z, const arb_t x, const arb_t y, long prec)

    void acb_mat_approx_add(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)
    void acb_mat_approx_sub(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec)
    void acb_mat_approx_norm(arb_t res, acb_mat_t x, long prec)
    void acb_mat_approx_scalar_mul(acb_mat_t res, const acb_mat_t A, const acb_t c, long prec)
    void acb_mat_approx_scalar_mul_arb(acb_mat_t res, const acb_mat_t A, const arb_t c, long prec)
    void acb_mat_approx_scalar_addmul(acb_mat_t res, acb_mat_t B, const acb_mat_t A, const acb_t c, long prec)
    void acb_mat_approx_dotc(acb_t res, acb_mat_t x, acb_mat_t y, long prec)

    void acb_approx_complex_sign(acb_t res, acb_t z, arb_t z_abs, long prec)
    void lartg(acb_t c, acb_t s, acb_t r, acb_t f, acb_t g, long prec)

def test_gmres():
    import numpy as np
    from iterative_solvers_dp import gmres_mgs_dp
    cdef int dimen = 100
    cdef int prec = 53
    RBF = RealBallField(prec)
    CBF = ComplexBallField(prec)
    tol = RBF(1e-5)

    A = np.random.rand(dimen,dimen)+np.random.rand(dimen,dimen)*1j
    for i in range(dimen):
        A[i,i] *= 10
    b = np.random.rand(dimen,1)+np.random.rand(dimen,1)*1j
    cdef ComplexBall cb_cast
    cdef Acb_Mat A_arb = Acb_Mat(dimen,dimen)
    for i in range(dimen):
        for j in range(dimen):
            cb_cast = CBF(A[i,j])
            acb_set(acb_mat_entry(A_arb.value,i,j), cb_cast.value)
    cdef Acb_Mat b_arb = Acb_Mat(dimen,1)
    for i in range(dimen):
        cb_cast = CBF(b[i][0])
        acb_set(acb_mat_entry(b_arb.value,i,0), cb_cast.value)
    cdef Acb_Mat x0_arb = Acb_Mat(dimen,1) #Zero vector

    x_gmres_arb_wrap = gmres_mgs_arb_wrap(A_arb, b_arb, x0_arb, prec, tol)
    x_gmres_dp = gmres_mgs_dp(A, b)
    x_gmres_arb_wrap[0].str(10)
    print(x_gmres_dp[0])

def test_givens():
    CBF = ComplexBallField(53)
    Q = []
    Q_entry = Acb_Mat(2,2)
    cdef ComplexBall cb
    cb = CBF(0.8910632 +0.j)
    acb_set(acb_mat_entry(Q_entry.value,0,0),cb.value)
    cb = CBF(0.23526728+0.38814389j)
    acb_set(acb_mat_entry(Q_entry.value,0,1),cb.value)
    cb = CBF(-0.23526728+0.38814389j)
    acb_set(acb_mat_entry(Q_entry.value,1,0),cb.value)
    cb = CBF(0.8910632 +0.j)
    acb_set(acb_mat_entry(Q_entry.value,1,1),cb.value)
    Q.append(Q_entry)
    v_source = Acb_Mat(6,1)
    v = v_source.get_window(0,0,6,6)
    cb = CBF(0.22712716-2.39280626j)
    acb_set(acb_mat_entry(v.value,0,0),cb.value)
    cb = CBF(4.68967002+2.80871354j)
    acb_set(acb_mat_entry(v.value,1,0),cb.value)
    cb = CBF(2.3997535 +0.j)
    acb_set(acb_mat_entry(v.value,2,0),cb.value)
    cb = CBF(0.        +0.j)
    acb_set(acb_mat_entry(v.value,3,0),cb.value)
    cb = CBF(0.        +0.j)
    acb_set(acb_mat_entry(v.value,4,0),cb.value)
    cb = CBF(0.        +0.j)
    acb_set(acb_mat_entry(v.value,5,0),cb.value)
    apply_givens(Q,v,1,53)
    for i in range(6):
        acb_printd(acb_mat_entry(v.value,i,0),10)
        print('')
    #Should return:
    # (0.2155255644 + 0.3489235562j)  +/-  (0, 0j)
    # (5.054109916 + 3.153848315j)  +/-  (0, 0j)
    # (2.3997535 + 0j)  +/-  (0, 0j)
    # (0 + 0j)  +/-  (0, 0j)
    # (0 + 0j)  +/-  (0, 0j)
    # (0 + 0j)  +/-  (0, 0j)


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
        acb_approx_set(acb_mat_entry(v.value,j,0),v1)
        acb_approx_set(acb_mat_entry(v.value,j+1,0),v2)
    acb_clear(v1)
    acb_clear(v2)
    acb_clear(t)

def gmres_mgs_arb_wrap(Acb_Mat A, Acb_Mat b, Acb_Mat x0, int prec, RealBall tol, restrt=None, maxiter=None, M=None):
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
    M : array, matrix, sparse matrix, LinearOperator
        n x n, inverted preconditioner, i.e. solve M A x = M b.
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
        - The LinearOperator class is in scipy.sparse.linalg.interface.
          Use this class if you prefer to define A or M as a mat-vec routine
          as opposed to explicitly constructing the matrix.  A.psolve(..) is
          still supported as a legacy.
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
    # Convert inputs to linear system, with error checking
    #A, M, x, b, postprocess = make_system(A, M, x0, b)
    cdef Acb_Mat x = Acb_Mat(dimen,1)
    acb_mat_set(x.value, x0.value)

    cdef arb_t normr, normv, arb_tmp
    arb_init(normr)
    arb_init(normv)
    arb_init(arb_tmp)
    cdef acb_t c, s, acb_tmp, acb_tmp2
    acb_init(c)
    acb_init(s)
    acb_init(acb_tmp)
    acb_init(acb_tmp2)

    cdef acb_mat_t aliasing_avoider_matrix


#     def axpy(x,y,n,a): #returns scalar a * vector x + vector y
#         return a*x+y
#     def dotu(x,y): #returns dot product
#         return np.dot(x,y)
#     def dotc(x,y): #returns dot product of conjugate(x) and y
#         return np.dot(np.conjugate(x),y)
#     def scal(a,x): #returns scalar a * vector x
#         return a*x

#     # Make full use of direct access to BLAS by defining own norm
#     def norm(z):
#         return np.sqrt(np.real(dotc(z, z)))

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

    #normr = norm(r)
    sig_on()
    acb_mat_approx_norm(normr,r.value,prec)
    sig_off()

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
        Q = []  # Givens Rotations
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

            # # New Search Direction
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
            # H[inner+1, inner] = normv
            sig_on()
            acb_mat_approx_norm(normv,v.value,prec)
            sig_off()
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
                    acb_mat_win_cast = g.get_window(inner, 0, inner+2, 1) #Might want to replace this with a custom version
                    acb_mat_init(aliasing_avoider_matrix, 2, 1) #Arb currently does not work for aliasing window-matrices
                    acb_mat_set(aliasing_avoider_matrix, acb_mat_win_cast.value)
                    sig_on()
                    acb_mat_approx_mul(acb_mat_win_cast.value, Qblock.value, aliasing_avoider_matrix, prec)
                    sig_off()
                    acb_mat_clear(aliasing_avoider_matrix)

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


    #     # Apply preconditioner
    #     r = np.ravel(M*r)
        # normr = norm(r)
        sig_on()
        acb_mat_approx_norm(normr,r.value,prec)
        sig_off()

        # test for convergence
        if arb_lt(normr, tol.value) == 1:
            return (x, 0)

    # end outer loop

    return (x, niter)
