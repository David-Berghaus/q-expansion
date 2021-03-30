#Based on: https://github.com/pyamg/pyamg/blob/main/pyamg/krylov/_gmres_mgs.py

#These functions are used as prototypes for the arb-version and should only be used for testing purposes

from __future__ import print_function
import numpy as np
import scipy as sp
import math
from scipy.sparse.linalg.isolve.utils import make_system
from scipy.sparse.sputils import upcast
from scipy.linalg import get_blas_funcs, get_lapack_funcs
from warnings import warn


__all__ = ['gmres_mgs']

def test_gmres():
    A = np.random.rand(10,10)+np.random.rand(10,10)*1j
    for i in range(10):
        A[i,i] *= 10
    b = np.random.rand(10,1)+np.random.rand(10,1)*1j
    res = gmres_mgs_dp(A,b)
    print("x: ", np.linalg.solve(A,b))
    return res

def apply_givens(Q, v, k):
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
    v1,v2 = 0.,0. #Should be complex!
    cdef int j
    for j in range(k):
        Qloc = Q[j]
        v1 = Qloc[0,0]*v[j]+Qloc[0,1]*v[j+1]
        v2 = Qloc[1,0]*v[j]+Qloc[1,1]*v[j+1]
        v[j] = v1
        v[j+1] = v2

def gmres_mgs_dp(A, b, x0=None, tol=1e-5, restrt=None, maxiter=None, M=None):
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
        initial guess, default is a vector of zeros
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
    cdef int inner = 0
    # Convert inputs to linear system, with error checking
    A, M, x, b, postprocess = make_system(A, M, x0, b)
    cdef int dimen = A.shape[0]

    def complex_sign(z,z_abs):
        if z==0:
            return 1
        else:
            return z/z_abs
    def lartg(f,g): #See Algorithm 1 of https://www.netlib.org/lapack/lawnspdf/lawn148.pdf
        if g == 0:
            c = 1
            s = 0
            r = f
        elif f == 0:
            g_abs = abs(g)
            c = 0
            s = complex_sign(np.conjugate(g),g_abs)
            r = g_abs
        else:
            f_abs, g_abs = abs(f), abs(g)
            fg_norm = math.hypot(f_abs,g_abs)
            f_sign = complex_sign(f,f_abs)
            c = f_abs/fg_norm
            s = f_sign*np.conjugate(g)/fg_norm
            r = f_sign*fg_norm
        return [c,s,r]
            
    def axpy(x,y,n,a): #returns scalar a * vector x + vector y
        return a*x+y
    def dotu(x,y): #returns dot product
        return np.dot(x,y)
    def dotc(x,y): #returns dot product of conjugate(x) and y
        return np.dot(np.conjugate(x),y)
    def scal(a,x): #returns scalar a * vector x
        return a*x

    # Make full use of direct access to BLAS by defining own norm
    def norm(z):
        return np.sqrt(np.real(dotc(z, z)))

    # Set number of outer and inner iterations
    if restrt:
        if maxiter:
            max_outer = maxiter
        else:
            max_outer = 1
        if restrt > dimen:
            warn('Setting number of inner iterations (restrt) to maximum\
                  allowed, which is A.shape[0] ')
            restrt = dimen
        max_inner = restrt
    else:
        max_outer = 1
        if maxiter > dimen:
            warn('Setting number of inner iterations (maxiter) to maximum\
                  allowed, which is A.shape[0] ')
            maxiter = dimen
        elif maxiter is None:
            maxiter = min(dimen, 40)
        max_inner = maxiter

    # Prep for method
    r = b - np.ravel(A*x)

    # Apply preconditioner
    r = np.ravel(M*r)
    normr = norm(r)

    # Use separate variable to track iterations.  If convergence fails, we
    # cannot simply report niter = (outer-1)*max_outer + inner.  Numerical
    # error could cause the inner loop to halt while the actual ||r|| > tol.
    niter = 0

    # Begin GMRES
    for outer in range(max_outer):

        # Preallocate for Givens Rotations, Hessenberg matrix and Krylov Space
        # Space required is O(dimen*max_inner).
        Q = []  # Givens Rotations
        # Upper Hessenberg matrix, which is then
        #   converted to upper tri with Givens Rots
        H = np.zeros((max_inner+1, max_inner+1), dtype=complex)
        V = np.zeros((dimen, max_inner+1), dtype=complex)  # Krylov Space
        # vs store the pointers to each column of V.
        #   This saves a considerable amount of time.
        vs = []
        # v = r/normr
        V[:, 0] = scal(1.0/normr, r)
        vs.append(V[:, 0])

        # This is the RHS vector for the problem in the Krylov Space
        g = np.zeros((dimen,), dtype=complex)
        g[0] = normr

        for inner in range(max_inner):

            # New Search Direction
            v = V[:, inner+1]
            v[:] = np.ravel(M*(A*vs[-1]))
            vs.append(v)
            normv_old = norm(v)

            #  Modified Gram Schmidt
            for k in range(inner+1):
                vk = vs[k]
                alpha = dotc(vk, v)
                H[k, inner] = alpha
                v[:] = axpy(vk, v, dimen, -alpha)

            normv = norm(v)
            H[inner+1, inner] = normv

            # Check for breakdown
            if H[inner+1, inner] != 0.0:
                v[:] = scal(1.0/H[inner+1, inner], v)

            # Apply previous Givens rotations to H
            if inner > 0:
                apply_givens(Q, H[:, inner], inner)

            # Calculate and apply next complex-valued Givens Rotation
            # ==> Note that if max_inner = dimen, then this is unnecessary
            # for the last inner
            #     iteration, when inner = dimen-1.
            if inner != dimen-1:
                if H[inner+1, inner] != 0:
                    [c, s, r] = lartg(H[inner, inner], H[inner+1, inner])
                    Qblock = np.array([[c, s], [-np.conjugate(s), c]],
                                      dtype=complex)
                    Q.append(Qblock)

                    # Apply Givens Rotation to g,
                    #   the RHS for the linear system in the Krylov Subspace.
                    g[inner:inner+2] = np.dot(Qblock, g[inner:inner+2])

                    # Apply effect of Givens Rotation to H
                    H[inner, inner] = dotu(Qblock[0, :],
                                           H[inner:inner+2, inner])
                    H[inner+1, inner] = 0.0

            niter += 1

            # Don't update normr if last inner iteration, because
            # normr is calculated directly after this loop ends.
            if inner < max_inner-1:
                normr = np.abs(g[inner+1])
                print(normr)
                if normr < tol:
                    break

        # end inner loop, back to outer loop

        # Find best update to x in Krylov Space V.  Solve inner x inner system.
        y = sp.linalg.solve(H[0:inner+1, 0:inner+1], g[0:inner+1]) #H is upper triangular
        update = np.ravel(V[:, :inner+1].dot(y.reshape(-1, 1)))
        x = x + update
        r = b - np.ravel(A*x)

        # Apply preconditioner
        r = np.ravel(M*r)
        normr = norm(r)

        # test for convergence
        if normr < tol:
            return (postprocess(x), 0)

    # end outer loop

    return (postprocess(x), niter)
