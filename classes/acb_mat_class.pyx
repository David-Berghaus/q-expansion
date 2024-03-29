# ****************************************************************************
#       Copyright (C) 2022 David Berghaus <berghaus@th.physik.uni-bonn.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_on, sig_str, sig_off
import numpy as np
cimport numpy as np

cdef extern from "complex.h":
    cdef double cimag(double complex)
    cdef double creal(double complex)
    cdef double cabs(double complex)
    cdef double complex cexp(double complex)

from sage.libs.arb.arb cimport *
from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_mat cimport *
from sage.matrix.matrix_complex_ball_dense cimport *
from sage.rings.complex_arb import ComplexBallField
from sage.matrix.matrix_space import MatrixSpace

from arblib_helpers.acb_approx cimport *

cdef class Acb_Mat():
    """
    This is a simple wrapper of acb_mat_t matrices into a python object which allows for more dynamical usage
    (for example to create instances inside a loop) as well as automatic deallocation.

    Why didn't we just use 'Matrix_complex_ball_dense' which is the builtin acb_mat_t wrapper of sage?
    'Matrix_complex_ball_dense' is constructed through a (python-declared) MatrixSpace. This means that the target precision is
    expected to be known during initialization. This does not reflect how acb_mat_t works internally and is inconvenient in our
    context because we often increase the precision of a matrix after initialization.
    Constructing a 'Matrix_complex_ball_dense' without specifying the precision during initialization does not seem to be possible.
    Consider the following example:
    cdef class Acb_Mat(Matrix_complex_ball_dense):
        def __cinit__(self,int nrows,int ncols):
            CBF = ComplexBallField(53) #53 bits serves as a placeholder here and does not restrict our arb precision
            parent = MatrixSpace(CBF,nrows,ncols)
            super().__init__(parent)
    This does not work with 'cdef-classes'...
    Therefore, during each initialization of a matrix we would need to construct a Matrix-Space with 'fake precision'
    (i.e., a precision that does not correspond to the actual internal precision) which is tedious and confusing.
    """
    def __cinit__(self, int nrows, int ncols):
        sig_str("Arb exception")
        acb_mat_init(self.value, nrows, ncols)
        sig_off()

    def __dealloc__(self):
        acb_mat_clear(self.value)

    cpdef Acb_Mat_Win get_window(self, int r1, int c1, int r2, int c2):
        return Acb_Mat_Win(self, r1, c1, r2, c2)

    cpdef Acb_Mat get_inv(self, int prec):
        nrows = acb_mat_nrows(self.value)
        ncols = acb_mat_ncols(self.value)
        cdef Acb_Mat V_inv = Acb_Mat(nrows, ncols)
        sig_on()
        acb_mat_approx_inv(V_inv.value, self.value, prec)
        sig_off()
        return V_inv

    def str(self, int digits): #Prints entries of matrix with specified precision
        nrows = acb_mat_nrows(self.value)
        ncols = acb_mat_ncols(self.value)
        cdef int i, j
        for i in range(nrows):
            print('['),
            for j in range(ncols):
                arf_printd(arb_midref(acb_realref(acb_mat_entry(self.value,i,j))), digits)
                print('+ '),
                arf_printd(arb_midref(acb_imagref(acb_mat_entry(self.value,i,j))), digits)
                print('j'),
                if j != ncols-1:
                    print(', '),
            print(']')
    
    def _set_mcbd(self, Matrix_complex_ball_dense B):
        """
        Constructs elements of Acb_Mat through Matrix_complex_ball_dense B
        """
        acb_mat_set(self.value, B.value)

    def _get_mcbd(self, int prec):
        """
        Returns Matrix_complex_ball_dense with entries identical to Acb_Mat
        """
        cdef int nrows, ncols
        nrows = acb_mat_nrows(self.value)
        ncols = acb_mat_ncols(self.value)
        CC = ComplexBallField(prec)
        cdef Matrix_complex_ball_dense B = MatrixSpace(CC, nrows, ncols).zero()
        acb_mat_set(B.value, self.value)
        return B

    cpdef nrows(self):
        return acb_mat_nrows(self.value)

    cpdef ncols(self):
        return acb_mat_ncols(self.value)

    # cpdef get_np(self):
    #     """
    #     Converts and returns np-matrix of self
    #     """
    #     nrows, ncols = self.nrows(), self.ncols()
    #     res = np.zeros(shape=(nrows,ncols), dtype=np.complex_)
    #     cdef double complex [:, ::1] res_view = res
    #     cdef int i, j
    #     for i in range(nrows):
    #         for j in range(ncols):
    #             res_view[i, j] = acb_to_dc(acb_mat_entry(self.value, i, j))
        
    #     return res

    # cpdef get_np_trunc(self, double tol):
    #     """
    #     Converts and returns np-matrix of self.
    #     All entries whose absolute value is less than tol are set to zero.
    #     """
    #     nrows, ncols = self.nrows(), self.ncols()
    #     res = np.zeros(shape=(nrows,ncols), dtype=np.complex_, order='F')
    #     cdef double complex [::1, :] res_view = res
    #     cdef int i, j
    #     trunc_col = ncols
    #     for j in range(ncols):
    #         above_tol = False #Set to True if value is larger than tol
    #         for i in range(nrows):
    #             res_view[i, j] = acb_to_dc(acb_mat_entry(self.value, i, j))
    #             if abs(creal(res_view[i,j])) > tol or abs(cimag(res_view[i,j])) > tol:
    #                 above_tol = True #Detected a large value so we are not finished yet
    #         if above_tol == False:
    #             trunc_col = j
    #             break
    #     return res[:, :trunc_col]

    cpdef set_np(self, A):
        """
        Sets self from np-matrix A
        """
        nrows, ncols = self.nrows(), self.ncols()
        cdef int i, j
        for i in range(nrows):
            for j in range(ncols):
                acb_set_d_d(acb_mat_entry(self.value,i,j), creal(A[i,j]), cimag(A[i,j]))

cdef class Acb_Mat_Win():
    """
    This is a simple wrapper of acb_mat_t windows into a python object which allows for more dynamical usage
    (for example to create instances inside a loop) as well as automatic deallocation.
    """
    def __cinit__(self, mat, int r1, int c1, int r2, int c2):
        """
        Initializes window to a window matrix into the submatrix of mat
        starting at the corner at row r1 and column c1 (inclusive) and ending at row r2 and column c2 (exclusive).
        """
        cdef Acb_Mat acb_mat_cast
        cdef Acb_Mat_Win acb_mat_win_cast
        if isinstance(mat, Acb_Mat):
            acb_mat_cast = mat
            sig_str("Arb exception")
            acb_mat_window_init(self.value, acb_mat_cast.value, r1, c1, r2, c2)
            sig_off()
            self.parents = [mat]
        elif isinstance(mat, Acb_Mat_Win):
            acb_mat_win_cast = mat
            sig_str("Arb exception")
            acb_mat_window_init(self.value, acb_mat_win_cast.value, r1, c1, r2, c2)
            sig_off()
            self.parents = mat.get_parents()[:] #Create copy
            (self.parents).append(mat)

    def __dealloc__(self):
        acb_mat_window_clear(self.value)

    def str(self, int digits): #Prints entries of matrix with specified precision
        nrows = acb_mat_nrows(self.value)
        ncols = acb_mat_ncols(self.value)
        cdef int i, j
        for i in range(nrows):
            print('['),
            for j in range(ncols):
                arf_printd(arb_midref(acb_realref(acb_mat_entry(self.value,i,j))), digits)
                print('+ '),
                arf_printd(arb_midref(acb_imagref(acb_mat_entry(self.value,i,j))), digits)
                print('j'),
                if j != ncols-1:
                    print(', '),
            print(']')
    
    def get_parents(self):
        return self.parents

    def _get_mcbd(self, int prec):
        """
        Returns Matrix_complex_ball_dense with entries identical to Acb_Mat_Win
        """
        cdef int nrows, ncols
        nrows = acb_mat_nrows(self.value)
        ncols = acb_mat_ncols(self.value)
        CC = ComplexBallField(prec)
        cdef Matrix_complex_ball_dense B = MatrixSpace(CC, nrows, ncols).zero()
        acb_mat_set(B.value, self.value)
        return B

    cpdef Acb_Mat_Win get_window(self, int r1, int c1, int r2, int c2):
        return Acb_Mat_Win(self, r1, c1, r2, c2)