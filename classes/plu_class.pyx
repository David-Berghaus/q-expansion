from cysignals.signals cimport sig_on, sig_str, sig_off
from cysignals.memory cimport sig_free,sig_malloc
import numpy as np
cimport numpy as np
from scipy.sparse import csc_matrix, linalg as sla

cdef extern from "complex.h":
    cdef double cimag(double complex)
    cdef double creal(double complex)
    cdef double cabs(double complex)
    cdef double complex cexp(double complex)

from sage.libs.arb.arb cimport *
from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_mat cimport *
from sage.matrix.matrix_complex_ball_dense cimport *

from arblib_helpers.acb_approx cimport *
from classes.acb_mat_class cimport Acb_Mat

cdef class PLU_Mat():
    """
    Computes and stores the LU decomposition of input matrix A at precision prec
    """
    def __cinit__(self, A, int prec):
        cdef int i, j, nrows, ncols
        cdef double tmp_double
        cdef Acb_Mat acb_mat_cast
        cdef double complex [::1, :] L, U #The result of scipy lu will be column major
        cdef int [:] piv

        if isinstance(A, Acb_Mat):
            acb_mat_cast = A
            nrows = acb_mat_nrows(acb_mat_cast.value)
            ncols = acb_mat_ncols(acb_mat_cast.value)
            sig_str("Arb exception")
            acb_mat_init(self.value, nrows, ncols)
            sig_off()
            self.P = <long*>sig_malloc(nrows*sizeof(long))
            sig_on()
            lu_successful = acb_mat_approx_lu(self.P, self.value, acb_mat_cast.value, prec)
            sig_off()
            if lu_successful == 0:
                raise ArithmeticError("LU-decomposition failed, please use higher precision!")
                
        elif isinstance(A, np.ndarray):
            nrows, ncols = A.shape
            A_sp = csc_matrix(A)
            lu_sp = sla.splu(A_sp, permc_spec="NATURAL", diag_pivot_thresh=0, options={"SymmetricMode":True}) #See: https://github.com/scipy/scipy/issues/7700#issuecomment-441131694

            sig_str("Arb exception")
            acb_mat_init(self.value, nrows, ncols)
            sig_off()
            self.P = <long*>sig_malloc(nrows*sizeof(long))
   
            L = lu_sp.L.A
            for i in range(1, nrows):
                for j in range(i):
                    tmp_double = creal(L[i,j])
                    if abs(tmp_double) > 1e-17: #If result is less than rounding error, we can set it to zero for faster multiplication later
                        arf_set_d(arb_midref(acb_realref(acb_mat_entry(self.value,i,j))), tmp_double)
                    tmp_double = cimag(L[i,j])
                    if abs(tmp_double) > 1e-17: #If result is less than rounding error, we can set it to zero for faster multiplication later
                        arf_set_d(arb_midref(acb_imagref(acb_mat_entry(self.value,i,j))), tmp_double)

            U = lu_sp.U.A
            for i in range(nrows):
                for j in range(i, ncols):
                    tmp_double = creal(U[i,j])
                    if abs(tmp_double) > 1e-17: #If result is less than rounding error, we can set it to zero for faster multiplication later
                        arf_set_d(arb_midref(acb_realref(acb_mat_entry(self.value,i,j))), tmp_double)
                    tmp_double = cimag(U[i,j])
                    if abs(tmp_double) > 1e-17: #If result is less than rounding error, we can set it to zero for faster multiplication later
                        arf_set_d(arb_midref(acb_imagref(acb_mat_entry(self.value,i,j))), tmp_double)

            P = self.P
            piv = lu_sp.perm_r
            for i in range(nrows):
                P[i] = piv[i]

    def __dealloc__(self):
        acb_mat_clear(self.value)
        sig_free(self.P)


