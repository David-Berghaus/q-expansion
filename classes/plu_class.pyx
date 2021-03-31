from cysignals.signals cimport sig_on, sig_str, sig_off
from cysignals.memory cimport sig_free,sig_malloc

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
    def __cinit__(self, Acb_Mat A, int prec):
        nrows = acb_mat_nrows(A.value)
        ncols = acb_mat_ncols(A.value)
        sig_str("Arb exception")
        acb_mat_init(self.value, nrows, ncols)
        sig_off()
        self.P = <long*>sig_malloc(nrows*sizeof(long))
        sig_on()
        lu_successful = acb_mat_approx_lu(self.P, self.value, A.value, prec)
        sig_off()
        if lu_successful == 0:
            raise ArithmeticError("LU-decomposition failed, please use higher precision!")

    def __dealloc__(self):
        acb_mat_clear(self.value)
        sig_free(self.P)


