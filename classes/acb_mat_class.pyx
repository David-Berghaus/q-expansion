from cysignals.signals cimport sig_on, sig_str, sig_off

from sage.libs.arb.arb cimport *
from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_mat cimport *
from sage.matrix.matrix_complex_ball_dense cimport *

from arblib_helpers.acb_mat_approx cimport *

cdef class Acb_Mat():
    """
    This is a simple wrapper of acb_mat_t matrices into a python object which allows for more dynamical usage
    (for example to create instances inside a loop) as well as automatic deallocation.
    """
    def __cinit__(self, int nrows, int ncols):
        sig_str("Arb exception")
        acb_mat_init(self.value, nrows, ncols)
        sig_off()

    def __dealloc__(self):
        acb_mat_clear(self.value)

    cpdef Acb_Mat_Win get_window(self, int r1, int c1, int r2, int c2):
        return Acb_Mat_Win(self, r1, c1, r2, c2)

    def str(self, int digits): #Prints entries of matrix with specified precision
        nrows = acb_mat_nrows(self.value)
        ncols = acb_mat_ncols(self.value)
        cdef int i, j
        for i in range(nrows):
            print('['),
            for j in range(ncols):
                #acb_printd(acb_mat_entry(self.value,i,j), digits)
                arf_printd(arb_midref(acb_realref(acb_mat_entry(self.value,i,j))), digits)
                print('+ '),
                arf_printd(arb_midref(acb_imagref(acb_mat_entry(self.value,i,j))), digits)
                print('j'),
                if j != ncols-1:
                    print(', '),
            print(']')
    
    def _set_mcbd(self, Matrix_complex_ball_dense B): #Sets elements of Acb_Mat to the ones specified by B
        cdef int i, j, nrows, ncols
        nrows = acb_mat_nrows(B.value)
        ncols = acb_mat_ncols(B.value)
        for i in range(nrows):
            for j in range(ncols):
                acb_set(acb_mat_entry(self.value,i,j), acb_mat_entry(B.value,i,j))


cdef class Acb_Mat_Win():
    """
    This is a simple wrapper of acb_mat_t windows into a python object which allows for more dynamical usage
    (for example to create instances inside a loop) as well as automatic deallocation.
    """
    def __cinit__(self, Acb_Mat mat, int r1, int c1, int r2, int c2):
        """
        Initializes window to a window matrix into the submatrix of mat
        starting at the corner at row r1 and column c1 (inclusive) and ending at row r2 and column c2 (exclusive).
        """
        sig_str("Arb exception")
        acb_mat_window_init(self.value, mat.value, r1, c1, r2, c2)
        sig_off()

    def __dealloc__(self):
        acb_mat_window_clear(self.value)

    def str(self, int digits): #Prints entries of matrix with specified precision
        nrows = acb_mat_nrows(self.value)
        ncols = acb_mat_ncols(self.value)
        cdef int i, j
        for i in range(nrows):
            print('['),
            for j in range(ncols):
                #acb_printd(acb_mat_entry(self.value,i,j), digits)
                arf_printd(arb_midref(acb_realref(acb_mat_entry(self.value,i,j))), digits)
                print('+ '),
                arf_printd(arb_midref(acb_imagref(acb_mat_entry(self.value,i,j))), digits)
                print('j'),
                if j != ncols-1:
                    print(', '),
            print(']')