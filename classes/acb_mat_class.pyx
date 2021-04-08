from cysignals.signals cimport sig_on, sig_str, sig_off

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