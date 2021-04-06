from cysignals.signals cimport sig_on, sig_off

from sage.libs.arb.acb_mat cimport *

from arblib_helpers.acb_approx cimport *
from classes.acb_mat_class cimport Acb_Mat, Acb_Mat_Win

cdef class Block_Factored_Element():
    def __cinit__(self, Acb_Mat J, Acb_Mat W):
        self.J = J
        self.W = W

cdef class Block_Factored_Mat():
    """
    This class stores a block-factored version of matrix V_tilde that can be used to compute
    the action of V_tilde on a vector without computing V_tilde explicitly.
    The class contains two variables:
    A : A two-dimensional dict of size ncusps which later gets filled with MxM matrices (J,W)
    diag : A a dict of size ncusps which later gets filled with Mx1 matrices that are diagonal elements
    """
    def __init__(self, int nc):
        A = [[None for cjj in range(nc)] for cii in range(nc)]
        diag = [None for cii in range(nc)]
        self.A = A
        self.diag = diag
        self.nc = nc
        self.max_len = 0

    cpdef Acb_Mat construct(self, int prec):
        """
        Explicitly performs the block-matrix multiplications at precision prec and returns the result
        """
        cdef int nc, cii, cjj, i, M
        nc = self.nc
        A = self.A
        diag = self.diag
        if diag[0] == None:
            raise NameError("Matrix is not properly initialized yet!")
        M = diag[0].nrows() #diag[0] cannot be None if matrix is initialized
        cdef Acb_Mat V_tilde = Acb_Mat(nc*M, nc*M)
        cdef Acb_Mat_Win V_view
        cdef Acb_Mat J, W, diag_cast
        cdef Block_Factored_Element block_factored_element
        
        for cii in range(nc):
            for cjj in range(nc):
                V_view = V_tilde.get_window(cii*M,cjj*M,(cii+1)*M,(cjj+1)*M)
                if A[cii][cjj] != None:
                    block_factored_element = A[cii][cjj]
                    J = block_factored_element.J
                    W = block_factored_element.W
                    sig_on()
                    acb_mat_approx_mul(V_view.value, J.value, W.value, prec)
                    sig_off()
                if cii == cjj:
                    diag_cast = diag[cii]
                    for i in range(M):
                        acb_approx_sub(acb_mat_entry(V_view.value,i,i),acb_mat_entry(V_view.value,i,i),acb_mat_entry(diag_cast.value,i,0),prec)   
        
        return V_tilde

    cpdef _get_max_len(self):
        """
        Returns the largest row-length of W (i.e. the highest amount of matching-points per cusp-pair)
        """
        cdef Block_Factored_Element block_factored_element
        if self.max_len == 0: #Not initialized yet so we need to compute it
            if self.diag[0] == None:
                raise NameError("Matrix is not properly initialized yet!")
            nc = self.nc
            A = self.A
            lengths = []

            for cii in range(nc):
                for cjj in range(nc):
                    if A[cii][cjj] != None:
                        block_factored_element = A[cii][cjj]
                        lengths.append((block_factored_element.W).nrows())

            self.max_len = max(lengths)
        
        return self.max_len


    cpdef act_on_vec(self, Acb_Mat b, Acb_Mat x, int prec):
        """
        Computes Block_Factored_Mat*x = b
        """
        cdef int nc, cii, cjj, M
        nc = self.nc
        A = self.A
        diag = self.diag
        if x.value == b.value:
            raise NameError("Aliasing not allowed here!")
        if diag[0] == None:
            raise NameError("Matrix is not properly initialized yet!")
        M = diag[0].nrows() #diag[0] cannot be None if matrix is initialized
        cdef Acb_Mat tmp = Acb_Mat(self._get_max_len(), 1)
        cdef Acb_Mat_Win tmp_view
        cdef Acb_Mat tmp2 = Acb_Mat(M, 1)
        cdef Acb_Mat J, W, diag_cast
        cdef Acb_Mat_Win b_cast, x_cast
        cdef Block_Factored_Element block_factored_element

        #Slice vectors x & b into blocks -> Maybe replace this with precomputed arguments later
        x_blocks = []
        b_blocks = []
        for cii in range(nc):
            x_blocks.append(x.get_window(cii*M,0,(cii+1)*M,1))
            b_blocks.append(b.get_window(cii*M,0,(cii+1)*M,1))

        for cii in range(nc):
            b_cast = b_blocks[cii]
            diag_cast = diag[cii]
            x_cast = x_blocks[cii]
            #b[cii] = -Diag[cii]*x[cii]
            sig_on()
            acb_mat_approx_left_mul_diag(b_cast.value, diag_cast.value, x_cast.value, prec)
            sig_off()
            sig_on()
            acb_mat_neg(b_cast.value, b_cast.value)
            sig_off()
            for cjj in range(nc):
                if A[cii][cjj] != None:
                    block_factored_element = A[cii][cjj]
                    J = block_factored_element.J
                    W = block_factored_element.W
                    x_cast = x_blocks[cjj]
                    tmp_view = tmp.get_window(0, 0, W.nrows(), 1)
                    #tmp = W[cii][cjj]*x[cjj]
                    sig_on()
                    acb_mat_approx_mul(tmp_view.value, W.value, x_cast.value, prec)
                    sig_off()
                    #tmp2 = J[cii][cjj]*tmp
                    sig_on()
                    acb_mat_approx_mul(tmp2.value, J.value, tmp_view.value, prec)
                    sig_off()
                    #b[cii] += tmp2
                    sig_on()
                    acb_mat_approx_add(b_cast.value, b_cast.value, tmp2.value, prec)
                    sig_off()

    cpdef act_on_vec_win(self, Acb_Mat_Win b, Acb_Mat_Win x, int prec):
        """
        Computes Block_Factored_Mat*x = b
        """
        cdef int nc, cii, cjj, M
        nc = self.nc
        A = self.A
        diag = self.diag
        if x.value == b.value:
            raise NameError("Aliasing not allowed here!")
        if diag[0] == None:
            raise NameError("Matrix is not properly initialized yet!")
        M = diag[0].nrows() #diag[0] cannot be None if matrix is initialized
        cdef Acb_Mat tmp = Acb_Mat(self._get_max_len(), 1)
        cdef Acb_Mat_Win tmp_view
        cdef Acb_Mat tmp2 = Acb_Mat(M, 1)
        cdef Acb_Mat J, W, diag_cast
        cdef Acb_Mat_Win b_cast, x_cast
        cdef Block_Factored_Element block_factored_element

        #Slice vectors x & b into blocks -> Maybe replace this with precomputed arguments later
        x_blocks = []
        b_blocks = []
        for cii in range(nc):
            x_blocks.append(x.get_window(cii*M,0,(cii+1)*M,1))
            b_blocks.append(b.get_window(cii*M,0,(cii+1)*M,1))

        for cii in range(nc):
            b_cast = b_blocks[cii]
            diag_cast = diag[cii]
            x_cast = x_blocks[cii]
            #b[cii] = -Diag[cii]*x[cii]
            sig_on()
            acb_mat_approx_left_mul_diag(b_cast.value, diag_cast.value, x_cast.value, prec)
            sig_off()
            sig_on()
            acb_mat_neg(b_cast.value, b_cast.value)
            sig_off()
            for cjj in range(nc):
                if A[cii][cjj] != None:
                    block_factored_element = A[cii][cjj]
                    J = block_factored_element.J
                    W = block_factored_element.W
                    x_cast = x_blocks[cjj]
                    tmp_view = tmp.get_window(0, 0, W.nrows(), 1)
                    #tmp = W[cii][cjj]*x[cjj]
                    sig_on()
                    acb_mat_approx_mul(tmp_view.value, W.value, x_cast.value, prec)
                    sig_off()
                    #tmp2 = J[cii][cjj]*tmp
                    sig_on()
                    acb_mat_approx_mul(tmp2.value, J.value, tmp_view.value, prec)
                    sig_off()
                    #b[cii] += tmp2
                    sig_on()
                    acb_mat_approx_add(b_cast.value, b_cast.value, tmp2.value, prec)
                    sig_off()

    cpdef nrows(self):
        """
        Returns full dimension of block-factored matrix
        """
        nc = self.nc
        diag = self.diag
        if diag[0] == None:
            raise NameError("Matrix is not properly initialized yet!")
        M = diag[0].nrows() #diag[0] cannot be None if matrix is initialized
        return M*nc



