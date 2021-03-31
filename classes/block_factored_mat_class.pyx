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

    cdef act_on_vec(self, Acb_Mat b, Acb_Mat x, int prec):
        """
        Computes Block_Factored_Mat*x = b
        """
        if x.value == b.value:
            raise NameError("Aliasing not allowed here!")
        cdef int nc, cii, cjj, M
        nc = self.nc
        A = self.A
        diag = self.diag
        M = len(diag[0]) #diag[0] cannot be None if matrix is intialized
        cdef Acb_Mat tmp = Acb_Mat(M, 1)
        cdef Acb_Mat tmp2 = Acb_Mat(M, 1)
        cdef Acb_Mat J, W, b_cast, x_cast, diag_cast

        #Slice vectors x & b into blocks -> Maybe replace this with precomputed arguments later
        x_blocks = []
        b_blocks = []
        for cii in range(nc):
            x_blocks.append(x.get_window(cii*M,0,(cii+1)*M,1))
            b_blocks.append(b.get_window(cii*M,0,(cii+1)*M,1))

        for cii in range(nc):
            b_cast = b[cii]
            diag_cast = diag[cii]
            x_cast = x_blocks[cii]
            #b[cii] = -Diag[cii]*x[cii] (array multiplication)
            acb_mat_approx_array_mul(b_cast.value, diag_cast.value, x_cast.value, prec)
            acb_mat_neg(b_cast.value, b_cast.value)
            for cjj in range(nc):
                if A[cii][cjj] != None:
                    J = (A[cii][cjj]).J
                    W = (A[cii][cjj]).W
                    x_cast = x_blocks[cjj]
                    #tmp = W[cii][cjj]*x[cjj]
                    acb_mat_approx_mul(tmp.value, W.value, x_cast.value, prec)
                    #tmp2 = J[cii][cjj]*tmp
                    acb_mat_approx_mul(tmp2.value, J.value, tmp.value, prec)
                    #b[cii] += tmp2
                    acb_mat_approx_add(b_cast.value, b_cast.value, tmp2.value, prec)



