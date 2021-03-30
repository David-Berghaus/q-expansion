from classes.acb_mat_class cimport Acb_Mat, Acb_Mat_Win

cdef class Block_Factored_element():
    cdef Acb_Mat J
    cdef Acb_Mat W

cdef class Block_Factored_Mat():
    """
    This class stores a block-factored version of matrix V_tilde that can be used to compute
    the action of V_tilde on a vector without computing V_tilde explicitly.
    The class contains two variables:
    A : A two-dimensional list of size ncusps consisting of MxM matrices (J,W)
    diag : A a list of size ncusps consisting of Mx1 matrices that are diagonal elements
    """