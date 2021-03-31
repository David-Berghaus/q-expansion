from classes.acb_mat_class cimport Acb_Mat

cdef class Block_Factored_Element():
    cdef Acb_Mat J
    cdef Acb_Mat W

cdef class Block_Factored_Matrix():
    cdef list A
    cdef list diag
    cdef int nc
    