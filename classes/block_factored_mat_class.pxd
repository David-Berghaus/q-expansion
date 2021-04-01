from classes.acb_mat_class cimport Acb_Mat

cdef class Block_Factored_Element():
    cdef Acb_Mat J
    cdef Acb_Mat W

cdef class Block_Factored_Mat():
    cdef list A
    cdef list diag
    cdef int nc
    cdef int max_len

    cpdef Acb_Mat construct(self, int prec)

    cpdef _get_max_len(self)

    cpdef act_on_vec(self, Acb_Mat b, Acb_Mat x, int prec)
    