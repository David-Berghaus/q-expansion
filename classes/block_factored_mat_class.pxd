from classes.acb_mat_class cimport Acb_Mat, Acb_Mat_Win

cdef class Block_Factored_Element():
    cdef Acb_Mat J
    cdef Acb_Mat W

cdef class Block_Factored_Mat():
    cdef list A
    cdef list diag
    cdef list diag_inv
    cdef int nc
    cdef int max_len

    cpdef Acb_Mat construct_non_sc(self, int prec)

    cpdef Acb_Mat construct_sc(self, int prec)

    cpdef Acb_Mat construct(self, int prec, is_scaled)

    cpdef _get_max_len(self)

    cpdef act_on_vec_non_sc(self, Acb_Mat b, Acb_Mat x, int prec)

    cpdef act_on_vec_sc(self, Acb_Mat b, Acb_Mat x, int prec)

    cpdef construct_sc_np(self)

    cpdef act_on_vec(self, Acb_Mat b, Acb_Mat x, int prec, is_scaled)

    cpdef act_on_vec_win_non_sc(self, Acb_Mat_Win b, Acb_Mat_Win x, int prec)

    cpdef act_on_vec_win_sc(self, Acb_Mat_Win b, Acb_Mat_Win x, int prec)

    cpdef act_on_vec_win(self, Acb_Mat_Win b, Acb_Mat_Win x, int prec, is_scaled)

    cpdef nrows(self)

    cpdef diag_scale_vec(self, Acb_Mat res, Acb_Mat x, int prec)

    cpdef diag_inv_scale_vec(self, Acb_Mat res, Acb_Mat x, int prec)
    