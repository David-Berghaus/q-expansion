cimport numpy as np

from sage.libs.arb.acb_mat cimport *

from classes.acb_mat_class cimport Acb_Mat, Acb_Mat_Win
from classes.acb_dft_class cimport Acb_DFT

cdef class J_class():
    cdef Acb_Mat J
    cdef bint use_FFT
    cdef bint is_initialized

    #These objects are only needed if use_FFT=True
    cdef Acb_Mat D_L
    cdef Acb_Mat D_R
    cdef object DFT_precomp
    cdef np.int64_t[::1] j_values
    cdef Acb_Mat two_Q_vec
    cdef Acb_Mat coord_len_vec

    cdef act_on_vec(self, acb_mat_t b, acb_mat_t x, int prec)

cdef class W_class():
    cdef Acb_Mat W
    cdef bint use_Horner
    cdef bint is_initialized
    cdef int _nrows

    cdef act_on_vec(self, acb_mat_t b, acb_mat_t x, int prec)

cdef class Block_Factored_Element():
    cdef J_class J
    cdef W_class W

cdef class Block_Factored_Mat():
    cdef list A
    cdef list diag
    cdef list diag_inv
    cdef int nc
    cdef int max_len
    cdef object S
    cdef object M
    cdef object Q
    cdef object Y
    cdef object normalization
    cdef object pb
    cdef object is_cuspform
    cdef bint parameters_for_dp_initialized

    cpdef _get_max_len(self)

    cpdef act_on_vec_sc(self, Acb_Mat b, Acb_Mat x, int prec)

    cpdef construct_sc_np(self)

    cpdef act_on_vec(self, Acb_Mat b, Acb_Mat x, int prec, is_scaled)

    cpdef act_on_vec_win_sc(self, Acb_Mat_Win b, Acb_Mat_Win x, int prec)

    cpdef act_on_vec_win(self, Acb_Mat_Win b, Acb_Mat_Win x, int prec, is_scaled)

    cpdef nrows(self)

    cpdef diag_scale_vec(self, Acb_Mat res, Acb_Mat x, int prec)

    cpdef diag_inv_scale_vec(self, Acb_Mat res, Acb_Mat x, int prec)
    