from arblib_helpers.acb_approx cimport *

cdef class Acb_DFT():
    cdef acb_dft_pre_t value
    cdef int guard_bits #The DFT uses interval-arithmetic so we can expect to loose some precision (although the precision-loss is quite low)