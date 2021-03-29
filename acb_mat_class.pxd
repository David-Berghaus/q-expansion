from sage.libs.arb.acb_mat cimport *

cdef class Acb_Mat():
    cdef acb_mat_t value

    cpdef Acb_Mat_Win get_window(self, int r1, int c1, int r2, int c2)

cdef class Acb_Mat_Win():
    cdef acb_mat_t value