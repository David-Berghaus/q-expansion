from sage.libs.arb.acb_mat cimport *

cdef class Acb_Mat():
    cdef acb_mat_t value

    cpdef Acb_Mat_Win get_window(self, int r1, int c1, int r2, int c2)

    cpdef nrows(self)

    cpdef ncols(self)

cdef class Acb_Mat_Win():
    cdef acb_mat_t value
    cdef Acb_Mat parent #We need to keep a reference to the parent object to prevent it from getting prematurely garbage-collected