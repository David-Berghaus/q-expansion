from sage.libs.arb.acb_mat cimport *

cdef class Acb_Mat():
    cdef acb_mat_t value

    cpdef Acb_Mat_Win get_window(self, int r1, int c1, int r2, int c2)

    cpdef Acb_Mat get_inv(self, int prec)

    cpdef nrows(self)

    cpdef ncols(self)

    # cpdef get_np(self)

    # cpdef get_np_trunc(self, double tol)

    cpdef set_np(self, A)

cdef class Acb_Mat_Win():
    cdef acb_mat_t value
    cdef list parents #We need to keep a reference to the parent objects to prevent them from getting prematurely garbage-collected

    cpdef Acb_Mat_Win get_window(self, int r1, int c1, int r2, int c2)