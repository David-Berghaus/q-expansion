from sage.libs.arb.acb_mat cimport *

cdef class PLU_Mat():
    cdef acb_mat_t value
    cdef long *P