from sage.rings.real_arb cimport RealBall

from classes.acb_mat_class cimport Acb_Mat

cpdef gmres_mgs_arb_wrap(A, Acb_Mat b, int prec, RealBall tol, x0=*, restrt=*, maxiter=*, M=*, PLU=*, is_scaled=*)

cdef mat_vec_mul(Acb_Mat b, A, Acb_Mat x, int prec, is_scaled)