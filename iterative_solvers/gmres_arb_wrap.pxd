from sage.rings.real_arb cimport RealBall

from classes.acb_mat_class cimport Acb_Mat

cpdef gmres_mgs_arb_wrap(Acb_Mat A, Acb_Mat b, Acb_Mat x0, int prec, RealBall tol, restrt=*, maxiter=*, M=*, PLU=*)