from sage.rings.real_arb cimport RealBall

from classes.acb_mat_class cimport Acb_Mat
from classes.plu_class cimport PLU_Mat
from classes.block_factored_mat_class cimport Block_Factored_Mat

cpdef iterative_refinement_arb_wrap(Block_Factored_Mat A, Acb_Mat b, int prec, RealBall tol, PLU_Mat PLU, x0=*, maxiter=*, mix_prec=*, is_scaled=*)