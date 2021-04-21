from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from point_matching.point_matching_arb_wrap import get_coefficients_gmres_arb_wrap, get_coefficients_ir_arb_wrap
from iterative_solvers.gmres_arb_wrap import test_gmres
from iterative_solvers.iterative_refinement_arb_wrap import test_ir

S = AutomorphicFormSpace(Gamma0(1),weight=12)