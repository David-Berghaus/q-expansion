from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from point_matching.point_matching_arb_wrap import get_coefficients_arb_wrap, get_coefficients_gmres_arb_wrap
from iterative_solvers.gmres_arb_wrap import test_gmres

S = AutomorphicFormSpace(Gamma0(1),weight=12)