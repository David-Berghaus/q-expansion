from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from point_matching.point_matching_dp import get_coefficients_dp
from iterative_solvers.iterative_solvers_arb_wrap import test_gmres
S = AutomorphicFormSpace(Gamma0(1),weight=12)