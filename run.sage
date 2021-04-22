import numpy as np

from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from point_matching.point_matching_arb_wrap import get_coefficients_gmres_arb_wrap, get_coefficients_ir_arb_wrap
from iterative_solvers.gmres_arb_wrap import test_gmres
from iterative_solvers.iterative_refinement_arb_wrap import test_ir
from point_matching.point_matching_dp import get_V_tilde_matrix_b_haupt_dp, get_coefficients_haupt_dp

S = AutomorphicFormSpace(Gamma0(1),weight=12)