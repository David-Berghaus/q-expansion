import numpy as np

from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from psage.groups.permutation_alg import MyPermutation
from psage.modform.arithgroup.mysubgroup import MySubgroup

from point_matching.point_matching_arb_wrap import get_coefficients_gmres_arb_wrap, get_coefficients_ir_arb_wrap, get_coefficients_haupt_ir_arb_wrap, get_coefficients_eisenstein_ir_arb_wrap
from iterative_solvers.gmres_arb_wrap import test_gmres
from iterative_solvers.iterative_refinement_arb_wrap import test_ir
from point_matching.point_matching_dp import get_V_tilde_matrix_b_haupt_dp
from classes.gamma_2_subgroup import Gamma_2_Subgroup
from classes.modform_class import ModForm
from point_matching.fft_test import test_fft

S = AutomorphicFormSpace(Gamma0(1),weight=12)

#These notations are part of the Fiori, Franc paper. Note that we do not get the same permT...
G1 = MySubgroup(o2='(1 2)(3 4)(5 6)(7)',o3='(1)(2 3 5)(4 6 7)')
H1 = MySubgroup(o2='(1 2)(3 4)(5 6)(7)',o3='(1)(2 3 5)(7 6 4)')
U1 = MySubgroup(o2='(1 2)(3 4)(5)(6 7)',o3='(1)(2 3 5)(4 6 7)')