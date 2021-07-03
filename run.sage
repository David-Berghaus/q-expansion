import numpy as np

from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from psage.groups.permutation_alg import MyPermutation
from psage.modform.arithgroup.mysubgroup import MySubgroup

from point_matching.point_matching_arb_wrap import get_coefficients_cuspform_ir_arb_wrap, get_coefficients_haupt_ir_arb_wrap, get_coefficients_modform_ir_arb_wrap, digits_to_bits, get_horo_height_arb_wrap
from classes.gamma_2_subgroup import Gamma_2_Subgroup
from classes.modform_class import ModForm
from point_matching.fft_test import test_fft

S = AutomorphicFormSpace(Gamma0(1),weight=12)

from classes.approx_modform import ApproxModForm, get_approxmodform_basis