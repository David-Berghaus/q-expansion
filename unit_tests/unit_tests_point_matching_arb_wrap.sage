import numpy as np

from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from psage.groups.permutation_alg import MyPermutation
from psage.modform.arithgroup.mysubgroup import MySubgroup

from point_matching.point_matching_arb_wrap import get_coefficients_arb_wrap, get_V_tilde_matrix_b_arb_wrap, get_V_tilde_matrix_arb_wrap, get_coefficients_haupt_ir_arb_wrap, get_coefficients_eisenstein_ir_arb_wrap

def run_unit_tests_point_matching_arb_wrap():
    test_get_coefficients_arb_wrap()
    test_get_V_tilde_matrix_b_arb_wrap()
    test_get_V_tilde_matrix_arb_wrap()
    test_get_coefficients_haupt_ir_arb_wrap()
    test_get_coefficients_eisenstein_ir_arb_wrap()

def test_get_coefficients_arb_wrap():
    S = AutomorphicFormSpace(Gamma0(2),weight=12) #Search for a multiplicity two old-form
    C = get_coefficients_arb_wrap(S,16)._get_mcbd(53)
    assert abs(C[0][0]-252)/abs(C[0][0]+252) < RBF(1e-10)
    assert abs(C[1][0]-(-2048))/abs(C[1][0]+(-2048)) < RBF(1e-10)
    print("test_get_coefficients_arb_wrap ok")

def test_get_V_tilde_matrix_b_arb_wrap():
    S = AutomorphicFormSpace(Gamma0(8),4) #This example contains some empty coordinate lists
    prec = 64
    RBF = RealBallField(prec)
    tmp = get_V_tilde_matrix_b_arb_wrap(S,66,RBF(0.09742785792574934),prec)
    V,b = tmp[0]._get_mcbd(prec), tmp[1]._get_mcbd(prec)
    C = V\b
    assert abs(C[0][0]) < RBF(1e-12)
    assert abs(C[1][0]-(-4))/abs(C[1][0]+(-4)) < RBF(1e-12)
    print("test_get_V_tilde_matrix_b_arb_wrap ok")

def test_get_V_tilde_matrix_arb_wrap():
    S = AutomorphicFormSpace(Gamma0(11),2)
    RBF = RealBallField(128)
    T = get_V_tilde_matrix_arb_wrap(S,82,RBF(0.07085662394599952),128)._get_mcbd(128)
    V = T[1:,1:]
    b = -T[1:,0]
    C = V\b
    assert abs(C[0][0]-(-2))/abs(C[0][0]+(-2)) < RBF(1e-12)
    assert abs(C[1][0]-(-1))/abs(C[1][0]+(-1)) < RBF(1e-12)
    print("test_get_V_tilde_matrix_arb_wrap ok")

def test_get_coefficients_haupt_ir_arb_wrap():
    S = AutomorphicFormSpace(Gamma0(4),weight=0)
    C = get_coefficients_haupt_ir_arb_wrap(S,16)._get_mcbd(53)
    assert abs(C[0][0]-20)/abs(C[0][0]+20) < RBF(1e-10)
    print("test_get_coefficients_haupt_ir_arb_wrap ok")

def test_get_coefficients_eisenstein_ir_arb_wrap():
    S = AutomorphicFormSpace(Gamma0(5),weight=2)
    C = get_coefficients_eisenstein_ir_arb_wrap(S,16)._get_mcbd(53)
    assert abs(C[0][0]-6)/abs(C[0][0]+6) < RBF(1e-10)
    print("test_get_coefficients_eisenstein_ir_arb_wrap ok")