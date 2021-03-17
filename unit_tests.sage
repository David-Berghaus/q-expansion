import numpy as np

from point_matching_dp import get_coefficients_dp, get_V_tilde_matrix_b_dp, get_V_tilde_matrix_dp
from point_matching_arb_wrap import get_coefficients_arb_wrap, get_V_tilde_matrix_b_arb_wrap, get_V_tilde_matrix_arb_wrap
from psage.modform.maass.automorphic_forms import AutomorphicFormSpace

def run_unit_tests():
    run_unit_tests_dp()
    run_unit_tests_arb_wrap()

def run_unit_tests_dp():
    test_get_coefficients_dp()
    test_get_V_tilde_matrix_b_dp()
    test_get_V_tilde_matrix_dp()

def test_get_coefficients_dp():
    S = AutomorphicFormSpace(Gamma0(2),weight=12) #Search for a multiplicity two old-form
    C = get_coefficients_dp(S)
    assert abs(C[0][0]-252)/abs(C[0][0]+252) < 1e-12
    assert abs(C[1][0]-(-2048))/abs(C[1][0]+(-2048)) < 1e-12
    print("test_get_coefficients_dp ok")

def test_get_V_tilde_matrix_b_dp():
    S = AutomorphicFormSpace(Gamma0(8),weight=4)
    V,b = get_V_tilde_matrix_b_dp(S,66,0.09742785792574934)
    C = np.linalg.solve(V,b)
    assert abs(C[0][0]) < 1e-12
    assert abs(C[1][0]-(-4))/abs(C[1][0]+(-4)) < 1e-12
    print("test_get_V_tilde_matrix_b_dp ok")

def test_get_V_tilde_matrix_dp():
    S = AutomorphicFormSpace(Gamma0(11),weight=2)
    T = get_V_tilde_matrix_dp(S,82,0.07085662394599952)
    V = T[1:,1:]
    b = -T[1:,0]
    C = np.linalg.solve(V,b)
    assert abs(C[0]-(-2))/abs(C[0]+(-2)) < 1e-12
    assert abs(C[1]-(-1))/abs(C[1]+(-1)) < 1e-12
    print("test_get_V_tilde_matrix_dp ok")

def run_unit_tests_arb_wrap():
    test_get_coefficients_arb_wrap()
    test_get_V_tilde_matrix_b_arb_wrap()
    test_get_V_tilde_matrix_arb_wrap()

def test_get_coefficients_arb_wrap():
    S = AutomorphicFormSpace(Gamma0(2),weight=12) #Search for a multiplicity two old-form
    C = get_coefficients_arb_wrap(S,16)
    assert abs(C[0][0]-252)/abs(C[0][0]+252) < RBF(1e-10)
    assert abs(C[1][0]-(-2048))/abs(C[1][0]+(-2048)) < RBF(1e-10)
    print("test_get_coefficients_arb_wrap ok")

def test_get_V_tilde_matrix_b_arb_wrap():
    S = AutomorphicFormSpace(Gamma0(8),4)
    RBF = RealBallField(128)
    V,b = get_V_tilde_matrix_b_arb_wrap(S,66,RBF(0.09742785792574934),128)
    C = V\b
    assert abs(C[0][0]) < RBF(1e-12)
    assert abs(C[1][0]-(-4))/abs(C[1][0]+(-4)) < RBF(1e-12)
    print("test_get_V_tilde_matrix_b_arb_wrap ok")

def test_get_V_tilde_matrix_arb_wrap():
    S = AutomorphicFormSpace(Gamma0(11),2)
    RBF = RealBallField(128)
    T = get_V_tilde_matrix_arb_wrap(S,82,RBF(0.07085662394599952),128)
    V = T[1:,1:]
    b = -T[1:,0]
    C = V\b
    assert abs(C[0][0]-(-2))/abs(C[0][0]+(-2)) < RBF(1e-12)
    assert abs(C[1][0]-(-1))/abs(C[1][0]+(-1)) < RBF(1e-12)
    print("test_get_V_tilde_matrix_arb_wrap ok")