import numpy as np

from point_matching_dp import get_coefficients_dp, get_V_tilde_matrix_b_dp, get_V_tilde_matrix_dp
from psage.modform.maass.automorphic_forms import AutomorphicFormSpace

def run_unit_tests():
    test_get_coefficients_dp()
    test_get_V_tilde_matrix_b_dp()
    test_get_V_tilde_matrix_dp()

def test_get_coefficients_dp():
    S = AutomorphicFormSpace(Gamma0(2))
    C = get_coefficients_dp(S,12,2) #Search for a multiplicity two old-form
    assert abs(C[0][0]-252)/abs(C[0][0]+252) < 1e-12
    assert abs(C[1][0]-(-2048))/abs(C[1][0]+(-2048)) < 1e-12
    print("test_get_coefficients_dp ok")

def test_get_V_tilde_matrix_b_dp():
    S = AutomorphicFormSpace(Gamma0(8))
    V,b = get_V_tilde_matrix_b_dp(S,66,0.09742785792574934,4,1)
    C = np.linalg.solve(V,b)
    assert abs(C[0][0]) < 1e-12
    assert abs(C[1][0]-(-4))/abs(C[1][0]+(-4)) < 1e-12
    print("test_get_V_tilde_matrix_b_dp ok")

def test_get_V_tilde_matrix_dp():
    S = AutomorphicFormSpace(Gamma0(11))
    T = get_V_tilde_matrix_dp(S,82,0.07085662394599952,2)
    V = T[1:,1:]
    b = -T[1:,0]
    C = np.linalg.solve(V,b)
    assert abs(C[0]-(-2))/abs(C[0]+(-2)) < 1e-12
    assert abs(C[1]-(-1))/abs(C[1]+(-1)) < 1e-12
    print("test_get_V_tilde_matrix_dp ok")