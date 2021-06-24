import numpy as np

from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from psage.groups.permutation_alg import MyPermutation
from psage.modform.arithgroup.mysubgroup import MySubgroup

from point_matching.point_matching_dp import get_coefficients_cuspform_dp, get_V_tilde_matrix_b_cuspform_dp, get_V_tilde_matrix_cuspform_dp, get_coefficients_haupt_dp
from classes.gamma_2_subgroup import Gamma_2_Subgroup
from classes.modform_class import ModForm

def run_unit_tests_point_matching_dp():
    test_get_coefficients_cuspform_dp()
    test_get_V_tilde_matrix_b_cuspform_dp()
    test_get_V_tilde_matrix_cuspform_dp()
    test_get_coefficients_haupt_dp()
    test_get_coefficients_cuspform_gamma_2_subgroup_dp()

def test_get_coefficients_cuspform_dp():
    r='(1 2 5)(3)(4 7 9)(6 10 8)'
    s='(1)(2)(3 4)(5 6)(7 8)(9 10)'
    G=MySubgroup(o2=s,o3=r) #This is a non-congruence subgroup
    S = AutomorphicFormSpace(G,weight=4) #Search for a multiplicity two new-form
    C = get_coefficients_cuspform_dp(S)
    assert abs(C[0][0]-(1.1779944322516185-1.177994432251737j))/abs(C[0][0]+(1.1779944322516185-1.177994432251737j)) < 1e-12
    assert abs(C[1][0]-(2.2534750215544057-5.440369959505566j))/abs(C[1][0]+(2.2534750215544057-5.440369959505566j)) < 1e-12
    print("test_get_coefficients_cuspform_dp ok")

def test_get_V_tilde_matrix_b_cuspform_dp():
    S = AutomorphicFormSpace(Gamma0(8),weight=4)
    V,b = get_V_tilde_matrix_b_cuspform_dp(S,66,0.09742785792574934)
    C = np.linalg.solve(V,b)
    assert abs(C[0][0]) < 1e-12
    assert abs(C[1][0]-(-4))/abs(C[1][0]+(-4)) < 1e-12
    print("test_get_V_tilde_matrix_b_cuspform_dp ok")

def test_get_V_tilde_matrix_cuspform_dp():
    S = AutomorphicFormSpace(Gamma0(11),weight=2)
    T = get_V_tilde_matrix_cuspform_dp(S,82,0.07085662394599952)
    V = T[1:,1:]
    b = -T[1:,0]
    C = np.linalg.solve(V,b)
    assert abs(C[0]-(-2))/abs(C[0]+(-2)) < 1e-12
    assert abs(C[1]-(-1))/abs(C[1]+(-1)) < 1e-12
    print("test_get_V_tilde_matrix_cuspform_dp ok")

def test_get_coefficients_haupt_dp():
    S = AutomorphicFormSpace(Gamma0(4),weight=0)
    C = get_coefficients_haupt_dp(S)
    assert abs(C[0][0]-20)/abs(C[0][0]+20) < 1e-12
    print("test_get_coefficients_haupt_dp ok")

def test_get_coefficients_cuspform_gamma_2_subgroup_dp(): #We only test Gamma(2) subgroups for dp since it should be the same for arb wrap
    o_0 = MyPermutation('(1 3 6 2 5 4)')
    o_inf = MyPermutation('(1 2 3 4 5 6)')
    G = Gamma_2_Subgroup(o_0,o_inf)
    S = ModForm(G,2)
    c1 = get_coefficients_cuspform_dp(S,Y=0.07)
    c2 = get_coefficients_cuspform_dp(S,Y=0.06)
    assert abs(c1[0][0]-c2[0][0])/abs(c1[0][0]+c2[0][0]) < 1e-10
    print("test_get_coefficients_cuspform_gamma_2_subgroup_dp ok")