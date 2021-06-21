from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from psage.groups.permutation_alg import MyPermutation
from psage.modform.arithgroup.mysubgroup import MySubgroup

from point_matching.point_matching_arb_wrap import get_coefficients_cuspform_ir_arb_wrap

def run_unit_tests_ir_arb_wrap():
    test_get_coefficients_cuspform_ir_arb_wrap()
    test_get_coefficients_cuspform_ir_arb_wrap2()

def test_get_coefficients_cuspform_ir_arb_wrap():
    S = AutomorphicFormSpace(Gamma0(8),4) #This example contains some empty coordinate lists
    C = get_coefficients_cuspform_ir_arb_wrap(S,35)._get_mcbd(53)
    assert abs(C[0][0]) < RBF(1e-12)
    assert abs(C[1][0]-(-4))/abs(C[1][0]+(-4)) < RBF(1e-12)

    print("test_get_coefficients_cuspform_ir_arb_wrap ok")

def test_get_coefficients_cuspform_ir_arb_wrap2():
    r='(1 2 5)(3)(4 7 9)(6 10 8)'
    s='(1)(2)(3 4)(5 6)(7 8)(9 10)'
    G=MySubgroup(o2=s,o3=r) #This is a non-congruence subgroup
    S = AutomorphicFormSpace(G,weight=4) #Search for a multiplicity two new-form
    C = get_coefficients_cuspform_ir_arb_wrap(S,35)._get_mcbd(53)
    assert abs(float(C[0][0].real())+float(C[0][0].imag())*1j-(1.1779944322516185-1.177994432251737j))/abs(float(C[0][0].real())+float(C[0][0].imag())*1j+(1.1779944322516185-1.177994432251737j)) < 1e-12
    print("test_get_coefficients_cuspform_ir_arb_wrap for multiplicity two non-congruence ok")