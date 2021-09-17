from psage.modform.maass.automorphic_forms import AutomorphicFormSpace

from classes.approx_modform import get_approxmodform_basis

load("../petersson_product.sage")

def run_unit_tests_eisenstein():
    test_compute_eisenstein_series()

def test_compute_eisenstein_series():
    S = AutomorphicFormSpace(Gamma0(2),8)
    cuspforms = get_approxmodform_basis(S,50)
    modforms = get_approxmodform_basis(S,50,modform_type="ModForm")
    eisforms = compute_eisenstein_series(cuspforms,modforms)
    M = ModularForms(2,8,prec=10).eisenstein_subspace()
    c_1_correct = M.q_echelon_basis()[1][2]
    c_1 = eisforms[1].get_cusp_expansion(Cusp(1,0),trunc_order=10)[2]
    assert abs(c_1-c_1_correct) < 1e-10
    print("test_compute_eisenstein_series ok")
