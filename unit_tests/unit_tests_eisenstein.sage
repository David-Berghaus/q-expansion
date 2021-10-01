from psage.modform.maass.automorphic_forms import AutomorphicFormSpace

from classes.approx_modform import get_approxmodform_basis

load("../eisenstein.sage")
load("../haberland.sage")
load("../nelson_collins.sage")

def run_unit_tests_eisenstein():
    test_compute_eisenstein_series()
    test_petersson_product()

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

def test_petersson_product():
    S = AutomorphicFormSpace(Gamma0(4),6)
    f, g = ApproxModForm(S,50), ApproxModForm(S,50,modform_type="ModForm")
    petersson_1 = compute_petersson_product_haberland(f,g)
    petersson_2 = compute_petersson_product_nelson_collins(f._get_cusp_expansions_dict(),g._get_cusp_expansions_dict())
    assert abs(petersson_1-petersson_2) < 1e-17
    print("test_petersson_product ok")