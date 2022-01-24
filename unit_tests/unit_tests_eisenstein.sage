from psage.modform.maass.automorphic_forms import AutomorphicFormSpace

from classes.fourier_expansion import get_cuspform_q_expansion_approx, get_modform_q_expansion_approx, get_cuspform_basis_approx, get_modform_basis_approx
from eisenstein.haberland import compute_petersson_product_haberland

load("../eisenstein/eisenstein.sage")
load("../eisenstein/nelson_collins.sage")

def run_unit_tests_eisenstein():
    test_compute_eisenstein_series()
    test_petersson_product()
    test_petersson_product2()

def test_compute_eisenstein_series():
    S = AutomorphicFormSpace(Gamma0(2),8)
    cuspforms = get_cuspform_basis_approx(S,50)
    modforms = get_modform_basis_approx(S,50)
    eisforms = compute_eisenstein_series(cuspforms,modforms)
    M = ModularForms(2,8,prec=10).eisenstein_subspace()
    c_1_correct = M.q_echelon_basis()[1][2]
    c_1 = eisforms[1].get_cusp_expansion(Cusp(1,0),trunc_order=10)[2]
    assert abs(c_1-c_1_correct) < 1e-10
    print("test_compute_eisenstein_series ok")

def test_petersson_product():
    S = AutomorphicFormSpace(Gamma0(4),6)
    f, g = get_cuspform_q_expansion_approx(S,50), get_modform_q_expansion_approx(S,50)
    petersson_1 = compute_petersson_product_haberland(f,g)
    petersson_2 = compute_petersson_product_nelson_collins(f,g)
    assert abs(petersson_1-petersson_2) < 1e-17
    print("test_petersson_product ok")

def test_petersson_product2(): #Also test petersson for noncongruence and cusp-width larger 1 at infinity
    U1 = MySubgroup(o2='(1 2)(3 4)(5)(6 7)',o3='(1)(2 3 5)(4 6 7)')
    S = AutomorphicFormSpace(U1,4)
    f, g = get_cuspform_q_expansion_approx(S,50), get_modform_q_expansion_approx(S,50)
    petersson_1 = compute_petersson_product_haberland(f,g)
    petersson_2 = compute_petersson_product_nelson_collins(f,g)
    assert abs(petersson_1-petersson_2) < 1e-21
    print("test_petersson_product2 ok")