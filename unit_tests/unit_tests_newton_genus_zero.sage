from psage.modform.maass.automorphic_forms import AutomorphicFormSpace

from belyi.newton_genus_zero import run_newton

def run_unit_tests_newton_genus_zero():
    test_newton_genus_zero()
    
def test_newton_genus_zero():
    S = AutomorphicFormSpace(Gamma0(6),0)
    factored_polynomials = run_newton(S,10,500,stop_when_coeffs_are_recognized=False)
    p = factored_polynomials[0].construct()
    assert ((p.coefficients()[0]).abs()-4102915888729).abs() < RBF(10)**(-163)
    print("test_newton_genus_zero ok")