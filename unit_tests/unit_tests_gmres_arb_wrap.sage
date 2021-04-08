from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from psage.groups.permutation_alg import MyPermutation
from psage.modform.arithgroup.mysubgroup import MySubgroup

from iterative_solvers.gmres_arb_wrap import gmres_mgs_arb_wrap
from point_matching.point_matching_arb_wrap import digits_to_bits, get_V_tilde_matrix_b_arb_wrap, get_coefficients_gmres_arb_wrap

def run_unit_tests_gmres_arb_wrap():
    test_gmres_non_factored()
    test_get_coefficients_gmres_arb_wrap()

def test_gmres_non_factored():
    """
    Tests GMRES for the case where the pm-matrix is not factored and where an approximate
    inverse is used as a preconditioner
    """
    digit_prec = 50
    S = AutomorphicFormSpace(Gamma0(1),weight=12)
    bit_prec = digits_to_bits(digit_prec)
    RBF = RealBallField(bit_prec)
    Y = RBF(S.group().minimal_height()*0.8)
    M = 30
    V, b = get_V_tilde_matrix_b_arb_wrap(S,M,Y,bit_prec)

    tol = RBF(10.0)**(-digit_prec)
    low_prec = 64
    V_inv = V.get_inv(low_prec)
    x_gmres_arb_wrap = gmres_mgs_arb_wrap(V, b, bit_prec, tol, M=V_inv, maxiter=10)

    res = x_gmres_arb_wrap[0]._get_mcbd(bit_prec)

    assert abs(res[0][0]+24) < RBF(1e-12)

    print("test_gmres_non_factored ok")

def test_get_coefficients_gmres_arb_wrap():
    S = AutomorphicFormSpace(Gamma0(8),4) #This example contains some empty coordinate lists
    C = get_coefficients_gmres_arb_wrap(S,35)._get_mcbd(53)
    assert abs(C[0][0]) < RBF(1e-12)
    assert abs(C[1][0]-(-4))/abs(C[1][0]+(-4)) < RBF(1e-12)

    print("test_get_coefficients_gmres_arb_wrap ok")
    