from sage.libs.arb.acb_mat cimport *
from sage.rings.complex_arb cimport *
from sage.libs.arb.acb_poly cimport *
from sage.rings.polynomial.polynomial_complex_arb cimport *
from sage.rings.polynomial.polynomial_complex_arb import Polynomial_complex_arb
from sage.rings.complex_arb import ComplexBallField
from sage.rings.real_arb import RealBallField

from arblib_helpers.acb_approx cimport *
from pullback.my_pullback cimport apply_moebius_transformation_arb_wrap
from classes.acb_mat_class cimport Acb_Mat, Acb_Mat_Win
from point_matching.point_matching_arb_wrap import get_pi_ball, get_coefficients_haupt_ir_arb_wrap, digits_to_bits

cpdef get_ell_2_maps(G):
    """
    Get maps to all elliptic points of order two of MySubgroup G as SL2Z_elt.
    """
    ell_2_maps = list()
    permS = G.permS
    fixed_cosets = permS.fixed_elements()
    coset_reps = G.coset_reps()
    for fixed_coset in fixed_cosets:
        ell_2_maps.append(coset_reps[fixed_coset-1]) #Note that the cosets in the permutation are numerated by 1,.. but the reps are stored in a list
    return ell_2_maps

cpdef get_ell_3_maps(G):
    """
    Get maps to all elliptic points of order three of MySubgroup G as SL2Z_elt.
    """
    ell_3_maps = list()
    permR = G.permR
    fixed_cosets = permR.fixed_elements()
    coset_reps = G.coset_reps()
    for fixed_coset in fixed_cosets:
        ell_3_maps.append(coset_reps[fixed_coset-1]) #Note that the cosets in the permutation are numerated by 1,.. but the reps are stored in a list
    return ell_3_maps

cpdef get_ell_2_points(G, bit_prec):
    """
    Gets approximations of all elliptic points of order two of MySubgroup G to bit_prec precision as a ComplexBallField.
    """
    ell_2_points = list()
    ell_2_maps = get_ell_2_maps(G)
    if len(ell_2_maps) != 0:
        CBF = ComplexBallField(bit_prec+16)
        fund_ell_2_point = CBF(0,1) #ell_2_point of modular group
        for ell_2_map in ell_2_maps:
            a, b, c, d = ell_2_map[0], ell_2_map[1], ell_2_map[2], ell_2_map[3]
            ell_2_points.append(apply_moebius_transformation_arb_wrap(fund_ell_2_point,a,b,c,d))
    return ell_2_points

cpdef get_ell_3_points(G, bit_prec):
    """
    Gets approximations of all elliptic points of order three of MySubgroup G to bit_prec precision as a ComplexBallField.
    """
    ell_3_points = list()
    ell_3_maps = get_ell_3_maps(G)
    if len(ell_3_maps) != 0:
        CBF = ComplexBallField(bit_prec+16)
        RBF = RealBallField(bit_prec+16)
        fund_ell_3_point = CBF(RBF(1)/2,RBF(3).sqrt()/2) #ell_3_point of modular group
        for ell_3_map in ell_3_maps:
            a, b, c, d = ell_3_map[0], ell_3_map[1], ell_3_map[2], ell_3_map[3]
            ell_3_points.append(apply_moebius_transformation_arb_wrap(fund_ell_3_point,a,b,c,d))
    return ell_3_points

cpdef get_non_inf_cusp_points(G, bit_prec):
    """
    Gets approximations of all cusp points that are not located at i*infinity.
    """
    cusp_points = list()
    cusps = G.cusps()
    if len(cusps) > 1:
        CBF = ComplexBallField(bit_prec+16)
        RBF = RealBallField(bit_prec+16)
        for i in range(1, len(cusps)):
            c = cusps[i]
            cusp_points.append(CBF(RBF(c.numerator())/c.denominator(),0))
    return cusp_points

cdef eval_hauptmodul(acb_mat_t coeffs, ComplexBall z, int bit_prec):
    """
    Evaluates hauptmodul at point z. Coeffs is given as a Nx1 matrix starting with c_1.
    We assume normalization J(z) = 1/q + 0 + c_1*q + ...
    """
    M = acb_mat_nrows(coeffs)
    CBF = ComplexBallField(bit_prec+16)
    RBF = RealBallField(bit_prec+16)
    cdef ComplexBall res = CBF(0,0)
    cdef ComplexBall q = (2*get_pi_ball(bit_prec)*CBF(0,1)*z).exp()
    cdef Acb_Mat q_vec = Acb_Mat(M,1)
    cdef int i
    acb_approx_set(acb_mat_entry(q_vec.value,0,0), q.value)
    for i in range(1,M):
        acb_approx_mul(acb_mat_entry(q_vec.value,i,0), acb_mat_entry(q_vec.value,i-1,0), q.value, bit_prec)
    acb_mat_approx_dot(res.value, coeffs, q_vec.value, bit_prec) #We could have also used Horner's method here but there is currently no approx implementation available
    res += 1/q
    return res

cpdef test(S, digit_prec):
    G = S.group()
    cdef Acb_Mat_Win c = get_coefficients_haupt_ir_arb_wrap(S,digit_prec)
    bit_prec = digits_to_bits(digit_prec)
    cdef ComplexBall z = get_ell_3_points(G, bit_prec)[0]
    res = eval_hauptmodul(c.value, z, bit_prec)

    CBF = ComplexBallField(bit_prec)
    cdef Polynomial_complex_arb x = Polynomial_complex_arb(CBF['x'], is_gen=True)
    cdef Polynomial_complex_arb p1 = x**2+CBF(1,1)
    cdef ComplexBall tmp = CBF(0,0)
    acb_poly_approx_evaluate_horner(tmp.value, p1.__poly, z.value, bit_prec)
    acb_printd(tmp.value,10)

    return res

