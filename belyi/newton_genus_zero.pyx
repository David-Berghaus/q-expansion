import math

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
from classes.factored_polynomial import Factored_Polynomial
from point_matching.point_matching_arb_wrap import get_pi_ball, get_coefficients_haupt_ir_arb_wrap, digits_to_bits

# Possible optimizations:
#     -use lower-level functions
#     -optimize construction of polynomials
#     -set approximate zeros in J to zero for more sparsity
#     -implement approximate polynomial functionality for potentially less precision loss

cpdef locate_coset_in_permT(int coset_index, permT):
    """
    Locates position of cycle that contains coset_index in permT.
    """
    for i in range(len(permT)):
        cycle = permT[i]
        if coset_index in cycle:
            return i

cpdef get_ell_2_points(G, bit_prec):
    """
    Computes approximations of all elliptic points of permutation of order two of MySubgroup G to bit_prec precision as a ComplexBallField.
    Returns a list of tuples where each tuple has entries ([(ell_point, nearest_cusp_index)], order)
    """
    cusps = G.cusps()
    coset_reps = G.coset_reps()
    CBF = ComplexBallField(bit_prec+16)
    permS = G.permS.cycle_tuples()
    permT = G.permT.cycle_tuples()
    fund_ell_2_point = CBF(0,1) #ell_2_point of modular group
    ell_points = []
    for perm in permS:
        order = len(perm)
        coset_index = perm[0] #We could try to locate the coset with the smallest cusp-width to get a bit more precision
        cusp_index = locate_coset_in_permT(coset_index, permT)
        ell_map = (G.cusp_normalizer(cusps[cusp_index])).inverse()*coset_reps[coset_index-1]
        ell_point = apply_moebius_transformation_arb_wrap(fund_ell_2_point,ell_map[0],ell_map[1],ell_map[2],ell_map[3])
        tuple_index = get_tuple_list_index(ell_points, 1, order)
        if tuple_index == None: #We havent considered an elliptic point of this order yet
            ell_points.append( ([(ell_point, cusp_index)],order) )
        else:
            ell_points[tuple_index][0].append((ell_point, cusp_index))
    return ell_points

cpdef get_ell_3_points(G, bit_prec):
    """
    Computes approximations of all elliptic points of permutation of order three of MySubgroup G to bit_prec precision as a ComplexBallField.
    Returns a list of tuples where each tuple has entries ([(ell_point, nearest_cusp_index)], order)
    """
    cusps = G.cusps()
    coset_reps = G.coset_reps()
    CBF = ComplexBallField(bit_prec+16)
    RBF = RealBallField(bit_prec+16)
    permR = G.permR.cycle_tuples()
    permT = G.permT.cycle_tuples()
    fund_ell_3_point = CBF(-1,RBF(3).sqrt())/2 #ell_3_point of modular group
    ell_points = []
    for perm in permR:
        order = len(perm)
        coset_index = perm[0] #We could try to locate the coset with the smallest cusp-width to get a bit more precision
        cusp_index = locate_coset_in_permT(coset_index, permT)
        ell_map = (G.cusp_normalizer(cusps[cusp_index])).inverse()*coset_reps[coset_index-1]
        ell_point = apply_moebius_transformation_arb_wrap(fund_ell_3_point,ell_map[0],ell_map[1],ell_map[2],ell_map[3])
        tuple_index = get_tuple_list_index(ell_points, 1, order)
        if tuple_index == None: #We havent considered an elliptic point of this order yet
            ell_points.append( ([(ell_point, cusp_index)],order) )
        else:
            ell_points[tuple_index][0].append((ell_point, cusp_index))
    return ell_points

cdef eval_principal_hauptmodul(acb_mat_t principal_coeffs, ComplexBall q, int bit_prec):
    """
    Evaluates hauptmodul at point z. principal_coeffs is given as a Mx1 matrix starting with c_1.
    We assume normalization J(z) = 1/q + 0 + c_1*q + ...
    """
    M = acb_mat_nrows(principal_coeffs)
    CBF = ComplexBallField(bit_prec+16)
    RBF = RealBallField(bit_prec+16)
    cdef ComplexBall res = CBF(0,0)
    cdef Acb_Mat q_vec = Acb_Mat(M,1)
    cdef int i = 0
    acb_approx_set(acb_mat_entry(q_vec.value,0,0), q.value)
    for i in range(1,M):
        acb_approx_mul(acb_mat_entry(q_vec.value,i,0), acb_mat_entry(q_vec.value,i-1,0), q.value, bit_prec)
    acb_mat_approx_dot(res.value, principal_coeffs, q_vec.value, bit_prec) #We could have also used Horner's method here but there is currently no approx implementation available
    res += 1/q
    acb_get_mid(res.value, res.value)
    return res

cdef eval_hauptmodul(acb_mat_t haupt_coeffs, ComplexBall q, int bit_prec):
    """
    Evaluates hauptmodul. haupt_coeffs is given as a Mx1 matrix starting with c_0.
    """
    M = acb_mat_nrows(haupt_coeffs)
    CBF = ComplexBallField(bit_prec+16)
    RBF = RealBallField(bit_prec+16)
    cdef ComplexBall res = CBF(0,0)
    cdef Acb_Mat q_vec = Acb_Mat(M,1)
    cdef int i = 0
    acb_one(acb_mat_entry(q_vec.value,0,0))
    acb_approx_set(acb_mat_entry(q_vec.value,1,0), q.value)
    for i in range(2,M):
        acb_approx_mul(acb_mat_entry(q_vec.value,i,0), acb_mat_entry(q_vec.value,i-1,0), q.value, bit_prec)
    acb_mat_approx_dot(res.value, haupt_coeffs, q_vec.value, bit_prec) #We could have also used Horner's method here but there is currently no approx implementation available
    acb_get_mid(res.value, res.value)
    return res

cdef get_non_inf_cusp_values(G, acb_mat_t coeffs, int M, int bit_prec):
    """
    Returns hauptmodul values at cusps that are not i*infinity by returning the values of c_0 of the corresponding cusps.
    Return format is a tuple denoting the value of the cusp and the cusp-width.
    """
    CBF = ComplexBallField(bit_prec+16)
    ncusps = G.ncusps()
    cusps = G.cusps()
    cusp_values = []
    cdef ComplexBall cb_cast
    for i in range(1,ncusps):
        cusp_value = CBF(0,0)
        cb_cast = cusp_value
        acb_set(cb_cast.value, acb_mat_entry(coeffs,i*M,0)) #Set cusp value
        acb_get_mid(cb_cast.value, cb_cast.value)
        width = G.cusp_width(cusps[i])
        pos = get_tuple_list_index(cusp_values, 1, width)
        if pos == None: #We havent considered a cusp with this width yet
            cusp_values.append( ([cusp_value],width) )
        else:
            cusp_values[pos][0].append(cusp_value)
    return cusp_values

cdef get_tuple_list_index(tuple_list, tuple_pos, value):
    """
    Given a list of tuples, return index of tuple whose entry at position 'tuple_pos' is equal to 'value'.
    If no such tuple exists, return 'None'.
    """
    for i in range(len(tuple_list)):
        if tuple_list[i][tuple_pos] == value:
            return i
    return None

cdef get_p3_haupt(G, x, Acb_Mat coeffs, int M, int bit_prec):
    """
    Returns polynomials corresponding to o3 in non-factored form as well as the exponents.
    The starting values are obtained by evaluating the hauptmodul at the elliptic points.
    """
    CBF = ComplexBallField(bit_prec+16)
    pi = get_pi_ball(bit_prec)
    cusps = G.cusps()
    ell_3_points = get_ell_3_points(G, bit_prec)
    haupt_values = []
    cdef Acb_Mat_Win acb_mat_win_cast
    for i in range(len(ell_3_points)):
        ell_3_point_tuple = ell_3_points[i]
        order = ell_3_point_tuple[1]
        haupt_values_order = [] #All values of the hauptmodul at current order
        for (ell_point, nearest_cusp_index) in ell_3_point_tuple[0]:
            cusp_width = G.cusp_width(cusps[nearest_cusp_index])
            q = (2*pi*CBF(0,1)*ell_point/cusp_width).exp()
            haupt_coeffs = coeffs.get_window(nearest_cusp_index*M,0,nearest_cusp_index*M+M,1)
            acb_mat_win_cast = haupt_coeffs
            if nearest_cusp_index == 0: #We are considering the principal cusp
                haupt_values_order.append(eval_principal_hauptmodul(acb_mat_win_cast.value, q, bit_prec))
            else:
                haupt_values_order.append(eval_hauptmodul(acb_mat_win_cast.value, q, bit_prec))
        haupt_values.append( (haupt_values_order, order) )
    p3 = Factored_Polynomial(x,haupt_values)
    return p3

cdef get_p2_haupt(G, x, Acb_Mat coeffs, int M, int bit_prec):
    """
    Returns polynomials corresponding to o2 in non-factored form as well as the exponents.
    The starting values are obtained by evaluating the hauptmodul at the elliptic points.
    """
    CBF = ComplexBallField(bit_prec+16)
    pi = get_pi_ball(bit_prec)
    cusps = G.cusps()
    ell_2_points = get_ell_2_points(G, bit_prec)
    haupt_values = []
    cdef Acb_Mat_Win acb_mat_win_cast
    for i in range(len(ell_2_points)):
        ell_2_point_tuple = ell_2_points[i]
        order = ell_2_point_tuple[1]
        haupt_values_order = [] #All values of the hauptmodul at current order
        for (ell_point, cusp_index) in ell_2_point_tuple[0]:
            cusp_width = G.cusp_width(cusps[cusp_index])
            q = (2*pi*CBF(0,1)*ell_point/cusp_width).exp()
            haupt_coeffs = coeffs.get_window(cusp_index*M,0,cusp_index*M+M,1)
            acb_mat_win_cast = haupt_coeffs
            if cusp_index == 0: #We are considering the principal cusp
                haupt_values_order.append(eval_principal_hauptmodul(acb_mat_win_cast.value, q, bit_prec))
            else:
                haupt_values_order.append(eval_hauptmodul(acb_mat_win_cast.value, q, bit_prec))
        haupt_values.append( (haupt_values_order, order) )
    p2 = Factored_Polynomial(x,haupt_values)
    return p2

cdef get_pc_haupt(G, x, Acb_Mat coeffs, int M, int bit_prec):
    """
    Returns polynomials corresponding to o2 in non-factored form as well as the exponents.
    The starting values are obtained by evaluating the hauptmodul at the cusps.
    """
    haupt_values = get_non_inf_cusp_values(G, coeffs.value, M, bit_prec)
    pc = Factored_Polynomial(x,haupt_values)
    return pc

cpdef get_f(factored_polynomials):
    """
    Compute f = p3-p2-1728*pc.
    factored_polynomials is a tuple (p_3, p_2, p_c) where p_i are objects of type Factored_Polynomial.
    """
    f = 0
    for poly_type in range(3): #Bad wording but poly_type refers to p_3, p_2 or p_c
        factored_polynomial = factored_polynomials[poly_type]
        p = factored_polynomial.construct()
        if poly_type == 1: #We are considering p_2
            p = -p
        elif poly_type == 2: #We are considering p_c
            p *= -1728
        f += p
    return f

cpdef get_jacobian(factored_polynomials, G):
    """
    Return Jacobi-matrix.
    factored_polynomials is a tuple (p_3, p_2, p_c) where p_i are objects of type Factored_Polynomial.
    """
    cdef int N = G.index()+1
    cdef Acb_Mat J = Acb_Mat(N,N)
    cdef Polynomial_complex_arb derivative
    cdef int i, j, index
    i, j, index = 0, 0, 0
    for poly_type in range(3): #Bad wording but poly_type refers to p_3, p_2 or p_c
        factored_polynomial = factored_polynomials[poly_type]
        for poly_index in range(len(factored_polynomial.factors)):
            amount_of_unknowns = factored_polynomial.factors[poly_index][0].degree() #Note that the last monomial is of the form 1*x^n, so we got n unknowns
            for coeff_index in range(amount_of_unknowns):
                derivative = factored_polynomial.derivative(poly_index, coeff_index)
                if poly_type == 1: #We are considering p_2
                    derivative = -derivative
                elif poly_type == 2: #We are considering p_c
                    derivative *= -1728
                for i in range(derivative.degree()+1):
                    acb_swap(acb_mat_entry(J.value,i,j), acb_poly_get_coeff_ptr(derivative.__poly,i)) #We dont need to keep poly so we can just swap
                if coeff_index == amount_of_unknowns-1: #Take care of boundary condition at infinity to fill last row
                    multiplicity = factored_polynomial.factors[poly_index][1]
                    if poly_type == 0: #We are considering p_3
                        acb_set_si(acb_mat_entry(J.value,N-1,j),multiplicity)
                    elif poly_type == 2: #We are considering p_c
                        acb_set_si(acb_mat_entry(J.value,N-1,j),-multiplicity)
                j += 1
    return J

cpdef get_unknown_coeff_column(factored_polynomials, int N):
    """
    Get a column vector containing the unknown coefficients.
    """
    cdef Polynomial_complex_arb poly_cast
    cdef Acb_Mat coeffs = Acb_Mat(N,1)
    cdef int i, j
    i = 0
    for factored_polynomial in factored_polynomials:
        for (poly_fact, _) in factored_polynomial.factors:
            poly_cast = poly_fact
            j = 0
            for j in range(acb_poly_degree(poly_cast.__poly)):
                acb_set(acb_mat_entry(coeffs.value,i,0), acb_poly_get_coeff_ptr(poly_cast.__poly,j))
                i += 1
    return coeffs

cdef newton_step(acb_mat_t res, factored_polynomials, G, int bit_prec):
    cdef Polynomial_complex_arb f = get_f(factored_polynomials)
    cdef Acb_Mat J = get_jacobian(factored_polynomials, G)
    N = acb_mat_nrows(J.value)
    cdef Acb_Mat x = get_unknown_coeff_column(factored_polynomials, N)
    cdef Acb_Mat f_x = Acb_Mat(N, 1)
    acb_mat_set_poly(f_x.value, f.__poly)

    #We are now in a position to take the boundary condition at infinity into account
    #to determine the last entry of f_x
    CBF = ComplexBallField(bit_prec)
    cdef ComplexBall condition_at_infinity = CBF(0,0)
    p3 = factored_polynomials[0]
    for (factor, order) in p3.factors:
        condition_at_infinity += order*factor.coefficients()[-2]
    pc = factored_polynomials[2]
    for (factor, order) in pc.factors:
        condition_at_infinity -= order*factor.coefficients()[-2]
    condition_at_infinity -= 744
    acb_set(acb_mat_entry(f_x.value,N-1,0),condition_at_infinity.value)

    cdef Acb_Mat update = Acb_Mat(N, 1)
    acb_mat_approx_solve(update.value, J.value, f_x.value, bit_prec)
    acb_mat_approx_sub(x.value, x.value, update.value, bit_prec)
    return x
            
cpdef test_jacobian(S, digit_prec):
    G = S.group()
    cdef Acb_Mat c
    c, M = get_coefficients_haupt_ir_arb_wrap(S,digit_prec,only_principal_expansion=False,return_M=True)
    bit_prec = digits_to_bits(2*digit_prec)
    CBF = ComplexBallField(bit_prec)
    cdef Polynomial_complex_arb x = Polynomial_complex_arb(CBF['x'], is_gen=True)
    # cdef Polynomial_complex_arb p1 = Polynomial_complex_arb(CBF['x'], [CBF(1,1), 0, 1])

    p2 = get_p2_haupt(G, x, c, M, bit_prec)
    p3 = get_p3_haupt(G, x, c, M, bit_prec)
    pc = get_pc_haupt(G, x, c, M, bit_prec)
    factored_polynomials = (p3, p2, pc)

    cdef Acb_Mat res = Acb_Mat(G.index()+1,1)
    newton_step(res.value, factored_polynomials, G, bit_prec)
    #J = get_jacobian(factored_polynomials, G)

cpdef test(S, digit_prec):
    G = S.group()
    cdef Acb_Mat c
    c, M = get_coefficients_haupt_ir_arb_wrap(S,digit_prec,only_principal_expansion=False,return_M=True)
    bit_prec = digits_to_bits(digit_prec)
    CBF = ComplexBallField(bit_prec)
    cdef Polynomial_complex_arb x = Polynomial_complex_arb(CBF['x'], is_gen=True)

    return get_p3_haupt(G, x, c, M, bit_prec)

