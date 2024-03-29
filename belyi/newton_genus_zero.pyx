# ****************************************************************************
#       Copyright (C) 2022 David Berghaus <berghaus@th.physik.uni-bonn.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import math

from sage.libs.arb.acb_mat cimport *
from sage.rings.complex_arb cimport *
from sage.libs.arb.acb_poly cimport *
from sage.rings.polynomial.polynomial_complex_arb cimport *
from sage.rings.polynomial.polynomial_complex_arb import Polynomial_complex_arb
from sage.rings.complex_arb import ComplexBallField
from sage.rings.real_arb import RealBallField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.modular.cusps import Cusp

from arblib_helpers.acb_approx cimport *
from pullback.my_pullback cimport apply_moebius_transformation_arb_wrap
from classes.acb_mat_class cimport Acb_Mat, Acb_Mat_Win
from belyi.number_fields import get_decimal_digit_prec, is_effectively_zero, get_numberfield_of_coeff
from belyi.expression_in_u_and_v import convert_from_Kv_to_Kw
from classes.factored_polynomial import Factored_Polynomial, recognize_coeffs_using_u, get_improved_choice_of_u_interior_Kv, get_updated_Kw_v_Kw, update_terms_linear_in_u_to_new_u
from classes.fourier_expansion import get_hauptmodul_q_expansion_approx
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

cpdef get_CBF_list(CBF, int N):
    """
    Returns a list of length N consisting of CBF(0,0).
    """
    cdef int i
    cbf_list = [CBF(0,0) for i in range(N)]
    return cbf_list

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

cdef get_non_inf_cusp_values(G, acb_mat_t coeffs, int M, int bit_prec, return_cusp_rep_values=False):
    """
    Returns hauptmodul values at cusps that are not i*infinity by returning the values of c_0 of the corresponding cusps.
    Return format is a tuple denoting the evaluation at the cusp and the cusp-width.
    If "return_cusp_rep_values==True" we also return the cusp representatives corresponding to the evaluations.
    This useful to identify the cusp corresponding to a root of pc later.
    """
    CBF = ComplexBallField(bit_prec+16)
    ncusps = G.ncusps()
    cusps = G.cusps()
    cusp_values = [] #This list collects cusp evaluations for each width
    cusp_rep_values = [] #This list contains entries of the form (cusp,cusp_evaluation)
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
        cusp_rep_values.append( (cusps[i],cusp_value) )

    if return_cusp_rep_values == False:
        return cusp_values
    else:
        return cusp_values, cusp_rep_values

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
    p3 = Factored_Polynomial(x,root_tuples=haupt_values)
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
    p2 = Factored_Polynomial(x,root_tuples=haupt_values)
    return p2

cdef get_pc_haupt(G, x, Acb_Mat coeffs, int M, int bit_prec, return_cusp_rep_values=False):
    """
    Returns polynomials corresponding to o2 in non-factored form as well as the exponents.
    The starting values are obtained by evaluating the hauptmodul at the cusps.
    """
    haupt_values, cusp_rep_values = get_non_inf_cusp_values(G, coeffs.value, M, bit_prec, return_cusp_rep_values=True)
    pc = Factored_Polynomial(x,root_tuples=haupt_values)
    if return_cusp_rep_values == False:
        return pc
    else:
        return pc, cusp_rep_values

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

def get_empty_coeff_tuples_list(factored_polynomials):
    """
    For each factored polynomial, return a list that is of the form [([c_0,c_1,...],order)] and set all c_i to "None" (except for the last which is always 1).
    """
    coeff_tuples_list = []
    for factored_polynomial in factored_polynomials:
        coeff_tuples = []
        for (poly_fact, order) in factored_polynomial.factors:
            amount_of_unknowns = poly_fact.degree() #Note that the last monomial is of the form 1*x^n, so we got n unknowns
            coeffs = [None for _ in range(amount_of_unknowns)]
            coeffs.append(1) #The highest order coefficient is always 1
            coeff_tuples.append( (coeffs,order) )
        coeff_tuples_list.append(coeff_tuples)
    return coeff_tuples_list

cdef get_coeff_tuples_from_coeff_column(acb_mat_t x, factored_polynomials, int bit_prec, coeff_tuples_list, swap=False):
    """
    Convert column-vector x, containing all unknown coefficients, into coeff_tuples for each polynomial.
    """
    CBF = ComplexBallField(bit_prec)
    cdef ComplexBall cb_cast
    if coeff_tuples_list == None:
        coeff_tuples_list = get_empty_coeff_tuples_list(factored_polynomials)
    x_index = 0
    for i in range(len(factored_polynomials)):
        coeff_tuples = coeff_tuples_list[i]
        for (coeffs, order) in coeff_tuples:
            for j in range(len(coeffs)-1): #Recall that the last coeff is 1
                if coeffs[j] == None or isinstance(coeffs[j].parent(),ComplexBallField) == True: #We need to update the value of this coefficient
                    cb_cast = CBF(0)
                    #It is important that we use approximate functions here because otherwise
                    #the error-bounds would prevent us from casting to higher precision
                    if swap == False:
                        acb_approx_set(cb_cast.value, acb_mat_entry(x,x_index,0))
                    else:
                        acb_approx_swap(cb_cast.value, acb_mat_entry(x,x_index,0))
                    coeffs[j] = cb_cast
                x_index += 1
    return coeff_tuples_list

cdef newton_step(factored_polynomials, G, int bit_prec, coeff_tuples_list=None):
    """
    Perform one step of the newton iteration.
    factored_polynomials is used to create a Jacobi matrix and to refine the precision of the coefficients.
    Returns a list consisting of the (refined) coeff_tuples of p_3, p_2, p_c. 
    """
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
        condition_at_infinity += order*factor[factor.degree()-1]
    pc = factored_polynomials[2]
    for (factor, order) in pc.factors:
        condition_at_infinity -= order*factor[factor.degree()-1]
    if G.cusp_width(Cusp(1,0)) == 1: #Otherwise the last equation is equal to zero so we don't have to change anything
        condition_at_infinity -= 744
    acb_set(acb_mat_entry(f_x.value,N-1,0),condition_at_infinity.value)

    cdef Acb_Mat update = Acb_Mat(N, 1)
    acb_mat_approx_solve(update.value, J.value, f_x.value, bit_prec)
    acb_mat_approx_sub(x.value, x.value, update.value, bit_prec)
    coeff_tuples_list = get_coeff_tuples_from_coeff_column(x.value, factored_polynomials, bit_prec, coeff_tuples_list, swap=True)
    return coeff_tuples_list

def get_simplest_non_zero_coeff(p3, p2, pc, digit_prec, coeff_shift=-1):
    """
    Tries to find a non-zero coefficient that is comparatively easy to identify which can be used to define the numberfield of the Belyi map.
    This is usually c_{N-1} (i.e., the second leading coefficient which is denoted by coeff_shift=-1) of the smallest factored polynomial (if it is unequal to zero).
    """
    p_smallest_deg = p3.get_smallest_degree_poly()
    c = p_smallest_deg[p_smallest_deg.degree()+coeff_shift]
    if p_smallest_deg.degree() != 1 or is_effectively_zero(c,digit_prec) == True: #Try to find a smaller poly or a coeff that is non-zero
        tmp = p2.get_smallest_degree_poly()
        c_tmp = tmp[tmp.degree()+coeff_shift]
        if is_effectively_zero(c,digit_prec) == True: #We cannot work with a zero coeff so we definitely want to switch
            p_smallest_deg = tmp
            c = c_tmp
        elif tmp.degree() < p_smallest_deg.degree() and is_effectively_zero(c_tmp,digit_prec) == False:
            p_smallest_deg = tmp
            c = c_tmp
    if p_smallest_deg.degree() != 1 or is_effectively_zero(c,digit_prec) == True: #Try to find a smaller poly or a coeff that is non-zero
        tmp = pc.get_smallest_degree_poly()
        if tmp != None: #If we only have one cusp then pc is empty
            c_tmp = tmp[tmp.degree()+coeff_shift]
            if is_effectively_zero(c,digit_prec) == True: #We cannot work with a zero coeff so we definitely want to switch
                p_smallest_deg = tmp
                c = c_tmp
            elif tmp.degree() < p_smallest_deg.degree() and is_effectively_zero(c_tmp,digit_prec) == False:
                p_smallest_deg = tmp
                c = c_tmp
    if is_effectively_zero(c,digit_prec) == True: #All coefficients with specified coeff_shift are zero
        return None
    return c

def get_expression_to_recognize(p3, p2, pc, digit_prec, cusp_width, index, G, max_extension_field_degree):
    """
    Return an expression that is linear in u and easy to recognize. 
    We usually choose the expression to be one of the second-leading order coefficients. 
    If these are all zero, we choose u to be a quotient of coefficients.
    """
    res = get_simplest_non_zero_coeff(p3,p2,pc,digit_prec,coeff_shift=-1)
    if res != None:
        return res
    if cusp_width == 1: #u is in QQ(v) so we can choose any coefficient
        for i in range(2,index):
            coeff_shift = -i
            res = get_simplest_non_zero_coeff(p3,p2,pc,digit_prec,coeff_shift=coeff_shift)
            if res != None:
                return res
        raise ArithmeticError("Found no suitable expression to recognize u!")
    else: #We have to construct quotients of coefficients
        for i in range(2,index):
            coeff_shift = -i
            den = get_simplest_non_zero_coeff(p3,p2,pc,digit_prec,coeff_shift=coeff_shift)
            if den != None:
                num = get_simplest_non_zero_coeff(p3,p2,pc,digit_prec,coeff_shift=coeff_shift-1)
                if num == None: #Consecutive coefficients are all zero which means that we cannot use a quotient to obtain an expression linear in u
                    #Try to recognize a coefficient that depends on u^2, otherwise give up
                    num = get_simplest_non_zero_coeff(p3,p2,pc,digit_prec,coeff_shift=-2)
                    if num != None:
                        res = num.sqrt()
                    else:
                        raise ArithmeticError("We have not considered this case yet!")
                else:
                    res = num/den
                return res
        if G.is_congruence() and max_extension_field_degree == 1: #We are safe here because we can set u = 1
            return 1
        raise ArithmeticError("Found no suitable expression to recognize u!")

def recognize_coeff_tuples_list(coeff_tuples_list, Kv, u, v_Kw, estimated_bit_prec, max_u_power=None):
    """
    Tries to recognize algebraic expressions of coeff_tuples_list.
    Note that the function is mutating the argument (i.e., the CBFs in coeff_tuples_list are replaced by the corresponding algebraic values).
    This function returns True if all coefficients up to max_u_power have been recognized and False otherwise.
    """
    have_all_coeffs_lists_been_recognized = True
    for coeff_tuples in coeff_tuples_list:
        for (coeffs,_) in coeff_tuples:
            coeffs, have_all_coeffs_been_recognized = recognize_coeffs_using_u(coeffs,Kv,u,v_Kw,estimated_bit_prec=estimated_bit_prec,max_u_power=max_u_power)
            if have_all_coeffs_been_recognized == False: #Found an invalid example so we need an additional newton iteration
                have_all_coeffs_lists_been_recognized = False
    return have_all_coeffs_lists_been_recognized

cpdef newton(factored_polynomials, G, int curr_bit_prec, int target_bit_prec, stop_when_coeffs_are_recognized, max_extension_field_degree=None):
    coeff_tuples_list = None
    u, Kv, Kw, v_Kw, u_interior_Kv = None, None, None, None, None #Variables that will be used if stop_when_coeffs_are_recognized == True
    has_u_been_improved = False

    while curr_bit_prec < target_bit_prec:
        coeff_tuples_list = newton_step(factored_polynomials, G, curr_bit_prec, coeff_tuples_list=coeff_tuples_list)
        if 2*curr_bit_prec < target_bit_prec:
            curr_bit_prec *= 2
            CBF = ComplexBallField(curr_bit_prec)
            polynomial_ring = PolynomialRing(CBF,"x")
            x = polynomial_ring.gen()
            p3 = Factored_Polynomial(x,coeff_tuples=coeff_tuples_list[0])
            p2 = Factored_Polynomial(x,coeff_tuples=coeff_tuples_list[1])
            pc = Factored_Polynomial(x,coeff_tuples=coeff_tuples_list[2])
            factored_polynomials = (p3, p2, pc)
        else: #If we are at our last iteration we don't need to increase the precision again before constructing p
            CBF = ComplexBallField(curr_bit_prec)
            polynomial_ring = PolynomialRing(CBF,"x")
            x = polynomial_ring.gen()
            p3 = Factored_Polynomial(x,coeff_tuples=coeff_tuples_list[0])
            p2 = Factored_Polynomial(x,coeff_tuples=coeff_tuples_list[1])
            pc = Factored_Polynomial(x,coeff_tuples=coeff_tuples_list[2])
            factored_polynomials = (p3, p2, pc)
            curr_bit_prec *= 2
        coeff_prec = get_coeff_min_precision(factored_polynomials,G.index()+1)
        coeff_bit_prec = digits_to_bits(coeff_prec)
        print("Estimated digit prec: ", coeff_prec)

        if stop_when_coeffs_are_recognized == True: #Try to recognize coefficients as algebraic numbers
            if max_extension_field_degree == None:
                raise ArithmeticError("Please specify the maximal degree of the extension field of the Belyi map!")
            if u == None: #No algebraic coefficients have been recognized yet
                principal_cusp_width = G.cusp_width(Cusp(1,0))
                coeff_fl = get_expression_to_recognize(p3,p2,pc,coeff_prec,principal_cusp_width,G.index(),G,max_extension_field_degree)
                tmp = get_numberfield_of_coeff(coeff_fl,max_extension_field_degree,principal_cusp_width,estimated_bit_prec=coeff_bit_prec)
                if tmp == None: #Failed to recognize coeffs as alg numbers
                    continue
                Kv, Kw, v_Kw, u_interior_Kv = tmp
                if principal_cusp_width == 1:
                    u = convert_from_Kv_to_Kw(u_interior_Kv,v_Kw)
                else:
                    u = Kw.gen()

            #If we are here, u has been recognized
            #Only try to recognize terms linear in u first
            have_all_coeffs_lists_been_recognized = recognize_coeff_tuples_list(coeff_tuples_list,Kv,u,v_Kw,coeff_bit_prec,max_u_power=1) #Note that this function might also mutate coeff_tuples_list
            #Afterwards use results to potentially improve u
            if have_all_coeffs_lists_been_recognized == True and has_u_been_improved == False:
                u_interior_Kv_updated = get_improved_choice_of_u_interior_Kv(coeff_tuples_list,u_interior_Kv,principal_cusp_width)
                if principal_cusp_width == 1:
                    u_interior_Kv = u_interior_Kv_updated
                    u = convert_from_Kv_to_Kw(u_interior_Kv,v_Kw)
                else:
                    Kw, v_Kw = get_updated_Kw_v_Kw(u_interior_Kv_updated,u_interior_Kv,Kw,v_Kw,principal_cusp_width)
                    update_terms_linear_in_u_to_new_u(coeff_tuples_list,u_interior_Kv_updated,u_interior_Kv,Kw,principal_cusp_width)
                    u_interior_Kv = u_interior_Kv_updated
                    u = Kw.gen()
                has_u_been_improved = True
            #Then continue and try to recognize remaining terms
            have_all_coeffs_lists_been_recognized = recognize_coeff_tuples_list(coeff_tuples_list,Kv,u,v_Kw,coeff_bit_prec) #Note that this function might also mutate coeff_tuples_list

            if have_all_coeffs_lists_been_recognized == True:
                polynomial_ring = PolynomialRing(Kw,"x")
                x = polynomial_ring.gen()
                p3 = Factored_Polynomial(x,coeff_tuples=coeff_tuples_list[0])
                p2 = Factored_Polynomial(x,coeff_tuples=coeff_tuples_list[1])
                pc = Factored_Polynomial(x,coeff_tuples=coeff_tuples_list[2])
                alg_factored_polynomials = (p3, p2, pc)
                return alg_factored_polynomials, v_Kw, u_interior_Kv
    
    if stop_when_coeffs_are_recognized == True:
        raise ArithmeticError("target_bit_prec was not sufficient to recognize coefficients as algebraic numbers!")

    return factored_polynomials

cpdef get_factored_polynomial_starting_values(S, digit_prec, return_cusp_rep_values=False):
    """
    Get first approximation of factored polynomial by computing the hauptmodul to digit_prec digits precision.
    """
    G = S.group()
    cdef Acb_Mat c
    c, M = get_coefficients_haupt_ir_arb_wrap(S,digit_prec,only_principal_expansion=False,return_M=True)
    bit_prec = digits_to_bits(2*digit_prec)
    CBF = ComplexBallField(bit_prec)
    cdef Polynomial_complex_arb x = Polynomial_complex_arb(CBF['x'], is_gen=True)

    p2 = get_p2_haupt(G, x, c, M, bit_prec)
    p3 = get_p3_haupt(G, x, c, M, bit_prec)
    pc, cusp_rep_values = get_pc_haupt(G, x, c, M, bit_prec, return_cusp_rep_values=True)
    factored_polynomials = (p3, p2, pc)

    if return_cusp_rep_values == False:
        return factored_polynomials
    else:
        j_G_hejhal = get_hauptmodul_q_expansion_approx(S,digit_prec,M_0=M,c_vec=c)
        return factored_polynomials, cusp_rep_values, j_G_hejhal

cpdef get_coeff_min_precision(factored_polynomials, int N):
    """
    Returns an estimate of the current coefficient precision (in digits) of factored_polynomials.
    """
    cdef Polynomial_complex_arb f = get_f(factored_polynomials) #This gets also computed in newton_step so for optimization we should re-use it
    cdef Acb_Mat f_x = Acb_Mat(N, 1)
    acb_mat_set_poly(f_x.value, f.__poly)
    cdef int i, digit_prec, real_prec, imag_prec
    cdef int smallest_digit_prec = 2147483647
    RBF = RealBallField(53) #The precision doesn't matter here
    cdef RealBall RB = RBF(0)
    f_x_mcbd = f_x._get_mcbd(53) #It is more convenient to use Sage's class here
    for i in range(N):
        real_prec, imag_prec = get_decimal_digit_prec(f_x_mcbd[i][0].real().mid()), get_decimal_digit_prec(f_x_mcbd[i][0].imag().mid())
        if real_prec < smallest_digit_prec:
            smallest_digit_prec = real_prec
        if imag_prec < smallest_digit_prec:
            smallest_digit_prec = imag_prec
    return smallest_digit_prec

cpdef run_newton(S, starting_digit_prec, target_digit_prec, max_extension_field_degree=None, stop_when_coeffs_are_recognized=True, return_cusp_rep_values=False):
    G = S.group()
    factored_polynomials, cusp_rep_values, j_G_hejhal = get_factored_polynomial_starting_values(S, starting_digit_prec, return_cusp_rep_values=True)
    curr_bit_prec = digits_to_bits(2*starting_digit_prec)
    target_bit_prec = digits_to_bits(target_digit_prec)

    if stop_when_coeffs_are_recognized == True:
        factored_polynomials, v_Kw, u_interior_Kv = newton(factored_polynomials,G,curr_bit_prec,target_bit_prec,stop_when_coeffs_are_recognized,max_extension_field_degree=max_extension_field_degree)
    else:
        factored_polynomials = newton(factored_polynomials,G,curr_bit_prec,target_bit_prec,stop_when_coeffs_are_recognized,max_extension_field_degree=max_extension_field_degree)
    if return_cusp_rep_values == False:
        return factored_polynomials
    else:
        if stop_when_coeffs_are_recognized == True:
            return factored_polynomials, cusp_rep_values, j_G_hejhal, v_Kw, u_interior_Kv
        return factored_polynomials, cusp_rep_values, j_G_hejhal