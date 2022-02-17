from sage.rings.polynomial.polynomial_complex_arb import Polynomial_complex_arb
from sage.rings.qqbar import QQbar
from sage.rings.rational_field import QQ
from sage.rings.complex_field import ComplexField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.number_field.number_field import NumberField, NumberField_generic
from sage.arith.misc import factor, prime_factors

from belyi.number_fields import to_K, get_numberfield_and_gen, is_effectively_zero, lindep
from belyi.expression_in_u_and_v import convert_from_Kv_to_Kw, factor_into_u_v
from point_matching.point_matching_arb_wrap import bits_to_digits

cpdef construct_poly_from_root_tuple(x, root_tuple):
    p = 1
    (roots, order) = root_tuple
    for root in roots:
        p *= (x-root) #For acb_poly we could use native implementations here for more performance...
    return [p,order]

cpdef construct_poly_from_coeff_tuple(x, coeff_tuple):
    (coeffs, order) = coeff_tuple
    if isinstance(x,Polynomial_complex_arb) == True: #We specifically work with acb_poly instead of polynomials over the acb-ring for C-access
        p = Polynomial_complex_arb(x.parent(), coeffs)
    else:
        polynomial_ring = x.parent()
        p = polynomial_ring(coeffs)
    return [p,order]

cpdef recognize_coeffs_using_u(coeffs, Kv, u, v_Kw, estimated_bit_prec=None, max_u_power=None):
    """
    Given a list of coefficients, try to recognize coefficients (those that are CBFs) as algebraic numbers by dividing them by a power of u.
    We assume that the numberfield and its embedding have already been identified.
    This function replaces all CBF-coefficients that have been recognized with its corresponding numberfield expressions.
    If all coeffs have been successfully recognized, the function also returns True, otherwise it returns False.
    If "max_u_power" is specified, only try to recognize coefficients with up to (and including) max_u_power.
    """
    if estimated_bit_prec == None: #We have not specified the precision so we use the full working precision
        bit_prec = coeffs[0].parent().precision()
    else:
        bit_prec = estimated_bit_prec
    CC = ComplexField(bit_prec)
    algebraic_coeffs = []
    p_degree = len(coeffs)-1
    if max_u_power == None:
        i_limit = -1
    else:
        i_limit = len(coeffs)-2-max_u_power
    for i in range(len(coeffs)-1,i_limit,-1): #Loop backwards through coefficients to get increasing powers of u
        if coeffs[i] == 1 or isinstance(coeffs[i].parent(),NumberField_generic) == False: #This coeff has not been recognized yet
            coeff_floating_approx = CC(coeffs[i]) #Because we cannot directly convert acbs to pari
            u_pow = p_degree-i #Exponent of u for the given coefficient
            expression_to_recognize = coeff_floating_approx/(CC(u)**u_pow)
            if expression_to_recognize.is_one() == True:
                recognized_expression = Kv(1)
            elif is_effectively_zero(expression_to_recognize,bits_to_digits(bit_prec)) == True:
                recognized_expression = Kv(0)
            else:
                recognized_expression = to_K(expression_to_recognize,Kv)

            if recognized_expression == None: #Found an invalid example, therefore precision is insufficient to recognize alg numbers
                have_all_coeffs_been_recognized = False
                return coeffs, have_all_coeffs_been_recognized #Return what have already been able to recognize so far

            recognized_expression_Kw = convert_from_Kv_to_Kw(recognized_expression, v_Kw)
            algebraic_expression = recognized_expression_Kw*(u**u_pow)
            coeffs[i] = algebraic_expression #Update list with recognized algebraic expression
    have_all_coeffs_been_recognized = True
    return coeffs, have_all_coeffs_been_recognized

def get_improved_choice_of_u_interior_Kv(coeff_tuples_list, u_interior_Kv, principal_cusp_width):
    """
    Given a list of coefficients for which all coefficients that are linear in u have been recognized,
    use these expressions to try to find an improved choice of u and mutate coeff_tuples_list with the updated choice of u.
    """
    linear_u_expressions = [] #list of all algebraic expressions that are linear in u
    for coeff_tuples in coeff_tuples_list:
        for (coeffs,_) in coeff_tuples:
            if len(coeffs) > 1:
                linear_u_expressions.append(coeffs[-2])
    if len(linear_u_expressions) < 2:
        return u_interior_Kv #Impossible to improve expressions
    if principal_cusp_width == 1:
        largest_denominator = max([factor.denominator() for coeff in linear_u_expressions for factor in list(coeff)])
        u_interior_Kv /= largest_denominator
    else:
        u_expressions_with_u_factored_out = [factor_into_u_v(linear_u_expression,1,u_interior_Kv,principal_cusp_width)[1] for linear_u_expression in linear_u_expressions]
        largest_denominator = max([factor.denominator() for coeff in u_expressions_with_u_factored_out for factor in list(coeff)])
        u_interior_Kv /= largest_denominator**principal_cusp_width
        #Should we also try to improve the numerator here? (for cusp-width one this does usually not seem to work)
        #Note also that for some single-cusp groups it seems like we also need to consider higher order u terms
        #in order to get the best choice of the denominator (which we have not implemented yet)
    return u_interior_Kv

def get_updated_Kw_v_Kw(u_interior_Kv_updated, u_interior_Kv_old, Kw, v_Kw, principal_cusp_width):
    """
    Suppose we have changed u to the new form u_interior_Kv_updated = update_factor*u_interior_Kv where update_factor is in QQ.
    This function returns the corresponding updates to Kw and v_Kw.
    """
    c = QQ(u_interior_Kv_updated/u_interior_Kv_old).nth_root(principal_cusp_width) #u_updated = c*u
    Kw_old_poly_coeffs = list(Kw.polynomial())
    Kw_new_coeffs = []
    for (i,coeff_old) in enumerate(Kw_old_poly_coeffs):
        Kw_new_coeffs.append(coeff_old*c**(len(Kw_old_poly_coeffs)-i-1))
    CC = ComplexField(1024) #Because embedding with QQbar can be very slow
    Kw_updated = NumberField(Kw.polynomial().parent()(Kw_new_coeffs),"w",embedding=CC(c*Kw.gen()))
    w_updated = Kw_updated.gen()
    v_Kw_updated = 0
    for (i,coeff_old) in enumerate(list(v_Kw)):
        v_Kw_updated += (coeff_old/(c**i))*w_updated**i
    return Kw_updated, v_Kw_updated

def update_terms_linear_in_u_to_new_u(coeff_tuples_list, u_interior_Kv_updated, u_interior_Kv_old, Kw_updated, principal_cusp_width):
    """
    Update all terms of coeff_tuples_list that are linear in u to the new choice of u (inplace).
    """
    c = QQ(u_interior_Kv_updated/u_interior_Kv_old).nth_root(principal_cusp_width) #u_updated = c*u
    w_updated = Kw_updated.gen()
    for coeff_tuples in coeff_tuples_list:
        for (coeffs,_) in coeff_tuples:
            if len(coeffs) > 1:
                coeff_updated = 0
                for (i,factor_old) in enumerate(list(coeffs[-2])):
                    coeff_updated += (factor_old/(c**i))*w_updated**i
                coeffs[-2] = coeff_updated

def get_factored_polynomial_in_u_v(factored_polynomial_in_Ku, u_interior_Kv, principal_cusp_width):
    """
    Given a factored polynomial that is defined over Ku, transform coefficients into expressions in u and v
    (which are internally expressed as polynomials in Kv over u).
    Although one can perform some arithmetic with these expressions, we mainly use them to print the results in a more convenient form.
    """
    coeff_tuples = []
    polygen = None
    for (p,order) in factored_polynomial_in_Ku.factors:
        p_deg = p.degree()
        coeffs_u_v = []
        for i in range(0,p_deg+1):
            u_pow = p_deg-i
            coeff_u_v = factor_into_u_v(p[i],u_pow,u_interior_Kv,principal_cusp_width)
            coeffs_u_v.append(coeff_u_v)
            if polygen == None and u_pow != 0:
                polynomial_ring = PolynomialRing(coeff_u_v.parent(),p.parent().variable_name())
                polygen = polynomial_ring.gen()
        coeff_tuples.append( (coeffs_u_v,order) )
    return Factored_Polynomial(polygen,coeff_tuples=coeff_tuples)

class Factored_Polynomial():
    """
    Class for working with polynomials that are factored like:
    p = p_1^(n_1)*...*p_N^(n_N)
    """
    def __init__(self, polygen, root_tuples=None, coeff_tuples=None):
        """
        root_tuples contains a list of tuples. 
        Each tuple is given by a list of ComplexBalls reflecting the roots as well as the order of each of the roots.
        Alternatively we can construct the polynomials directly from their coefficients.
        """
        self.polygen = polygen
        factors = []
        if (root_tuples == None and coeff_tuples == None) or (root_tuples != None and coeff_tuples != None):
            raise ArithmeticError("Invalid construction. Construct either through root_tuples or coeff_tuples!")
        if root_tuples != None: #Construct polynomial from root tuples. This usually happens during the first iteration where the haupt-values are passed
            for root_tuple in root_tuples:
                factors.append(construct_poly_from_root_tuple(polygen,root_tuple))
        if coeff_tuples != None:
            for coeff_tuple in coeff_tuples:
                factors.append(construct_poly_from_coeff_tuple(polygen,coeff_tuple))
        self.factors = factors
    
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        res = ""
        for (p,order) in self.factors:
            res += "("
            res += p.__str__()
            res += ")"
            if order != 1:
                res += "^" + str(order)
        return res

    def construct(self):
        """
        Explicitly constructs self as a polynomial.
        """
        factors = self.factors
        res = self.polygen.parent().one()
        for (factor, order) in factors:
            res *= factor**order
        return res

    def derivative(self, poly_index, coeff_index):
        """
        Returns a polynomial that corresponds to the derivative with respect to the 'coeff_index'th coefficient of polynomial at 'poly_index'.
        """
        factors = self.factors
        outer_derivative = self.polygen.parent().one()
        for i in range(len(factors)):
            p, order = factors[i]
            if i == poly_index:
                order -= 1
            outer_derivative *= p**order
        inner_derivative = (factors[poly_index][1])*(self.polygen)**coeff_index
        derivative = inner_derivative*outer_derivative #this multiplication by a monomial can certainly be optimized
        return derivative
    
    def get_smallest_degree_poly(self):
        """
        Return the polynomial p_i that has the smallest degree.
        This routine is useful because the coefficients of the smallest degree polynomial are usually the easiest to recognize.
        """
        factors = self.factors
        if len(factors) == 0:
            return None
        p_smallest_deg = factors[0][0]
        for i in range(1,len(factors)):
            p,_ = factors[i]
            if p.degree() < p_smallest_deg.degree():
                p_smallest_deg = p
            if p_smallest_deg.degree() == 1: #Cannot get smaller
                break
        return p_smallest_deg