from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.laurent_series_ring import LaurentSeriesRing

def convert_from_Kv_to_Ku(expression_in_Kv, v_Ku):
    """
    Given an expression in Kv, convert the expression efficiently to Ku by plugging in v(u).
    """
    if expression_in_Kv == 0:
        return 0
    coeffs = list(expression_in_Kv.polynomial())
    res = coeffs[-1]
    for i in range(len(coeffs)-2,-1,-1): #Horner's method
        res = res*v_Ku+coeffs[i]
    return res

def convert_from_Ku_to_Kv(expression_in_Ku, u_interior_Kv, principal_cusp_width):
    """
    Given an expression in Ku, convert the expression to Kv by plugging in u_interior_Kv if possible.
    Otherwise use Sage to convert between the numberfields.
    """
    Ku, Kv = expression_in_Ku.parent(), u_interior_Kv.parent()
    if is_u_minpoly_maximal_degree(Ku.polynomial(),Kv.degree(),principal_cusp_width) == True: #We can efficiently factor expression by plugging in u_interior_Kv 
        coeffs = list(expression_in_Ku)
        res = coeffs[-1]
        for i in range(len(coeffs)-2,-1,-1): #Horner's method
            if coeffs[i] != 0:
                if i%principal_cusp_width != 0:
                    raise ArithmeticError("Unable to factor expression into u and v!")
                res = res*u_interior_Kv+coeffs[i]
    else: #We have to use Sage for this computation... Can we do it more efficiently?
        res = Kv(expression_in_Ku)
    return res

def is_u_minpoly_maximal_degree(u_minpoly, extension_field_degree, principal_cusp_width):
    """
    Detect if u_minpoly is of maximal degree (i.e., of degree extension_field_degree*principal_cusp_width).
    If u_minpoly is of maximal degree, then its monomials are given by [x**(N*principal_cusp_width),x**((N-1)*principal_cusp_width),...]
    which means that some operations can be simplified.
    """
    return u_minpoly.degree() == extension_field_degree*principal_cusp_width

def factor_into_u_v(expression_in_Ku, u_pow, u_interior_Kv, principal_cusp_width):
    """
    Given an expression in Ku that can be factored into (expression_in_Kv)*u**u_pow, perform this factorization.
    The result will be represented as a polynomial over Kv in u.
    """
    Ku, Kv = expression_in_Ku.parent(), u_interior_Kv.parent()
    Ku_gen = Ku.gen()
    expression_shifted = expression_in_Ku/(Ku_gen**u_pow) #This expression can be written as an element in v
    expression_shifted_Kv = convert_from_Ku_to_Kv(expression_shifted,u_interior_Kv,principal_cusp_width)
    P = PolynomialRing(Kv,"u")
    u = P.gen()
    return expression_shifted_Kv*u**u_pow

def factor_q_expansion_into_u_v(q_expansion, u_interior_Kv, principal_cusp_width):
    """
    Given a q_expansion that is represented as a Laurent series, factor coefficients into u and v.
    """
    Kv = u_interior_Kv.parent()
    Pu = PolynomialRing(Kv,"u")
    Pq = LaurentSeriesRing(Pu,q_expansion.variable())
    q = Pq.gen()
    leading_order_exponent = q_expansion.valuation()
    res, u_pow = 0, 0
    for i in range(leading_order_exponent,q_expansion.prec()):
        res += factor_into_u_v(q_expansion[i],u_pow,u_interior_Kv,principal_cusp_width)*q**i
        u_pow += 1
    return res.O(q_expansion.prec())