from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.laurent_series_ring import LaurentSeriesRing

def convert_from_Kv_to_Kw(expression_in_Kv, v_Kw):
    """
    Given an expression in Kv, convert the expression efficiently to Kw by plugging in v(w).
    """
    if expression_in_Kv == 0:
        return 0
    coeffs = list(expression_in_Kv.polynomial())
    res = coeffs[-1]
    for i in range(len(coeffs)-2,-1,-1): #Horner's method
        res = res*v_Kw+coeffs[i]
    return res

def convert_from_Kw_to_Kv(expression_in_Kw, u_interior_Kv, principal_cusp_width):
    """
    Given an expression in Kw, convert the expression to Kv.
    For cusp_width one, this is done by replacing w with v.
    For the general case, this is done by plugging in u_interior_Kv if possible.
    Otherwise use Sage to convert between the numberfields.
    """
    Kw, Kv = expression_in_Kw.parent(), u_interior_Kv.parent()
    if principal_cusp_width == 1: #In this case Kw == Kv (up to different variable names)
        v = Kv.gen()
        Kw_expr_poly = expression_in_Kw.polynomial()
        res = 0
        for i in range(Kw_expr_poly.degree()+1):
            res += Kw_expr_poly[i]*v**i
    else:
        if is_u_minpoly_maximal_degree(Kw.polynomial(),Kv.degree(),principal_cusp_width) == True: #We can efficiently factor expression by plugging in u_interior_Kv 
            coeffs = list(expression_in_Kw)
            res = coeffs[-1]
            for i in range(len(coeffs)-2,-1,-1): #Horner's method
                if coeffs[i] != 0:
                    if i%principal_cusp_width != 0:
                        raise ArithmeticError("Unable to factor expression into u and v!")
                    res = res*u_interior_Kv+coeffs[i]
        else: #We have to use Sage for this computation... Can we do it more efficiently?
            res = Kv(expression_in_Kw)
    return res

def is_u_minpoly_maximal_degree(u_minpoly, extension_field_degree, principal_cusp_width):
    """
    Detect if u_minpoly is of maximal degree (i.e., of degree extension_field_degree*principal_cusp_width).
    If u_minpoly is of maximal degree, then its monomials are given by [x**(N*principal_cusp_width),x**((N-1)*principal_cusp_width),...]
    which means that some operations can be simplified.
    """
    return u_minpoly.degree() == extension_field_degree*principal_cusp_width

def factor_into_u_v(expression_in_Kw, u_pow, u_interior_Kv, principal_cusp_width):
    """
    Given an expression in Kw that can be factored into (expression_in_Kv)*u**u_pow, perform this factorization.
    The result will be represented as a polynomial over Kv in u.
    """
    Kw, Kv = expression_in_Kw.parent(), u_interior_Kv.parent()
    w = Kw.gen()
    if principal_cusp_width == 1: #In this case Kw == Kv (up to different variable names)
        expression_shifted = expression_in_Kw/(u_interior_Kv**u_pow) #This expression can be written as an element in v
    else:
        expression_shifted = expression_in_Kw/(w**u_pow) #This expression can be written as an element in v
    expression_shifted_Kv = convert_from_Kw_to_Kv(expression_shifted,u_interior_Kv,principal_cusp_width)
    P = PolynomialRing(Kv,"u")
    u = P.gen()
    return expression_shifted_Kv*u**u_pow

def factor_q_expansion_into_u_v(q_expansion, u_interior_Kv, principal_cusp_width, trunc_order):
    """
    Given a q_expansion that is represented as a Laurent series, factor coefficients into u and v.
    """
    Kv = u_interior_Kv.parent()
    Pu = PolynomialRing(Kv,"u")
    Pq = LaurentSeriesRing(Pu,q_expansion.variable())
    q = Pq.gen()
    leading_order_exponent = q_expansion.valuation()
    res, u_pow = 0, 0
    if trunc_order == None:
        trunc_order = q_expansion.prec()
    for i in range(leading_order_exponent,trunc_order):
        expression_in_Ku = q_expansion[i]
        if expression_in_Ku.polynomial().degree() == 0: #This expression can be defined over QQ (and hence independently of u)
            res += factor_into_u_v(expression_in_Ku,0,u_interior_Kv,principal_cusp_width)*q**i
        else:
            res += factor_into_u_v(expression_in_Ku,u_pow,u_interior_Kv,principal_cusp_width)*q**i
        u_pow += 1
    return res.O(trunc_order)