from sage.rings.polynomial.polynomial_complex_arb import Polynomial_complex_arb
from sage.rings.qqbar import QQbar
from sage.rings.complex_field import ComplexField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from belyi.number_fields import to_QQbar

cpdef construct_poly_from_root_tuple(x, root_tuple):
    p = 1
    (roots, order) = root_tuple
    for root in roots:
        p *= (x-root) #For acb_poly we could use native implementations here for more performance...
    return [p,order]

cpdef construct_poly_from_coeff_tuple(x, coeff_tuple):
    (coeffs, order) = coeff_tuple
    if isinstance(x,Polynomial_complex_arb) == True: #We specifically work with acb_poly instead of polynomials over the acb-ring for performance and C-access
        p = Polynomial_complex_arb(x.parent(), coeffs)
    else:
        polynomial_ring = PolynomialRing(x.parent(),"x")
        p = polynomial_ring(coeffs)
    return [p,order]

cpdef get_algebraic_poly_coeffs(p, max_extension_field_degree):
    """
    Given a polynomial p, try to recognize coefficients as algebraic numbers. 
    Return False if this does not succeed.
    Otherwise return polynomial over QQbar.
    """
    bit_prec = p[0].parent().precision()
    CC = ComplexField(bit_prec)
    algebraic_coeffs = []
    for i in range(p.degree()+1):
        algebraic_expression = to_QQbar(CC(p[i]),max_extension_field_degree+1)
        numberfield = algebraic_expression.minpoly()
        if numberfield[max_extension_field_degree+1] != 0: #Found an invalid example, therefore precision is insufficient to recognize alg numbers
            return False
        else:
            algebraic_coeffs.append(algebraic_expression)
    var_name = p.variable_name()
    polynomial_ring = PolynomialRing(QQbar,var_name)

    return polynomial_ring(algebraic_coeffs)

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
        res = 1
        for (factor, order) in factors:
            res *= factor**order
        return res

    def derivative(self, poly_index, coeff_index):
        """
        Returns a polynomial that corresponds to the derivative with respect to the 'coeff_index'th coefficient of polynomial at 'poly_index'.
        """
        factors = self.factors
        outer_derivative = 1
        for i in range(len(factors)):
            p, order = factors[i]
            if i == poly_index:
                order -= 1
            outer_derivative *= p**order
        inner_derivative = (factors[poly_index][1])*(self.polygen)**coeff_index
        derivative = inner_derivative*outer_derivative #this multiplication by a monomial can certainly be optimized
        return derivative
    
    def get_algebraic_expressions(self, max_extension_field_degree, reduce_numberfields=False):
        """
        Tries to recognize coefficients of factor polynomials as algebraic numbers.
        If this succeeds (which we only verify empirically here), return instance of Factored_Polynomial over algebraic numbers.
        Otherwise return False.
        """
        algebraic_factors = []
        for (p,order) in self.factors:
            p_algebraic = get_algebraic_poly_coeffs(p, max_extension_field_degree)
            if p_algebraic == False:
                return False
            else:
                algebraic_factors.append([p_algebraic,order])
        
        if reduce_numberfields == True: #For spurious numberfields this might take very long so it is important to call this last when it seems likely that the numberfields are correct
            raise NotImplementedError("This functionality has not been added yet!")

        polygen = algebraic_factors[0][0][0].parent().gen()
        algebraic_factored_polynomial = Factored_Polynomial(polygen,coeff_tuples=algebraic_factors)
        return algebraic_factored_polynomial