from sage.rings.polynomial.polynomial_complex_arb import Polynomial_complex_arb
from sage.rings.qqbar import QQbar
from sage.rings.complex_field import ComplexField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from belyi.number_fields import to_QQbar, get_numberfield_and_gen

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

cpdef get_algebraic_poly_coeffs(p, gen, extension_field_degree, principal_cusp_width, estimated_bit_prec=None):
    """
    Given a polynomial p, try to recognize coefficients as algebraic numbers.
    We assume that the generators of the numberfield have already been identified.
    Return False if this does not succeed.
    Otherwise return polynomial over QQbar.
    """
    if estimated_bit_prec == None: #We have not specified the precision so we use the full working precision
        bit_prec = p[0].parent().precision()
    else:
        bit_prec = estimated_bit_prec
    CC = ComplexField(bit_prec)
    algebraic_coeffs = []
    for i in range(p.degree()+1):
        coeff_floating_approx = CC(p[i]) #Because we cannot directly convert acbs to pari
        expression_to_recognize = coeff_floating_approx**principal_cusp_width 
        if expression_to_recognize.is_one() == True:
            recognized_expression = QQbar(1)
        else:
            recognized_expression = to_QQbar(expression_to_recognize,gen,extension_field_degree)

        if recognized_expression == False: #Found an invalid example, therefore precision is insufficient to recognize alg numbers
            return False
        #We need to recognize the correct root. Is there a better way for this?
        potential_algebraic_expressions = recognized_expression.nth_root(principal_cusp_width,all=True)
        diffs = [(potential_algebraic_expression-coeff_floating_approx).abs() for potential_algebraic_expression in potential_algebraic_expressions]
        algebraic_expression = potential_algebraic_expressions[diffs.index(min(diffs))]
        algebraic_coeffs.append(algebraic_expression)
    var_name = p.variable_name()
    polynomial_ring = PolynomialRing(QQbar,var_name)
    p_algebraic = polynomial_ring(algebraic_coeffs)

    return p_algebraic

def get_numberfield_of_poly(p, max_extension_field_degree, principal_cusp_width, estimated_bit_prec=None):
    """
    Try to recognize the numberfield of (one of the coefficients) of p by trying to express the first non-trivial coefficient as an algebraic number.
    Note that we define the numberfield to be the numberfield of c**principal_cusp_width.
    If this succeeds, return the (potentially reduced) numberfield and its generator, otherwise return False.
    """
    if estimated_bit_prec == None: #We have not specified the precision so we use the full working precision
        bit_prec = p[0].parent().precision()
    else:
        bit_prec = estimated_bit_prec
    CC = ComplexField(bit_prec)
    coeff_floating_approx = CC(p[p.degree()-1]) #The second leading coefficient is usually easiest to identify
    expression_to_recognize = coeff_floating_approx**principal_cusp_width
    res = get_numberfield_and_gen(expression_to_recognize, max_extension_field_degree)
    return res

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
    
    def get_algebraic_expressions(self, gen, extension_field_degree, principal_cusp_width, estimated_bit_prec=None):
        """
        Tries to recognize coefficients of factor polynomials as algebraic numbers defined over a numberfield with generator gen.
        If this succeeds (which we only verify empirically here), return instance of Factored_Polynomial over algebraic numbers.
        Otherwise return False.
        """
        if len(self.factors) == 0:
            return self #The empty class is already (somewhat) algebraic
        algebraic_factors = []
        for (p,order) in self.factors:
            p_algebraic = get_algebraic_poly_coeffs(p, gen, extension_field_degree, principal_cusp_width, estimated_bit_prec=estimated_bit_prec)
            if p_algebraic == False:
                return False
            else:
                algebraic_factors.append([p_algebraic,order])

        polygen = algebraic_factors[0][0][0].parent().gen()
        algebraic_factored_polynomial = Factored_Polynomial(polygen,coeff_tuples=algebraic_factors)
        return algebraic_factored_polynomial
    
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