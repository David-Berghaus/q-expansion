from sage.rings.polynomial.polynomial_complex_arb import Polynomial_complex_arb
from sage.rings.qqbar import QQbar
from sage.rings.rational_field import QQ
from sage.rings.complex_field import ComplexField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.arith.misc import factor, prime_factors

from belyi.number_fields import to_QQbar, get_numberfield_and_gen, is_effectively_zero, lindep
from point_matching.point_matching_arb_wrap import bits_to_digits

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

cpdef get_algebraic_poly_coeffs_nth_power(p, gen, extension_field_degree, principal_cusp_width, estimated_bit_prec=None):
    """
    Given a polynomial p, try to recognize coefficients as algebraic numbers by raising them to the n-th power,
    where n denots the principal_cusp_width.
    We assume that the generators of the numberfield have already been identified.
    Return None if this does not succeed.
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
        elif is_effectively_zero(expression_to_recognize,bits_to_digits(bit_prec)) == True:
            recognized_expression = QQbar(0)
        else:
            recognized_expression = to_QQbar(expression_to_recognize,gen,extension_field_degree)

        if recognized_expression == None: #Found an invalid example, therefore precision is insufficient to recognize alg numbers
            return None
        #We need to recognize the correct root. Is there a better way for this?
        potential_algebraic_expressions = recognized_expression.nth_root(principal_cusp_width,all=True)
        diffs = [(potential_algebraic_expression-coeff_floating_approx).abs() for potential_algebraic_expression in potential_algebraic_expressions]
        algebraic_expression = potential_algebraic_expressions[diffs.index(min(diffs))]
        algebraic_coeffs.append(algebraic_expression)
    var_name = p.variable_name()
    polynomial_ring = PolynomialRing(QQbar,var_name)
    p_algebraic = polynomial_ring(algebraic_coeffs)

    return p_algebraic

cpdef get_algebraic_poly_coeffs_u(p, gen, extension_field_degree, u, estimated_bit_prec=None):
    """
    Given a polynomial p, try to recognize coefficients as algebraic numbers by dividing them by a power of u.
    We assume that the generators of the numberfield have already been identified.
    Return None if this does not succeed.
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
        u_pow = p.degree()-i #Exponent of u for the given coefficient
        expression_to_recognize = coeff_floating_approx/(u**u_pow)
        if expression_to_recognize.is_one() == True:
            recognized_expression = QQbar(1)
        elif is_effectively_zero(expression_to_recognize,bits_to_digits(bit_prec)) == True:
            recognized_expression = QQbar(0)
        else:
            recognized_expression = to_QQbar(expression_to_recognize,gen,extension_field_degree)

        if recognized_expression == None: #Found an invalid example, therefore precision is insufficient to recognize alg numbers
            return None
        algebraic_expression = recognized_expression*(u**u_pow)
        algebraic_coeffs.append(algebraic_expression)
    var_name = p.variable_name()
    polynomial_ring = PolynomialRing(QQbar,var_name)
    p_algebraic = polynomial_ring(algebraic_coeffs)

    return p_algebraic

def get_numberfield_of_coeff(x, max_extension_field_degree, principal_cusp_width, estimated_bit_prec=None):
    """
    Try to recognize the numberfield of one of the coefficients by trying to express it as an algebraic number.
    Note that we define the numberfield to be the numberfield of x**principal_cusp_width.
    If this succeeds, return the (potentially reduced) numberfield, its generator and an algebraic expression for u.
    Otherwise return None.
    """
    if estimated_bit_prec == None: #We have not specified the precision so we use the full working precision
        bit_prec = x.parent().precision()
    else:
        bit_prec = estimated_bit_prec
    CC = ComplexField(bit_prec)
    coeff_floating_approx = CC(x)
    expression_to_recognize = coeff_floating_approx**principal_cusp_width
    res = get_numberfield_and_gen(expression_to_recognize,max_extension_field_degree)
    if res == None:
        return None
    nf, gen = res
    c = get_u_factor(expression_to_recognize,gen,principal_cusp_width,nf.degree())
    recognized_expression = to_QQbar(expression_to_recognize,gen,nf.degree())
    if recognized_expression == None: #Although this should never happen because the same computation was done in get_numberfield_and_gen...
        return None
    potential_algebraic_expressions = recognized_expression.nth_root(principal_cusp_width,all=True)
    diffs = [(potential_algebraic_expression-coeff_floating_approx).abs() for potential_algebraic_expression in potential_algebraic_expressions]
    algebraic_expression = potential_algebraic_expressions[diffs.index(min(diffs))]
    u = algebraic_expression/c
    return nf, gen, u

def get_u_factor(x, gen, principal_cusp_width, extension_field_degree):
    """
    Given an expression x which can be written as x = (c*u)^principal_cusp_width, get an expression for c.
    This is done by factoring out common primes whose power is larger than principal_cusp_width.
    """
    LLL_basis = [x]
    gen_approx = x.parent()(gen)
    for i in range(extension_field_degree):
        LLL_basis.append(gen_approx**i)
    LLL_res = lindep(LLL_basis,check_if_result_is_invalid=True)
    if LLL_res == None:
        return None
    numerator_primes_pows = [list(factor(i)) for i in LLL_res[1:]] #Primes and powers
    numerator_primes = [prime_factors(i) for i in LLL_res[1:]] #Only primes
    print("numerator_primes_pows: ", numerator_primes_pows)
    print("denominator: ", factor(LLL_res[0]))
    numerator_factors = [] #Primes and powers of the numerator that can be factored out
    #Now look for common prime factors which have a power that is larger than principal_cusp_width
    for (prime,power) in numerator_primes_pows[0]:
        if power < principal_cusp_width: #This prime is cannot be factored out
            continue
        powers = [power] #Powers of the given prime for all coefficients
        for i in range(1,len(numerator_primes_pows)):
            if prime in numerator_primes[i]:
                powers.append(numerator_primes_pows[i][numerator_primes[i].index(prime)][1]) #Ugly!
            else:
                break
        if len(powers) != len(numerator_primes_pows): #The above loop breaked early which means that prime is not a factor in all coeffs
            continue
        min_pow = min(powers)
        if min_pow < principal_cusp_width:
            continue
        numerator_factors.append( (prime,min_pow-min_pow%principal_cusp_width) )
    print("numerator_factors: ", numerator_factors)
    c = QQ(1)
    for (prime,power) in numerator_factors:
        c *= prime**(power//principal_cusp_width)
    return c

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
    
    def get_algebraic_expressions(self, gen, extension_field_degree, principal_cusp_width=None, u=None, estimated_bit_prec=None):
        """
        Tries to recognize coefficients of factor polynomials as algebraic numbers defined over a numberfield with generator gen.
        If principal_cusp_width != None, the coefficients are raised to the n-th power, where n denots the principal_cusp_width.
        If u != None, this is done by dividing the coefficients by a power of u.
        If this succeeds (which we only verify empirically here), return instance of Factored_Polynomial over algebraic numbers.
        Otherwise return None.
        """
        if principal_cusp_width == None and u == None:
            raise ArithmeticError("Please provide either princial_cusp_width or u!")
        if len(self.factors) == 0:
            return self #The empty class is already (somewhat) algebraic
        algebraic_factors = []
        for (p,order) in self.factors:
            if u == None:
                p_algebraic = get_algebraic_poly_coeffs_nth_power(p, gen, extension_field_degree, principal_cusp_width, estimated_bit_prec=estimated_bit_prec)
            else:
                p_algebraic = get_algebraic_poly_coeffs_u(p, gen, extension_field_degree, u, estimated_bit_prec=estimated_bit_prec)
            if p_algebraic == None:
                return None
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