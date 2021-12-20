from sage.rings.polynomial.polynomial_complex_arb import Polynomial_complex_arb
from sage.rings.qqbar import QQbar
from sage.rings.rational_field import QQ
from sage.rings.complex_field import ComplexField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.number_field.number_field import NumberField
from sage.arith.misc import factor, prime_factors

from belyi.number_fields import to_K, get_numberfield_and_gen, is_effectively_zero, lindep, convert_from_Kv_to_Ku
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

cpdef get_algebraic_poly_coeffs_nth_power(p, K, principal_cusp_width, estimated_bit_prec=None):
    """
    Given a polynomial p, try to recognize coefficients as algebraic numbers by raising them to the n-th power,
    where n denots the principal_cusp_width.
    We assume that the numberfield and its embedding have already been identified.
    Return None if this does not succeed.
    Otherwise return polynomial over K.
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
            recognized_expression = to_K(expression_to_recognize,K)

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

cpdef get_algebraic_poly_coeffs_u(p, Kv, u, v_Ku, estimated_bit_prec=None):
    """
    Given a polynomial p, try to recognize coefficients as algebraic numbers by dividing them by a power of u.
    We assume that the numberfield and its embedding have already been identified.
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
            recognized_expression = Kv(1)
        elif is_effectively_zero(expression_to_recognize,bits_to_digits(bit_prec)) == True:
            recognized_expression = Kv(0)
        else:
            recognized_expression = to_K(expression_to_recognize,Kv)

        if recognized_expression == None: #Found an invalid example, therefore precision is insufficient to recognize alg numbers
            return None

        recognized_expression_Ku = convert_from_Kv_to_Ku(recognized_expression, v_Ku)
        algebraic_expression = recognized_expression_Ku*(u**u_pow)
        algebraic_coeffs.append(algebraic_expression)
    var_name = p.variable_name()
    Ku = v_Ku.parent()
    polynomial_ring = PolynomialRing(Ku,var_name)
    p_algebraic = polynomial_ring(algebraic_coeffs)

    return p_algebraic

def get_numberfield_of_coeff(x, max_extension_field_degree, principal_cusp_width, estimated_bit_prec=None):
    """
    Try to recognize the numberfield of one of the coefficients by trying to express it as an algebraic number.
    Note that we define the numberfield Kv to be the numberfield of x**principal_cusp_width.
    If this succeeds, return the (potentially reduced) numberfield Kv with specified embedding and the numberfield of u 
    (as well as expressions to convert between Ku and Kv).
    Otherwise return None.
    """
    if estimated_bit_prec == None: #We have not specified the precision so we use the full working precision
        bit_prec = x.parent().precision()
    else:
        bit_prec = estimated_bit_prec
    CC = ComplexField(bit_prec)
    coeff_floating_approx = CC(x)
    expression_to_recognize = coeff_floating_approx**principal_cusp_width
    Kv = get_numberfield_and_gen(expression_to_recognize,max_extension_field_degree)
    if Kv == None:
        return None
    tmp = get_u_factor(expression_to_recognize,Kv.gen(),principal_cusp_width,Kv.degree())
    if tmp == None: #Although this should never happen because the same computation was done in get_numberfield_and_gen...
        return None
    c, recognized_expression = tmp
    u_interior_Kv = recognized_expression/(c**principal_cusp_width) #This corresponds to u(v)**princial_cusp_width, i.e., u(v) = (u_interior_Kv)**(1/princial_cusp_width)
    Ku = get_numberfield_of_u(c,u_interior_Kv,coeff_floating_approx,principal_cusp_width,Kv.degree())
    v_Ku = get_v_Ku(Kv,Ku,principal_cusp_width,bit_prec)

    return Kv, Ku, v_Ku, u_interior_Kv

def get_numberfield_of_u(c, u_interior_Kv, coeff_floating_approx, principal_cusp_width, extension_field_degree):
    """
    Let u be a principal_cusp_widths root of an expression in Kv.
    Return the NumberField containing u with specified embedding.
    """
    u_floating_approx = coeff_floating_approx/c
    potential_algebraic_expressions = QQbar(u_interior_Kv).nth_root(principal_cusp_width,all=True)
    diffs = [(potential_algebraic_expression-u_floating_approx).abs() for potential_algebraic_expression in potential_algebraic_expressions]
    u_alg = potential_algebraic_expressions[diffs.index(min(diffs))]
    u_minpoly = get_u_minpoly(u_alg,principal_cusp_width,extension_field_degree,coeff_floating_approx.prec())
    Ku = NumberField(u_minpoly,"u",embedding=u_floating_approx) #Embedding with u_alg can apparently become very slow
    return Ku

def get_u_minpoly(u, principal_cusp_width, extension_field_degree, starting_bit_prec):
    """
    Return the minimal polynomial of u by using the LLL algorithm.
    We could of course simply use u.minpoly() which can however become very slow.
    """
    LLL_res = None
    bit_prec = starting_bit_prec
    while LLL_res == None:
        CC = ComplexField(bit_prec)
        LLL_basis = []
        u_approx_pow = CC(u)**principal_cusp_width
        for i in range(extension_field_degree+1):
            LLL_basis.append(u_approx_pow**i)
        LLL_res = lindep(LLL_basis,check_if_result_is_invalid=True)
        bit_prec *= 2
        if LLL_res == None:
            print("Unable to recognize minpoly of u, maybe optimize bitprec further!")
    P = PolynomialRing(QQ,"x")
    p = P(LLL_res)
    x = P.gen()
    p = p.subs({x:x**principal_cusp_width})
    if p.is_irreducible() == False:
        CC = ComplexField(bit_prec)
        u_approx = CC(u)
        factors = p.factor()
        min_factor_root_diffs = [] #Stores the minimal diff of root of factor
        for (factor,multiplicity) in factors:
            factor_roots = factor.roots(ring=QQbar,multiplicities=False)
            diffs = [(factor_root-u_approx).abs() for factor_root in factor_roots]
            min_factor_root_diffs.append(min(diffs))
        factor_index = min_factor_root_diffs.index(min(min_factor_root_diffs)) #Index of the factor for which u is a root
        p = factors[factor_index][0]
    return p

def get_v_Ku(Kv, Ku, principal_cusp_width, starting_bit_prec):
    """
    Express v as an element in Ku. This is in principle possible by using v_Ku = Ku(v) which can however compute for very long.
    Instead we express v in terms of u by using the LLL algorithm.
    """
    v, u = Kv.gen(), Ku.gen()
    LLL_res = None
    bit_prec = starting_bit_prec
    while LLL_res == None:
        CC = ComplexField(bit_prec)
        LLL_basis = [CC(v)]
        u_approx_pow = CC(u)**principal_cusp_width
        for i in range(Kv.degree()):
            LLL_basis.append(u_approx_pow**i)
        LLL_res = lindep(LLL_basis,check_if_result_is_invalid=True)
        bit_prec *= 2
        if LLL_res == None:
            print("Unable to recognize v_Ku, maybe optimize bitprec further!")
    v_Ku = 0
    for i in range(1,len(LLL_res)):
        power = (i-1)*principal_cusp_width
        v_Ku += LLL_res[i]*u**power
    v_Ku /= -LLL_res[0]
    return v_Ku

def get_u_factor(x, gen, principal_cusp_width, extension_field_degree):
    """
    Given an expression x which can be written as x = (c*u)^principal_cusp_width, get an expression for c and QQbar(x).
    This is done by factoring out common primes whose power is larger than principal_cusp_width.
    """
    LLL_basis = [x]
    gen_approx = x.parent()(gen)
    for i in range(extension_field_degree):
        LLL_basis.append(gen_approx**i)
    LLL_res = lindep(LLL_basis,check_if_result_is_invalid=True)
    if LLL_res == None:
        return None
    x_alg = 0
    for i in range(1,len(LLL_res)):
        x_alg += LLL_res[i]*gen**(i-1)
    x_alg /= -LLL_res[0]
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
    return c, x_alg

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
    
    def get_algebraic_expressions(self, K, principal_cusp_width=None, u=None, v_Ku=None, estimated_bit_prec=None):
        """
        Tries to recognize coefficients of factor polynomials as algebraic numbers defined over K.
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
                p_algebraic = get_algebraic_poly_coeffs_nth_power(p, K, principal_cusp_width, estimated_bit_prec=estimated_bit_prec)
            else:
                p_algebraic = get_algebraic_poly_coeffs_u(p, K, u, v_Ku, estimated_bit_prec=estimated_bit_prec)
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