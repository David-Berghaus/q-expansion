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

def get_numberfield_of_coeff(floating_expression_linear_in_u, max_extension_field_degree, principal_cusp_width, estimated_bit_prec=None):
    """
    Try to recognize the numberfield of one of the coefficients by trying to express it as an algebraic number.
    Note that we define the numberfield Kv to be the numberfield of x**principal_cusp_width.
    If this succeeds, return the (potentially reduced) numberfield Kv with specified embedding and the numberfield of w 
    (as well as expressions to convert between Kw and Kv).
    Otherwise return None.
    """
    if estimated_bit_prec == None: #We have not specified the precision so we use the full working precision
        bit_prec = floating_expression_linear_in_u.parent().precision()
    else:
        bit_prec = estimated_bit_prec
    CC = ComplexField(bit_prec)
    coeff_floating_approx = CC(floating_expression_linear_in_u)
    expression_to_recognize = coeff_floating_approx**principal_cusp_width
    Kv = get_numberfield_and_gen(expression_to_recognize,max_extension_field_degree)
    if Kv == None:
        return None
    v = Kv.gen()
    tmp = get_u(coeff_floating_approx,v,principal_cusp_width,Kv.degree())
    if tmp == None:
        return None
    u_interior_Kv, u_embedding = tmp
    Kw = get_Kw(u_interior_Kv,u_embedding,principal_cusp_width,Kv.degree())
    v_Kw = get_v_Kw(Kv,Kw,principal_cusp_width,bit_prec)

    return Kv, Kw, v_Kw, u_interior_Kv

def get_Kw(u_interior_Kv, u_embedding, principal_cusp_width, extension_field_degree):
    """
    Get Kw which we define to be a numberfield that contains all coefficients of the Belyi map.
    We choose Kw to be equal to Kv iff the principal cusp width is equal to one.
    Otherwise we choose Kw to be equivalent to Ku.
    """
    if principal_cusp_width == 1:
        Kv = u_interior_Kv.parent()
        CC = u_embedding.parent()
        Kw = NumberField(Kv.polynomial(),"w",embedding=CC(Kv.gen()))
    else:
        Kw = get_numberfield_of_u(u_interior_Kv,u_embedding,principal_cusp_width,extension_field_degree,"w")
    return Kw

def get_numberfield_of_u(u_interior_Kv, u_embedding, principal_cusp_width, extension_field_degree, gen_str):
    """
    Let u be a principal_cusp_widths root of an expression in Kv.
    Return the NumberField containing u with specified embedding.
    """
    potential_algebraic_expressions = QQbar(u_interior_Kv).nth_root(principal_cusp_width,all=True)
    diffs = [(potential_algebraic_expression-u_embedding).abs() for potential_algebraic_expression in potential_algebraic_expressions]
    u_alg = potential_algebraic_expressions[diffs.index(min(diffs))]
    u_minpoly = get_u_minpoly(u_alg,principal_cusp_width,extension_field_degree,u_embedding.prec())
    Ku = NumberField(u_minpoly,gen_str,embedding=u_embedding) #Embedding with u_alg can apparently become very slow
    return Ku

def get_u_minpoly(u_alg, principal_cusp_width, extension_field_degree, starting_bit_prec):
    """
    Return the minimal polynomial of u by using the LLL algorithm.
    We could of course simply use u.minpoly() which can however become very slow.
    """
    LLL_res = None
    bit_prec = starting_bit_prec
    while LLL_res == None:
        CC = ComplexField(bit_prec)
        LLL_basis = []
        u_approx_pow = CC(u_alg)**principal_cusp_width
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
        u_approx = CC(u_alg)
        factors = p.factor()
        min_factor_root_diffs = [] #Stores the minimal diff of root of factor
        for (factor,multiplicity) in factors:
            factor_roots = factor.roots(ring=QQbar,multiplicities=False)
            diffs = [(factor_root-u_approx).abs() for factor_root in factor_roots]
            min_factor_root_diffs.append(min(diffs))
        factor_index = min_factor_root_diffs.index(min(min_factor_root_diffs)) #Index of the factor for which u is a root
        p = factors[factor_index][0]
    return p

def get_v_Kw(Kv, Kw, principal_cusp_width, starting_bit_prec):
    """
    Express v as an element in Kw. This is in principle possible by using v_Kw = Kw(v) which can however compute for very long.
    Instead we express v in terms of w by using the LLL algorithm.
    """
    v, w = Kv.gen(), Kw.gen()
    if principal_cusp_width == 1: #In this case Kw == Kv (up to different variable names)
        v_Kw = w #In this scenario it is trivial
        return v_Kw
    LLL_res = None
    bit_prec = starting_bit_prec
    while LLL_res == None:
        CC = ComplexField(bit_prec)
        LLL_basis = [CC(v)]
        w_approx_pow = CC(w)**principal_cusp_width
        for i in range(Kv.degree()):
            LLL_basis.append(w_approx_pow**i)
        LLL_res = lindep(LLL_basis,check_if_result_is_invalid=True)
        bit_prec *= 2
        if LLL_res == None:
            print("Unable to recognize v_Kw, maybe optimize bitprec further!")
    v_Kw = 0
    for i in range(1,len(LLL_res)):
        power = (i-1)*principal_cusp_width
        v_Kw += LLL_res[i]*w**power
    v_Kw /= -LLL_res[0]
    return v_Kw

def get_u(floating_expression_linear_in_u, v, principal_cusp_width, extension_field_degree):
    """
    Determine u which is a factor that should cancel out common factors in the coefficients.
    This function returns u_interior_Kv which is an expression in Kv such that u = (u_interior_Kv)**(1/principal_cusp_width)
    and the embedding of u into CC.
    If principal_cusp_width == 1 we choose u to be in QQ (and hence independently of Kv).
    """
    if principal_cusp_width == 1:
        u_interior_Kv = QQ(1)*v**0
        CC = floating_expression_linear_in_u.parent()
        u_embedding = CC(u_interior_Kv)
    else:
        x = floating_expression_linear_in_u**principal_cusp_width
        tmp = get_u_factor(x,v,principal_cusp_width,extension_field_degree)
        if tmp == None:
            return None
        c, recognized_expression = tmp
        u_interior_Kv = recognized_expression/(c**principal_cusp_width) #This corresponds to u(v)**princial_cusp_width, i.e., u(v) = (u_interior_Kv)**(1/princial_cusp_width)
        u_embedding = floating_expression_linear_in_u/c
    return u_interior_Kv, u_embedding

def get_u_factor(x, v, principal_cusp_width, extension_field_degree):
    """
    Given a floating expression which can be written as x = (c*u)^principal_cusp_width, get an expression for c and QQbar(x).
    This is done by factoring out common primes whose power is larger than principal_cusp_width.
    """
    LLL_basis = [x]
    v_approx = x.parent()(v)
    for i in range(extension_field_degree):
        LLL_basis.append(v_approx**i)
    LLL_res = lindep(LLL_basis,check_if_result_is_invalid=True)
    if LLL_res == None:
        return None
    x_alg = 0
    for i in range(1,len(LLL_res)):
        x_alg += LLL_res[i]*v**(i-1)
    x_alg /= -LLL_res[0]
    numerator_primes_pows = [] #Primes and powers
    numerator_primes = [] #Only primes
    for i in LLL_res[1:]:
        if i == 0: #Factorization of zero is not defined
            numerator_primes_pows.append([(0,1)])
            numerator_primes.append([0])
        else:
            i_primes_pows = list(factor(i))
            i_primes = [prime for (prime,_) in i_primes_pows]
            numerator_primes_pows.append(i_primes_pows)
            numerator_primes.append(i_primes)
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
        #Should we also try to improvement the numerator here (for cusp-width one this does usually not seem to work)  
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