from copy import copy

from sage.rings.complex_arb import ComplexBallField
from sage.rings.complex_mpc import MPComplexField
from sage.rings.real_mpfr import RealField
from sage.rings.qqbar import QQbar
from sage.rings.rational_field import QQ
from sage.rings.complex_field import ComplexField
from sage.interfaces.gp import gp
from sage.interfaces.gp import pari
from sage.symbolic.constants import pi
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.number_field.number_field import NumberField
from sage.arith.misc import factor, prime_factors

from point_matching.point_matching_arb_wrap import digits_to_bits, bits_to_digits

cpdef get_decimal_digit_prec(epsilon):
    """
    Given a real number epsilon, return the decimal digit precision.
    """
    if epsilon == 0:
        return 2147483647
    else:
        RF = RealField(53)
        eps_approx = RF(epsilon).abs() #There is no reason to perform this computation at full precision
        return -int(eps_approx.log10())

def is_effectively_zero(x, estimated_digit_prec):
    """
    Given a number x, test if x is effectively zero (up to the estimated precision).
    """
    return get_decimal_digit_prec(x.abs()) > estimated_digit_prec

cpdef to_K(x, K):
    """
    Given a floating-point number x, try to express x as an element in K using LLL.
    If this does not succeed, return None.
    """
    gen, extension_field_degree = K.gen(), K.degree()
    gen_approx = x.parent(gen)
    LLL_basis = [x]
    for i in range(extension_field_degree):
        LLL_basis.append(gen_approx**i)
    LLL_res = lindep(LLL_basis,check_if_result_is_invalid=True)
    if LLL_res == None: #Result is invalid
        return None
    
    x_alg = 0
    for i in range(extension_field_degree):
        x_alg += LLL_res[i+1]*gen**i
    x_alg /= -LLL_res[0]
    return x_alg

def get_numberfield_and_gen(x, max_extension_field_degree, reduce_numberfield=True):
    """
    Try to express x in a numberfield over QQbar with specified max_extension_field_degree.
    This function returns a Numberfield with specified embedding or False.
    If reduce_numberfield == True then try to reduce the numberfield using pari.polredabs() and find the new generator.
    """
    LLL_basis = [x.parent().one(),x]
    LLL_res = lindep(LLL_basis,check_if_result_is_invalid=True)
    if LLL_res == None: #Current basis is not large enough to find solution
        for i in range(2,max_extension_field_degree+1): #Look for larger degree extension fields
            LLL_basis.append(x**i) #Increase basis size
            LLL_res = lindep(LLL_basis,check_if_result_is_invalid=True)
            if LLL_res != None: #We found a solution that seems to be valid
                break
        if LLL_res == None: #Found no solution even for the maximal basis size
            return None
    
    P = PolynomialRing(ZZ,"x")
    p = P(LLL_res)
    if reduce_numberfield == False:
        #Now we need to recognize to which root our expression corresponds. Is there a better way for this?
        roots = p.roots(ring=QQbar,multiplicities=False)
        diffs = [(root-x).abs() for root in roots]
        root_index = diffs.index(min(diffs))
        gen = roots[root_index]
        K = NumberField(p,"v",embedding=gen)
        return K
    else:
        p_red = P(pari.polredabs(p))
        print("Identified numberfield: ", p_red)
        if p_red.degree() == 1:
            K = NumberField(p_red-1,"v")
            return K
        else:
            #Now we need to recognize to which root our expression corresponds. Is there a better way for this?
            roots = p_red.roots(ring=QQbar,multiplicities=False)
            for root in roots:
                root_approx = x.parent(root)
                LLL_basis = [x]
                for i in range(p_red.degree()):
                    LLL_basis.append(root_approx**i)
                LLL_res = lindep(LLL_basis,check_if_result_is_invalid=True)
                if LLL_res != None: #This seems to be the correct generator
                    K = NumberField(p_red,"v",embedding=root)
                    return K
            raise ArithmeticError("We should not get here!")

def lindep(L, check_if_result_is_invalid=True):
    """
    Given a list L, use PARI to return a list of factors F_i s.t. sum_i F_i*L_i = 0.
    If check_if_result_is_invalid==True, we perform several checks to test if the result is invalid.
    If solution is invalid, return None.
    """
    F = pari(L).lindep().list()
    res = [f.sage() for f in F]
    if check_if_result_is_invalid == True or len(res) == 0: #len(res) == 0 happens when len(L) == 2 because 1 and a complex number are independent over QQ...
        #We add additional constants to our list to check if LLL detects that these are not part of the solution
        parent = L[-1].parent()
        L_copy_results = []
        for c in [parent(pi),parent(pi,pi)]:
            L_copy = copy(L)
            L_copy.append(c)
            F = pari(L_copy).lindep().list()
            L_copy_res = [f.sage() for f in F]
            if L_copy_res[-1] != 0 or L_copy_res[:len(res)] != res: #Found an invalid example
                return None
            if len(L_copy_results) != 0 and L_copy_results[-1] != L_copy_res: #Found an invalid example
                return None
            L_copy_results.append(L_copy_res)
        return L_copy_results[0][:len(L)]
    return res

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

def get_u(floating_expression_linear_in_u, v, principal_cusp_width, extension_field_degree):
    """
    Determine u which is a factor that should cancel out common factors in the coefficients.
    This function returns u_interior_Kv which is an expression in Kv such that u = (u_interior_Kv)**(1/principal_cusp_width)
    and the embedding of u into CC.
    If principal_cusp_width == 1 we choose u to be in QQ (and hence independently of Kv).
    """
    x = floating_expression_linear_in_u**principal_cusp_width
    tmp = get_u_factor(x,v,principal_cusp_width,extension_field_degree)
    if tmp == None:
        return None
    c, recognized_expression = tmp
    if principal_cusp_width == 1: #In this case don't factor out "c"
        largest_denominator = max([Q.denominator() for Q in list(recognized_expression.polynomial())])
        u_interior_Kv = (QQ(1)/largest_denominator)*v**0
        CC = floating_expression_linear_in_u.parent()
        u_embedding = CC(u_interior_Kv)
    else:
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
