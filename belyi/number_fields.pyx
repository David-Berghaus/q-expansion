from copy import copy

from sage.rings.complex_arb import ComplexBallField
from sage.rings.complex_mpc import MPComplexField
from sage.rings.real_mpfr import RealField
from sage.rings.qqbar import QQbar
# from sage.rings.complex_mpfr import ComplexField #This only works in Sage/9.3
from sage.interfaces.gp import gp
from sage.interfaces.gp import pari
from sage.symbolic.constants import pi
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.number_field.number_field import NumberField

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
