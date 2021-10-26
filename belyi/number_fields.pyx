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

from point_matching.point_matching_arb_wrap import digits_to_bits

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

# cpdef to_QQbar(x, max_extension_field_degree, check_if_poly_has_max_degree=False):
#     """
#     Tries to express x as an element of QQbar (this expression might not be correct).
#     If check_if_poly_has_max_degree == True then also return if the polynomial is of maximal order which
#     makes it likely that the expression is not a true solution.
#     """
#     numberfield = x.algebraic_dependency(max_extension_field_degree)
#     roots = numberfield.roots(ring=QQbar,multiplicities=False)

#     #Now we need to recognize to which root our expression corresponds. Is there a better way for this?
#     diffs = [(root-x).abs() for root in roots]
#     root_index = diffs.index(min(diffs))

#     if check_if_poly_has_max_degree == False:
#         return roots[root_index]
#     else:
#         return roots[root_index], numberfield.degree()==max_extension_field_degree

cpdef to_QQbar(x, gen, extension_field_degree, check_if_result_is_invalid=True):
    """
    Try to express x as an algebraic number in terms of gen.

    If check_if_result_is_invalid == True then also place pi into LLL-basis (therefore, if the coefficient of pi is non-zero
    then the detected solution is not valid) and if result is invalid, return False.
    """
    gen_approx = x.parent(gen)
    LLL_basis = [x]
    for i in range(extension_field_degree):
        LLL_basis.append(gen_approx**i)
    if check_if_result_is_invalid == True:
        LLL_basis.append(x.parent(pi))
    LLL_res = lindep(LLL_basis)
    if check_if_result_is_invalid == True and LLL_res[-1] != 0: #Result cannot be correct because pi is not part of true basis
        return False
    
    x_alg = 0
    for i in range(extension_field_degree):
        x_alg += LLL_res[i+1]*gen**i
    x_alg /= -LLL_res[0]
    return x_alg   

def get_numberfield_and_gen(x, max_extension_field_degree, check_if_result_is_invalid=True, reduce_numberfield=True):
    """
    Try to express x in a numberfield over QQbar with specified max_extension_field_degree.
    This function returns a polynomial over which the numberfield is defined as well as the generator (expressed as an algebraic number).

    If check_if_result_is_invalid == True then also place pi into LLL-basis (therefore, if the coefficient of pi is non-zero
    then the detected solution is not valid) and if result is invalid, return False.

    If reduce_numberfield == True then try to reduce the numberfield and return the generator of this numberfield.
    It is advised to only use this option if one is certain that the result is valid, because this function might take very long for
    spurious numberfields.
    """
    LLL_basis = [x**i for i in range(max_extension_field_degree+1)]
    if check_if_result_is_invalid == True:
        LLL_basis.append(x.parent(pi))
    LLL_res = lindep(LLL_basis)
    if check_if_result_is_invalid == True and LLL_res[-1] != 0: #Result cannot be correct because pi is not part of true basis
        return False
    
    P = PolynomialRing(ZZ,"x")
    numberfield = P(LLL_res)
    if reduce_numberfield == False:
        #Now we need to recognize to which root our expression corresponds. Is there a better way for this?
        roots = numberfield.roots(ring=QQbar,multiplicities=False)
        diffs = [(root-x).abs() for root in roots]
        root_index = diffs.index(min(diffs))
        gen = roots[root_index]
        return numberfield, gen
    else:
        numberfield_red = P(pari.polredabs(numberfield))
        if numberfield_red.degree() == 1:
            return numberfield_red, 1
        else:
            #Now we need to recognize to which root our expression corresponds. Is there a better way for this?
            roots = numberfield_red.roots(ring=QQbar,multiplicities=False)
            for root in roots:
                root_approx = x.parent(root)
                LLL_basis = [x]
                for i in range(numberfield_red.degree()):
                    LLL_basis.append(root_approx**i)
                LLL_basis.append(x.parent(pi))
                LLL_res = lindep(LLL_basis)
                if LLL_res[-1] == 0: #This seems to be the correct generator
                    return numberfield_red, root
            raise ArithmeticError("We should not get here!")

def lindep(L):
    """
    Given a list L, use PARI to return a list of factors F_i s.t. sum_i F_i*L_i = 0.
    """
    F = pari(L).lindep().list()
    res = [f.sage() for f in F]
    return res
