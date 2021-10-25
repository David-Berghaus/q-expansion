from sage.rings.complex_arb import ComplexBallField
from sage.rings.complex_mpc import MPComplexField
from sage.rings.real_mpfr import RealField
from sage.rings.qqbar import QQbar
# from sage.rings.complex_mpfr import ComplexField #This only works in Sage/9.3
from sage.interfaces.gp import gp

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

cpdef to_QQbar(x, max_extension_field_degree, check_if_poly_has_max_degree=False):
    """
    Tries to express x as an element of QQbar (this expression might not be correct).
    If check_if_poly_has_max_degree == True then also return if the polynomial is of maximal order which
    makes it likely that the expression is not a true solution.
    """
    numberfield = x.algebraic_dependency(max_extension_field_degree)
    roots = numberfield.roots(ring=QQbar,multiplicities=False)

    #Now we need to recognize to which root our expression corresponds. Is there a better way for this?
    diffs = [(root-x).abs() for root in roots]
    root_index = diffs.index(min(diffs))

    if check_if_poly_has_max_degree == False:
        return roots[root_index]
    else:
        return roots[root_index], numberfield.degree()==max_extension_field_degree

cpdef gp_lindep(x, int correct_digits, int max_order): #Not sure if this function will be useful
    bit_prec = digits_to_bits(correct_digits)
    CF = MPComplexField(bit_prec)
    x_cf = CF(x.real(),x.imag())
    x_str = x_cf.str()
    gp("x = " + x_str)
    gp_command = "lindep([1"
    for i in range(1,max_order+1):
        gp_command += ",x^" + str(i)
    gp_command += "])"
    lindep_res = gp(gp_command).sage()
    print("Warning, we have not set default(realprecision) correctly!")
    return lindep_res
