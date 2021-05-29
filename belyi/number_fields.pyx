from sage.rings.complex_arb import ComplexBallField
from sage.rings.complex_mpc import MPComplexField
from sage.rings.real_mpfr import RealField
# from sage.rings.complex_mpfr import ComplexField #This only works in Sage/9.3
from sage.interfaces.gp import gp

from point_matching.point_matching_arb_wrap import digits_to_bits

cpdef get_decimal_digit_prec(epsilon):
    """
    Given a real number epsilon, return the decimal digit precision.
    """
    RF = RealField(53)
    eps_approx = RF(epsilon) #There is no reason to perform this computation at full precision
    return -int(epsilon.log10())

cpdef algebraic_dependency(x, int correct_digits, int max_order):
    bit_prec = digits_to_bits(correct_digits)
    CF = MPComplexField(bit_prec)
    x_cf = CF(x.real(),x.imag())
    return x_cf.algebraic_dependency(max_order)

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
    return lindep_res
