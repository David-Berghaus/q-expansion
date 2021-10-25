from sage.rings.qqbar import QQbar
from sage.modular.cusps import Cusp
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.power_series_poly import PowerSeries_poly
from sage.modular.modform.j_invariant import j_invariant_qexp
from sage.rings.complex_field import ComplexField
from sage.rings.laurent_series_ring import LaurentSeriesRing

from psage.modform.maass.automorphic_forms import AutomorphicFormSpace

from belyi.newton_genus_zero import run_newton
from point_matching.point_matching_arb_wrap import get_coefficients_haupt_ir_arb_wrap, digits_to_bits

def get_zero_multiplicity(zero, zeros):
    """
    Get multiplity of zero in a list of zeros that is of the form [(zero,multiplicity)].
    """
    for (zero_entry, multiplicity) in zeros:
        if zero == zero_entry:
            return multiplicity
    return 0

def get_klein_invariant_derivative(trunc_order):
    """
    Returns 1/(2*pi*i) * d/dtau j(tau).
    """
    j = j_invariant_qexp(trunc_order)
    q = j.parent().gen()
    j_prime = -1/q #Derivative of q^-1
    for n in range(1,j.degree()+1): #The derivative of the q^0 term is zero
        j_prime += n*j[n]*q**n
    return j_prime

def get_j_Gamma(S,digit_prec,trunc_order=None):
    """
    Get q-expansion of j_Gamma (we currently use Hejhal's method for this).
    We are currently only considering the cusp at infinity.
    """
    G = S.group()
    if G.cusp_width(Cusp(1,0)) != 1:
        raise NotImplementedError("We have not considered this case yet!")
    c, M = get_coefficients_haupt_ir_arb_wrap(S,digit_prec,only_principal_expansion=True,return_M=True)
    bit_prec = digits_to_bits(digit_prec)
    c = c._get_mcbd(bit_prec)
    CC = ComplexField(bit_prec)
    q = LaurentSeriesRing(CC,"q").gen()
    if trunc_order == None:
        trunc_order = M
    elif trunc_order > M:
        raise ArithmeticError("")
    
    j_Gamma = 1/q
    j_Gamma = j_Gamma.O(trunc_order)
    for j in range(1,trunc_order): #j_Gamma is normalized to have a zero constant term.
        j_Gamma += CC(c[j-1,0])*q**j
    
    return j_Gamma

def get_cuspform_q_expansion_from_hauptmodul(S, digit_prec):
    weight_half = S.weight()//2
    G = S.group()
    j_Gamma = get_j_Gamma(AutomorphicFormSpace(G,0),digit_prec)
    j_prime = get_klein_invariant_derivative(j_Gamma.degree())
    B = BelyiMap(AutomorphicFormSpace(G,0))

    R = B.get_cuspform_rational_functions(S.weight())[0]
    R_q_expansion = R.subs(x=j_Gamma)

    res = j_prime**weight_half*R_q_expansion
    return res

class BelyiMap():
    """
    Class for working with Belyi maps that are expressed in terms of Factored_Polynomials with algebraic coefficients.
    """
    def __init__(self, S, starting_digit_prec=42, target_digit_prec=50000, max_extension_field_degree=None):
        """
        Compute Belyi map from scratch. This is done by first using Hejhal's method to get approximate values at the elliptic points and cusps.
        Afterwards, Newton iterations are used to refine the solution. In the last step, the result is verified.
        """
        G = S.group()
        if G.genus() != 0:
            raise ArithmeticError("This function only works for genus zero subgroups!")
        if S.weight() != 0:
            raise ArithmeticError("This function only works for weight zero!")
        if max_extension_field_degree == None:
            max_extension_field_degree = G.index() #For the groups that we are considering the conjugacy class size is <= the index
        
        (p3, p2, pc) = run_newton(S, starting_digit_prec, target_digit_prec, stop_when_coeffs_are_recognized=True)
        self.G = G
        self.p3, self.p2, self.pc = p3, p2, pc
        self.p3_constructed, self.p2_constructed, self.pc_constructed = p3.construct(), p2.construct(), pc.construct()

        #If these 4 steps succeed, then the result is verified
        self._ell_2_point_evaluations = self._get_elliptic_two_point_evaluations()
        self._ell_3_point_evaluations = self._get_elliptic_three_point_evaluations()
        self._cusp_evaluations = self._get_cusp_evaluations()
        self.verify_polynomial_equation()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return self.p3.__str__() + " / " + self.pc.__str__()
    
    def _get_elliptic_two_point_evaluations(self):
        """
        Returns algebraic evaluations of hauptmodul at elliptic points of order two.
        """
        ell_2_point_evaluations = []
        for (p,multiplicity) in self.p2.factors:
            roots = p.roots(ring=QQbar)
            for (root, order) in roots:
                if order != 1:
                    raise ArithmeticError("Something went wrong, we need distinct roots here!")
                else:
                    ell_2_point_evaluations.append([root,multiplicity])
        return ell_2_point_evaluations

    def _get_elliptic_three_point_evaluations(self):
        """
        Returns algebraic evaluations of hauptmodul at elliptic points of order three.
        """
        ell_3_point_evaluations = []
        for (p,multiplicity) in self.p3.factors:
            roots = p.roots(ring=QQbar)
            for (root, order) in roots:
                if order != 1:
                    raise ArithmeticError("Something went wrong, we need distinct roots here!")
                else:
                    ell_3_point_evaluations.append([root,multiplicity])
        return ell_3_point_evaluations
    
    def _get_cusp_evaluations(self):
        """
        Returns algebraic evaluations of hauptmodul at cusps.
        The cusp at infinity gets treated separately.
        """
        cusp_evaluations = []
        for (p,multiplicity) in self.pc.factors:
            roots = p.roots(ring=QQbar)
            for (root, order) in roots:
                if order != 1:
                    raise ArithmeticError("Something went wrong, we need distinct roots here!")
                else:
                    cusp_evaluations.append([root,multiplicity])
        return cusp_evaluations

    def verify_polynomial_equation(self):
        p3_constructed, p2_constructed, pc_constructed = self.p3_constructed, self.p2_constructed, self.pc_constructed
        if p3_constructed-p2_constructed-1728*pc_constructed != 0: #Verify result
            raise ArithmeticError("Verification of polynomial_equation failed!")
        if pc_constructed.degree()+self.G.cusp_width(Cusp(1,0)) != self.G.index():
            raise ArithmeticError("Wrong behavior at infinity!")
    
    def get_rational_function(self):
        """
        Returns the quotient p3_constructed/pc_constructed.
        """
        return self.p3_constructed/self.pc_constructed
    
    def get_rational_function_derivative(self):
        """
        Returns the derivative of the quotient p3_constructed/pc_constructed.
        """
        R = self.get_rational_function()
        return R.derivative()
    
    # def _get_B_elliptic(self, weight):
        # """
        # Returns product (x-elliptic_evaluation)^beta_e for all elliptic points.
        # beta_e = weight/2 * (elliptic_order-1)
        # """
        # weight_half = weight//2
        # B_elliptic = 1
        # x = self.pc_constructed.parent().gen()
        # ell_2_point_evaluations, ell_3_point_evaluations = self._ell_2_point_evaluations, self._ell_3_point_evaluations
        
    #     for (ell_2_point_evaluation,multiplicity) in ell_2_point_evaluations:
    #         beta_e = weight_half*(2-1)
    #         power = multiplicity*beta_e
    #         B_elliptic *= (x-ell_2_point_evaluation)**power
    #     for (ell_3_point_evaluation,multiplicity) in ell_3_point_evaluations:
    #         beta_e = weight_half*(3-1)
    #         power = multiplicity*beta_e
    #         B_elliptic *= (x-ell_3_point_evaluation)**power
        
    #     return B_elliptic
    
    # def _get_B_cusp(self, weight):
    #     """
    #     Returns product (x-cusp_evaluation)^alpha_c for all cusps not equal infinity.
    #     alpha_c = weight/2 * ord(j_Gamma_prime(c))
    #     """
    #     weight_half = weight//2
    #     R_prime = self.get_rational_function_derivative()
    #     B_cusp = R_prime.denominator()**weight_half
    #     return B_cusp

    #     # #Note that j_Gamma_prime = j_prime/R_prime so we need to get the valuations of R_prime at the cusps
    #     # R_prime = self.get_rational_function_derivative()
    #     # #Somehow we have to call roots instead of zeros here...
    #     # num_zeros, den_zeros = R_prime.numerator().roots(ring=QQbar), R_prime.denominator().roots(ring=QQbar)

    #     # cusp_evaluations = self._cusp_evaluations
    #     # j_Gamma_prime_valuations = []
    #     # for (cusp_evaluation,cusp_evaluation_multiplicity) in cusp_evaluations:
    #     #     num_zero_order, den_zero_order = 0, 0
    #     #     for (num_zero,multiplicity) in num_zeros:
    #     #         if num_zero == cusp_evaluation:
    #     #             num_zero_order = multiplicity
    #     #     for (den_zero,multiplicity) in den_zeros:
    #     #         if den_zero == cusp_evaluation:
    #     #             den_zero_order = multiplicity
    #     #     cusp_evaluation_order = -1-(num_zero_order-den_zero_order) #The -1 comes from the pole of j_prime
    #     #     if cusp_evaluation_order < 0:
    #     #         raise ArithmeticError("This should not happen...")
    #     #     else:
    #     #         j_Gamma_prime_valuations.append([cusp_evaluation,cusp_evaluation_order])
        
    #     # weight_half = weight//2
    #     # B_cusp = 1
    #     # x = self.pc_constructed.parent().gen()
    #     # for (j_Gamma_prime_evaluation,order) in j_Gamma_prime_valuations:
    #     #     alpha_c = weight_half*order
    #     #     power = alpha_c
    #     #     B_cusp *= (x-j_Gamma_prime_evaluation)**power

    #     # return B_cusp
    
    def _get_B_elliptic(self, weight):
        """
        Returns product (x-elliptic_evaluation)^beta_e for all elliptic points.
        beta_e = weight/2 * (elliptic_order-1)
        """
        weight_half = weight//2
        B_elliptic = 1
        x = self.p2_constructed.parent().gen()
        ell_2_point_evaluations, ell_3_point_evaluations = self._ell_2_point_evaluations, self._ell_3_point_evaluations

        R_prime = self.get_rational_function_derivative()
        #Somehow we have to call roots instead of zeros here...
        numerator_zeros = R_prime.numerator().roots(ring=QQbar)

        for (ell_2_point_evaluation,multiplicity) in ell_2_point_evaluations:
            beta_e = weight_half*((2-1) - get_zero_multiplicity(ell_2_point_evaluation,numerator_zeros))
            if beta_e < 0:
                raise ArithmeticError("This should not happen...")
            B_elliptic *= (x-ell_2_point_evaluation)**beta_e
        for (ell_3_point_evaluation,multiplicity) in ell_3_point_evaluations:
            beta_e = weight_half*((3-1) - get_zero_multiplicity(ell_3_point_evaluation,numerator_zeros))
            if beta_e < 0:
                raise ArithmeticError("This should not happen...")
            B_elliptic *= (x-ell_3_point_evaluation)**beta_e
        
        return B_elliptic

    def _get_B_cusp(self, weight):
        """
        Returns product (x-cusp_evaluation)^alpha_c for all cusps not equal infinity.
        alpha_c = weight/2 * ord(j_Gamma_prime(c))
        """
        weight_half = weight//2
        B_cusp = 1
        x = self.pc_constructed.parent().gen()
        cusp_evaluations = self._cusp_evaluations

        R_prime = self.get_rational_function_derivative()
        #Somehow we have to call roots instead of zeros here...
        denominator_zeros = R_prime.denominator().roots(ring=QQbar)

        for (cusp_evaluation,multiplicity) in cusp_evaluations:
            alpha_c = weight_half*(get_zero_multiplicity(cusp_evaluation,denominator_zeros) - 1)
            if alpha_c < 0:
                raise ArithmeticError("This should not happen...")
            B_cusp *= (x-cusp_evaluation)**alpha_c
        
        return B_cusp
        
    def _get_B(self, weight):
        """
        Returns B_cusp*B_elliptic
        """
        B_cusp = self._get_B_cusp(weight)
        B_elliptic = self._get_B_elliptic(weight)
        B = B_cusp*B_elliptic
        return B

    def _get_P_list_cuspform(self, weight):
        """
        Return list of polynomials of the form prod (x-cusp_evaluation) that allows to prescribe zeros of given order at the cusps.
        We assume that P vanishes to degree 1 at all cusps not equal infinity.
        For the cusp at infinity we assume an order of vanishing in reduced-row echelon form.
        """
        weight_half = weight//2
        x = self.pc_constructed.parent().gen()
        cusp_evaluations = self._cusp_evaluations
        cuspform_dim = self.G.dimension_cusp_forms(weight)
        if cuspform_dim == 0:
            raise ArithmeticError("The dimension for cuspforms is zero for this space!")
        B = self._get_B(weight)
        P_list = [1 for _ in range(cuspform_dim)]

        for n in range(1,cuspform_dim+1):
            P_degree = -weight_half+B.degree()-n
            for (cusp_evaluation,_) in cusp_evaluations:
                P_list[n-1] *= (x-cusp_evaluation)
            #Now we need to get the correct order of vanishing at infinity
            power = P_degree-len(cusp_evaluations)
            if power < 0:
                raise ArithmeticError("This should not happen...")
            P_list[n-1] *= x**power
        
        return P_list
    
    def get_cuspform_rational_functions(self, weight):
        """
        Return rational functions R, such that (j')^(weight/2)*R gives a basis of cuspforms.
        """
        weight_half = weight//2
        P_list = self._get_P_list_cuspform(weight)
        R_prime_pow_weight_half = self.get_rational_function_derivative()**weight_half
        B = self._get_B(weight)
        cuspform_rational_functions = []

        tmp = 1/(R_prime_pow_weight_half*B)
        for P in P_list:
            R = tmp*P
            cuspform_rational_functions.append(R)

        return cuspform_rational_functions

    def get_x_expansion(self, trunc_order):
        """
        Get expansion of rational function in terms of x (the hauptmodul).
        """
        power_series_ring = PowerSeriesRing(QQbar, "x")
        O = PowerSeries_poly(power_series_ring,prec=trunc_order) #This corresponds to O(x**trunc_order) in sage syntax
        p3_series = self.p3_constructed+O
        pc_series = self.pc_constructed+O
        R_series = p3_series/pc_series
        return R_series

    def get_q_expansion(self, trunc_order):
        raise NotImplementedError("This functionality has not been implemented yet!")
    
    def get_approx_q_expansion(self, trunc_order, digit_prec):
        raise NotImplementedError("This functionality has not been implemented yet!")
