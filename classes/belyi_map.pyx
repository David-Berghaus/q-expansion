from math import ceil

from sage.libs.arb.acb_poly cimport *
from sage.rings.qqbar import QQbar
from sage.modular.cusps import Cusp
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.power_series_poly import PowerSeries_poly
from sage.modular.modform.j_invariant import j_invariant_qexp
from sage.rings.complex_field import ComplexField
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.real_mpfr import RealField
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.polynomial.polynomial_complex_arb import Polynomial_complex_arb
from sage.rings.complex_arb import ComplexBallField

from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from psage.modform.arithgroup.mysubgroup import MySubgroup

from belyi.newton_genus_zero import run_newton
from point_matching.point_matching_arb_wrap import get_coefficients_haupt_ir_arb_wrap, digits_to_bits
from classes.fourier_expansion import FourierExpansion, to_reduced_row_echelon_form

def get_j_Gamma_hejhal(S,digit_prec,trunc_order=None):
    """
    Get q-expansion of j_Gamma (we currently use Hejhal's method for this).
    We are currently only considering the cusp at infinity.
    """
    G = S.group()
    print("Note that we absorb the cusp-width into q!")
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

def get_n_th_root_of_1_over_j(trunc_order,n): #We work over QQ because it is usually not slower than approx and one does not need to worry about conditioning
    """
    Returns (1/j)^(1/n) up to trunc_order terms.
    The result is a power series in q_n where q_n = exp(2*pi*I/n).
    """
    N = int(ceil(float(trunc_order)/n)) #q^N = q_n^trunc_order
    j = j_invariant_qexp(N)
    j_inv = j.inverse()
    var_name = "q_" + str(n)
    L = PowerSeriesRing(j[0].parent(),var_name)
    q_n = L.gen()
    #.subs() seems quite slow, maybe try to write faster cython code
    tmp = L(j_inv.subs(q=q_n**n)) #Because we currently cannot work with Puiseux series in Sage.
    res = tmp.nth_root(n)
    return res

def get_B_factored_degree(B_factored):
    """
    Given B_factored, as returned by _get_B_factored(), return the degree in x (i.e. the sum of all orders).
    """
    degree = 0
    for (_, order) in B_factored:
        degree += order
    return degree

def identify_cusp_from_pc_root(pc_root, cusp_rep_values):
    """
    Given an algebraic expression of a root of pc, try to recognize the cusp corresponding to this root by 
    comparing it to the numerical output that we got from Hejhal's method.
    """
    diffs = []
    for (_,cusp_rep_value) in cusp_rep_values:
        diffs.append((pc_root-cusp_rep_value).abs())
    cusp,_ = cusp_rep_values[diffs.index(min(diffs))]
    return cusp

class BelyiMap():
    """
    Class for working with Belyi maps that are expressed in terms of Factored_Polynomials with algebraic coefficients.
    """
    def __init__(self, G, starting_digit_prec=42, target_digit_prec=50000, max_extension_field_degree=None):
        """
        Compute Belyi map from scratch. This is done by first using Hejhal's method to get approximate values at the elliptic points and cusps.
        Afterwards, Newton iterations are used to refine the solution. In the last step, the result is verified.
        """
        if G.genus() != 0:
            raise ArithmeticError("This function only works for genus zero subgroups!")
        if max_extension_field_degree == None:
            print("We are using the index as max_extension_field_degree which is not correct in general!")
            max_extension_field_degree = G.index()
        
        G = MySubgroup(G)
        S = AutomorphicFormSpace(G,0)
        (p3, p2, pc), cusp_rep_values = run_newton(S,starting_digit_prec,target_digit_prec,stop_when_coeffs_are_recognized=True,return_cusp_rep_values=True,max_extension_field_degree=max_extension_field_degree)
        self.G = G
        self.p3, self.p2, self.pc = p3, p2, pc
        self.p3_constructed, self.p2_constructed, self.pc_constructed = p3.construct(), p2.construct(), pc.construct()
        self.princial_cusp_width = G.cusp_width(Cusp(1,0))

        self._e2_valuations = self._get_e2_fixed_point_valuations()
        self._e3_valuations = self._get_e3_fixed_point_valuations()
        self._cusp_valuations = self._get_cusp_valuations(cusp_rep_values)
        self.verify_polynomial_equation()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return self.p3.__str__() + " / " + self.pc.__str__()
    
    def _get_e2_fixed_point_valuations(self):
        """
        Returns algebraic valuation of hauptmodul at elliptic fixed points of order two.
        """
        for (p,multiplicity) in self.p2.factors:
            if multiplicity == 1: #We are only interested in the fixed points
                e2_valuations = []
                roots = p.roots(ring=QQbar)
                for (root, order) in roots:
                    if order != 1:
                        raise ArithmeticError("Something went wrong, we need distinct roots here!")
                    e2_valuations.append(root)
                return e2_valuations
        return []

    def _get_e3_fixed_point_valuations(self):
        """
        Returns algebraic valuation of hauptmodul at elliptic fixed points of order three.
        """
        for (p,multiplicity) in self.p3.factors:
            if multiplicity == 1: #We are only interested in the fixed points
                e3_valuations = []
                roots = p.roots(ring=QQbar)
                for (root, order) in roots:
                    if order != 1:
                        raise ArithmeticError("Something went wrong, we need distinct roots here!")
                    e3_valuations.append(root)
                return e3_valuations
        return []
    
    def _get_cusp_valuations(self, cusp_rep_values):
        """
        Returns the cusp representatives and the algebraic valuations of hauptmodul attached to these.
        The cusp at infinity gets treated separately.
        The parameter "cusp_rep_values" contains contains floating point approximations for each cusp which we try to use
        to attach each cusp to a root of pc.
        """
        G = self.G
        cusp_valuations = dict()
        for (p,multiplicity) in self.pc.factors:
            roots = p.roots(ring=QQbar)
            for (root, order) in roots:
                if order != 1:
                    raise ArithmeticError("Something went wrong, we need distinct roots here!")
                cusp = identify_cusp_from_pc_root(root,cusp_rep_values)
                if G.cusp_width(cusp) != multiplicity:
                    raise ArithmeticError("This should not happen!")
                cusp_valuations[cusp] = root

        if len(cusp_valuations) != G.ncusps()-1: #Some cusps have multiple values attached to them which is invalid
            raise ArithmeticError("This should not happen!")

        return cusp_valuations

    def verify_polynomial_equation(self):
        p3_constructed, p2_constructed, pc_constructed = self.p3_constructed, self.p2_constructed, self.pc_constructed
        if p3_constructed-p2_constructed-1728*pc_constructed != 0: #Verify result
            raise ArithmeticError("Verification of polynomial_equation failed!")
        if pc_constructed.degree()+self.princial_cusp_width != self.G.index():
            raise ArithmeticError("Wrong behavior at infinity!")
    
    def _get_B_elliptic_factored(self, weight):
        """
        Returns factors (x-elliptic_evaluation)^beta_e for all elliptic fixed points.
        beta_e is chosen to cancel the zeros of j_G'^weight_half.
        """
        weight_half = weight//2
        x = self.p2_constructed.parent().gen()
        B_elliptic_factors = []
        e2_valuations, e3_valuations = self._e2_valuations, self._e3_valuations
        for e2_valuation in e2_valuations:
            beta_e = weight_half*(2-1)//2 #We need to divide because (x-j_G(e_2)) has a double zero
            B_elliptic_factors.append([x-e2_valuation,beta_e])
        for e3_valuation in e3_valuations:
            beta_e = weight_half*(3-1)//3 #We need to divide because (x-j_G(e_3)) has a triple zero
            B_elliptic_factors.append([x-e3_valuation,beta_e])
        return B_elliptic_factors

    def _get_B_cusp_factored(self, weight):
        """
        Returns factors (x-cusp_evaluation)^alpha_c for all cusps not equal infinity.
        alpha_c is chosen to cancel the zeros of j_G'^weight_half.
        """
        weight_half = weight//2
        B_cusp_factors = []
        x = self.pc_constructed.parent().gen()
        cusp_valuations = self._cusp_valuations
        for (_,cusp_valuation) in cusp_valuations.items():
            alpha_c = weight_half
            B_cusp_factors.append([x-cusp_valuation,alpha_c])
        return B_cusp_factors
        
    def _get_B_factored(self, weight):
        """
        Returns B_cusp*B_elliptic in factored form.
        """
        B_cusp_factors = self._get_B_cusp_factored(weight)
        B_elliptic_factors = self._get_B_elliptic_factored(weight)
        B = B_cusp_factors+B_elliptic_factors #This is the syntax for merging lists in python
        return B

    def _get_p_list_cuspform(self, weight, B_factored):
        """
        This function returns a list of polynomials p such that for each p in list, (j_G'^weight_half)*p/B corresponds to a cuspform.
        p is written in factors of (x-cusp_evaluation) which prescribe zeros of given order at the cusps.
        We assume that p vanishes to degree 1 at all cusps not equal infinity.
        For the cusp at infinity we assume orders of vanishing up to the dimension of the space.
        """
        weight_half = weight//2
        x = self.pc_constructed.parent().gen()
        cusp_valuations = self._cusp_valuations
        cuspform_dim = self.G.dimension_cusp_forms(weight)
        if cuspform_dim == 0:
            raise ArithmeticError("The dimension of cuspforms is zero for this space!")
        B_factored_degree = get_B_factored_degree(B_factored)
        p_list = [[] for _ in range(cuspform_dim)]

        for n in range(1,cuspform_dim+1):
            for (_,cusp_valuation) in cusp_valuations.items():
                p_list[n-1].append([x-cusp_valuation,1]) #Prescribe zeros of order one at all cusps != infty
            #Now we need to get the correct order of vanishing at infinity
            p_degree = -weight_half+B_factored_degree-n
            power = p_degree-len(cusp_valuations) #Exponent of x s.t. cuspform vanishes to degree n at infty.
            if power < 0:
                raise ArithmeticError("This should not happen...")
            p_list[n-1].append([x,power])
        
        return p_list
    
    def _get_p_list_modform(self, weight, B_factored):
        """
        This function returns a list of polynomials p such that for each p in list, (j_G'^weight_half)*p/B corresponds to a holomorphic modform.
        p is written in factors of (x-cusp_evaluation) which prescribe zeros of given order at the cusps.
        We assume that p is constant and non-zero at all cusps not equal infinity.
        For the cusp at infinity we assume orders of vanishing up to the dimension of the space.
        """
        weight_half = weight//2
        x = self.pc_constructed.parent().gen()
        modform_dim = self.G.dimension_modular_forms(weight)
        if modform_dim == 0:
            raise ArithmeticError("The dimension of modforms is zero for this space!")
        B_factored_degree = get_B_factored_degree(B_factored)
        p_list = [[] for _ in range(modform_dim)]

        for n in range(1,modform_dim+1):
            #Now we need to get the correct order of vanishing at infinity
            p_degree = -weight_half+B_factored_degree-n+1
            power = p_degree #Exponent of x s.t. modform vanishes to degree n-1 at infty.
            if power < 0:
                raise ArithmeticError("This should not happen...")
            p_list[n-1].append([x,power])
        
        return p_list
    
    def get_cuspforms(self, weight, trunc_order, digit_prec=None):
        """
        Use Hauptmodul to return a basis of cuspforms with specified weight in reduced row-echelon form.
        If "digit_prec" is given, use approximate arithmetic.
        """
        cusp = Cusp(1,0) #Add functionality for other cusps later
        if digit_prec == None:
            j_G = self.get_hauptmodul_q_expansion(cusp,trunc_order) #We could precompute this
        else:
            j_G = self.get_hauptmodul_q_expansion_approx(cusp,trunc_order,digit_prec) #We could precompute this
        B_factored = self._get_B_factored(weight)
        F = self._get_regularized_modular_form_q_expansion(weight,j_G,B_factored) #We could re-use this for modforms of the same weight
        p_list = self._get_p_list_cuspform(weight,B_factored)
        ring = j_G[1].parent()
        M = MatrixSpace(ring,len(p_list),trunc_order)
        A = []

        for p in p_list:
            cuspform = F
            for (factor, order) in p:
                x = factor.parent().gen()
                cuspform *= factor.subs(x=j_G)**order
            A.append([cuspform[i] for i in range(trunc_order)]) #Generate matrix of coefficients
        
        A = M(A).echelon_form()
        P = PowerSeriesRing(ring,j_G.variable())
        cuspforms = []
        for i in range(len(p_list)):
            cuspform = P(A[i,:].list())
            cuspform = cuspform.O(trunc_order)
            cuspforms.append(cuspform)

        return cuspforms
    
    def get_modforms(self, weight, trunc_order, digit_prec=None):
        """
        Use Hauptmodul to return a basis of modforms with specified weight in reduced row-echelon form.
        If "digit_prec" is given, use approximate arithmetic.
        """
        cusp = Cusp(1,0) #Add functionality for other cusps later
        if digit_prec == None:
            j_G = self.get_hauptmodul_q_expansion(cusp,trunc_order) #We could precompute this
        else:
            j_G = self.get_hauptmodul_q_expansion_approx(cusp,trunc_order,digit_prec) #We could precompute this
        B_factored = self._get_B_factored(weight)
        F = self._get_regularized_modular_form_q_expansion(weight,j_G,B_factored) #We could re-use this for cuspforms of the same weight
        p_list = self._get_p_list_modform(weight,B_factored)
        ring = j_G[1].parent()
        M = MatrixSpace(ring,len(p_list),trunc_order)
        A = []

        for p in p_list:
            modform = F
            for (factor, order) in p:
                x = factor.parent().gen()
                modform *= factor.subs(x=j_G)**order
            A.append([modform[i] for i in range(trunc_order)]) #Generate matrix of coefficients
        
        A = M(A).echelon_form()
        P = PowerSeriesRing(ring,j_G.variable())
        modforms = []
        for i in range(len(p_list)):
            modform = P(A[i,:].list())
            modform = modform.O(trunc_order)
            modforms.append(modform)

        return modforms

    def get_hauptmodul_q_expansion(self, trunc_order):
        """
        Return q-expansion of hauptmodul at all cusps using rigorous arithmetic.
        We return the result as an instance of FourierExpansion.
        """
        cusp_expansions = dict()
        for cusp in self.G.cusps():
            if cusp == Cusp(1,0):
                cusp_expansion = self._get_hauptmodul_q_expansion_infinity(trunc_order)
            else:
                cusp_expansion = self._get_hauptmodul_q_expansion_non_infinity(cusp,trunc_order)
            cusp_expansions[cusp] = cusp_expansion
        return FourierExpansion(self.G,0,cusp_expansions,"Hauptmodul")
    
    def get_hauptmodul_q_expansion_approx(self, trunc_order, digit_prec, try_to_overcome_ill_conditioning=True):
        """
        Return q-expansion of hauptmodul at all cusps using floating point arithmetic.
        We return the result as an instance of FourierExpansion.
        Note that not all coefficients need to be correct up to the specified precision!
        """
        cusp_expansions = dict()
        for cusp in self.G.cusps():
            if cusp == Cusp(1,0):
                cusp_expansion = self._get_hauptmodul_q_expansion_infinity_approx(trunc_order,digit_prec,try_to_overcome_ill_conditioning=try_to_overcome_ill_conditioning)
            else:
                cusp_expansion = self._get_hauptmodul_q_expansion_non_infinity_approx(cusp,trunc_order,digit_prec,try_to_overcome_ill_conditioning=try_to_overcome_ill_conditioning)
            cusp_expansions[cusp] = cusp_expansion
        return FourierExpansion(self.G,0,cusp_expansions,"Hauptmodul")

    def _get_hauptmodul_q_expansion_infinity(self, trunc_order):
        """
        Computes the q-expansion at the principal cusp which is normalized to j_Gamma = 1/q_N + 0 + c_1*q_N + ...
        """
        parent = self.p2_constructed[0].parent()
        L = LaurentSeriesRing(parent,"x")
        x = L.gen()
        princial_cusp_width = self.princial_cusp_width
        s = (L(self.pc_constructed).subs({x:1/x}).O(trunc_order)/L(self.p3_constructed).subs({x:1/x}).O(trunc_order)).power_series().nth_root(self.princial_cusp_width)
        r = s.reverse().inverse()
        n_sqrt_j_inverse = get_n_th_root_of_1_over_j(trunc_order,princial_cusp_width)
        j_G = r.subs(x=n_sqrt_j_inverse)
        return j_G

    def _get_hauptmodul_q_expansion_non_infinity(self, cusp, trunc_order):
        """
        Computes the q-expansion at a non-principal cusp which is normalized to j_Gamma = c_0 + c_1*q_N + ...
        """
        cusp_evaluation = self._cusp_valuations[cusp]
        cusp_width = self.G.cusp_width(cusp)
        parent = self.p2_constructed[0].parent()
        L = LaurentSeriesRing(parent,"x")
        x = L.gen()
        s = (L(self.pc_constructed.subs(x=x+cusp_evaluation)).O(trunc_order)/L(self.p3_constructed.subs(x=x+cusp_evaluation)).O(trunc_order)).power_series().nth_root(cusp_width)
        r = s.reverse()
        n_sqrt_j_inverse = get_n_th_root_of_1_over_j(trunc_order,cusp_width)
        j_G = cusp_evaluation + r.subs(x=n_sqrt_j_inverse)
        return j_G

    # def _get_hauptmodul_q_expansion_infinity_approx_sage(self, trunc_order, digit_prec):
    #     """
    #     Computes q-expansion of hauptmodul using Sage implementations.
    #     These are significantly slower than the ones provided by Arb and it is hence advised to use "get_hauptmodul_approx_q_expansion" instead.
    #     """
    #     bit_prec = digits_to_bits(digit_prec)
    #     CC = ComplexField(bit_prec)
    #     L = LaurentSeriesRing(CC,"x")
    #     x = L.gen()
    #     princial_cusp_width = self.princial_cusp_width
    #     s = (L(self.pc_constructed).subs({x:1/x}).O(trunc_order)/L(self.p3_constructed).subs({x:1/x}).O(trunc_order)).power_series().nth_root(self.princial_cusp_width)
    #     r = s.reverse().inverse()
    #     n_sqrt_j_inverse = get_n_th_root_of_1_over_j(trunc_order,princial_cusp_width)
    #     j_G = r.subs({x:n_sqrt_j_inverse})
    #     return j_G #Note that some of the higher coefficients might be wrong due to rounding errors

    def _get_r_for_laurent_expansion(self, trunc_order, bit_prec):
        """
        Get the reversed series of the reciprocal of the Belyi map expanded in 1/x.
        We need to compute this for "get_hauptmodul_q_expansion_infinity_approx".
        """
        CC = ComplexField(bit_prec)
        CBF = ComplexBallField(bit_prec)
        L = LaurentSeriesRing(CC,"x")
        x = L.gen()
        #We need to perform these computations with CC because somehow nth_root does not work with CBF,
        #although I would have expected these routines to be generic...
        #Maybe add these functions in the future so that we can output a result with rigorous error bounds (and potentially performance gains).
        s = (L(self.pc_constructed).subs({x:1/x}).O(trunc_order)/L(self.p3_constructed).subs({x:1/x}).O(trunc_order)).power_series().nth_root(self.princial_cusp_width)
        s_prec = s.prec() #Exponent of the O-term
        s_arb_reverted = s.polynomial().change_ring(CBF).revert_series(s_prec) #Perform the reversion in arb because it is expensive
        r = L(s_arb_reverted).O(s_prec).inverse()
        return r

    def _get_hauptmodul_q_expansion_infinity_approx(self, trunc_order, digit_prec, try_to_overcome_ill_conditioning=True):
        """
        Compute approximation of the coefficients of the q-expansion of the hauptmodul truncated to "trunc_order" and with "digit_prec"
        working precision. We make use of fast Arb implementations for performance.
        Because of the large coefficients involved, the arithmetic might become ill-conditioned. 
        If "try_to_overcome_ill_conditioning" == True, we try to detect these cases and increase the
        working precision if required (still, the higher coefficients will in general not have the full displayed precision).
        """
        CC_res = ComplexField(digits_to_bits(digit_prec))
        princial_cusp_width = self.princial_cusp_width
        n_sqrt_j_inverse = get_n_th_root_of_1_over_j(trunc_order,princial_cusp_width)
        if try_to_overcome_ill_conditioning == True:
            #Now guess the minimal precision required to get the correct order of magnitude of the last coefficient
            #We do this by constructing "r" to low precision to get the size of its largest exponent
            r_low_prec = self._get_r_for_laurent_expansion(trunc_order,64)
            required_prec = int(round(r_low_prec[r_low_prec.degree()].abs().log10()))
            working_prec = max(digit_prec,required_prec)
            if working_prec > digit_prec:
                print("Used higher digit precision during Hauptmodul q-expansion computation: ", working_prec)
        else:
            working_prec = digit_prec

        working_bit_prec = digits_to_bits(working_prec)
        CBF = ComplexBallField(working_bit_prec)
        r = self._get_r_for_laurent_expansion(trunc_order,working_bit_prec)

        n_sqrt_j_inverse_CBF = n_sqrt_j_inverse.polynomial().change_ring(CBF)
        r_pos_degree = r[0:r.degree()+1].power_series().polynomial().change_ring(CBF) #We treat the 1/x term later because arb only supports polynomials
        tmp = r_pos_degree.compose_trunc(n_sqrt_j_inverse_CBF,trunc_order) #Perform composition in arb because it is expensive
        P = PowerSeriesRing(CBF,n_sqrt_j_inverse_CBF.variable_name())
        one_over_x_term = P(n_sqrt_j_inverse_CBF).O(n_sqrt_j_inverse_CBF.degree()+1).inverse() #Treat 1/x term separately because it is not supported by arb
        j_G = one_over_x_term + P(tmp)

        return j_G.change_ring(CC_res)

    def _get_r_for_taylor_expansion(self, cusp, trunc_order, bit_prec):
        """
        Get the reversed series of the reciprocal of the Belyi map expanded in x.
        We need to compute this for "get_hauptmodul_q_expansion_non_infinity_approx".
        """
        #We need to perform these computations with CC because somehow nth_root does not work with CBF,
        #although I would have expected these routines to be generic...
        #Maybe add these functions in the future so that we can output a result with rigorous error bounds (and potentially performance gains).
        cusp_evaluation = self._cusp_valuations[cusp]
        cusp_width = self.G.cusp_width(cusp)
        CC = ComplexField(bit_prec)
        CBF = ComplexBallField(bit_prec)
        L = LaurentSeriesRing(CC,"x")
        x = L.gen()
        cusp_evaluation_CC = CC(cusp_evaluation)
        pc_shifted, p3_shifted = L(self.pc_constructed).subs({x:x+cusp_evaluation_CC}).O(trunc_order), L(self.p3_constructed).subs({x:x+cusp_evaluation_CC}).O(trunc_order)
        
        #It is very important that the leading order terms of pc_shifted are truely zero, otherwise we can get very large rounding errors
        #We therefore set the coefficients that are effectively zero to true zeros
        pc_shifted_coeffs = []
        for i in range(pc_shifted.degree()+1):
            if i < cusp_width:
                pc_shifted_coeffs.append(CC(0))
            else:
                pc_shifted_coeffs.append(pc_shifted[i])
        pc_shifted = L(pc_shifted_coeffs).O(trunc_order)

        s = (pc_shifted/p3_shifted).power_series().nth_root(cusp_width)
        s_prec = s.prec() #Exponent of the O-term
        s_arb_reverted = s.polynomial().change_ring(CBF).revert_series(s_prec) #Perform the reversion in arb because it is expensive
        s_arb_reverted_prec = s_arb_reverted.prec() #Exponent of the O-term
        r = L(s_arb_reverted).O(s_arb_reverted_prec)

        return r

    def _get_hauptmodul_q_expansion_non_infinity_approx(self, cusp, trunc_order, digit_prec, try_to_overcome_ill_conditioning=True):
        """
        Compute approximation of the coefficients of the q-expansion of the hauptmodul at a non-principal cusp truncated to "trunc_order" 
        and with "digit_prec" working precision. We make use of fast Arb implementations for performance.
        Because of the large coefficients involved, the arithmetic might become ill-conditioned. 
        If "try_to_overcome_ill_conditioning" == True, we try to detect these cases and increase the
        working precision if required (still, the higher coefficients will in general not have the full displayed precision).
        """
        CC_res = ComplexField(digits_to_bits(digit_prec))

        if try_to_overcome_ill_conditioning == True:
            #Now guess the minimal precision required to get the correct order of magnitude of the last coefficient
            #We do this by constructing "r" to low precision to get the size of its largest exponent
            r_low_prec = self._get_r_for_taylor_expansion(cusp,trunc_order,64)
            required_prec = int(round(r_low_prec[r_low_prec.degree()].abs().log10()))
            working_prec = max(digit_prec,required_prec)
            if working_prec > digit_prec:
                print("Used higher digit precision during Hauptmodul q-expansion computation: ", working_prec)
        else:
            working_prec = digit_prec

        working_bit_prec = digits_to_bits(working_prec)
        CBF = ComplexBallField(working_bit_prec)
        cusp_evaluation = self._cusp_valuations[cusp]
        cusp_width = self.G.cusp_width(cusp)
        r = self._get_r_for_taylor_expansion(cusp,trunc_order,working_bit_prec)

        n_sqrt_j_inverse = get_n_th_root_of_1_over_j(trunc_order,cusp_width).polynomial().change_ring(CBF)
        r_pos_degree = r.power_series().polynomial().change_ring(CBF)
        tmp = r_pos_degree.compose_trunc(n_sqrt_j_inverse,r.degree()+1) #Perform composition in arb because it is expensive
        P = PowerSeriesRing(CBF,n_sqrt_j_inverse.variable_name())
        j_G = CBF(cusp_evaluation) + P(tmp).O(trunc_order)
        return j_G.change_ring(CC_res)
    
    def _get_hauptmodul_q_expansion_derivative(self, j_G):
        """
        Returns 1/(2*pi*i) * d/dtau j_Gamma(tau).
        This function works for both rigorous and approximate arithmetic.
        """
        q = j_G.parent().gen()
        j_G_prime = 0
        for n in range(1,j_G.prec()): #The derivative of the q^0 term is zero
            j_G_prime += n*j_G[n]*q**n
        if j_G[-1] != 0: #j_G starts with 1/q instead of a constant term
            j_G_prime += -1/q #Derivative of q^-1
        return j_G_prime.O(j_G.prec())

    def _get_regularized_modular_form_q_expansion(self, weight, j_G, B_factored):
        """
        Returns a (non-holomorphic!) modular form of G that has no poles outside infinity and no zeros.
        This form is given by (j'_Gamma)^weight_half/B.
        """
        weight_half = weight//2
        j_G_prime = self._get_hauptmodul_q_expansion_derivative(j_G)
        num = j_G_prime**weight_half
        #It is probably best to avoid divisions of PowerSeries so we first build the denominator through multiplication
        den = 1
        for (B_factor,order) in B_factored:
            x = B_factor.parent().gen()
            factor_q_expansion = B_factor.subs(x=j_G)
            den *= factor_q_expansion**order #Working with powers of these factors should generally be faster than constructing p(x) and substituting with Horner
        return num/den
