from sage.rings.qqbar import QQbar
from sage.modular.cusps import Cusp
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.power_series_poly import PowerSeries_poly
from sage.modular.modform.j_invariant import j_invariant_qexp
from sage.rings.complex_field import ComplexField
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.real_mpfr import RealField
from sage.matrix.matrix_space import MatrixSpace

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

def get_j_Gamma_hejhal(S,digit_prec,trunc_order=None):
    """
    Get q-expansion of j_Gamma (we currently use Hejhal's method for this).
    We are currently only considering the cusp at infinity.
    """
    G = S.group()
    if G.cusp_width(Cusp(1,0)) != 1: #Here we would need to work with q_n
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

def get_n_th_root_of_1_over_j(trunc_order,n):
    """
    Returns (1/j)^(1/n) up to trunc_order terms.
    The result is a power series in q_n where q_n = exp(2*pi*I/n).
    """
    j = j_invariant_qexp(trunc_order)
    j_inv = j.inverse()
    var_name = "q_" + str(n)
    L = PowerSeriesRing(j[0].parent(),var_name)
    q_n = L.gen()
    #.subs() seems quite slow, maybe try to write faster cython code
    tmp = L(j_inv.power_series().subs(q=q_n**n)) #Because we currently cannot work with Puiseux series in Sage.
    res = tmp.nth_root(n)
    return res

def get_approx_n_th_root_of_1_over_j(trunc_order,n,digit_prec):
    """
    Returns (1/j)^(1/n) up to trunc_order terms with digit_prec precision.
    The result is a power series in q_n where q_n = exp(2*pi*I/n).
    """
    bit_prec = digits_to_bits(digit_prec)
    RR = RealField(bit_prec)
    j = j_invariant_qexp(trunc_order)
    j_inv = j.inverse() #This seems to be in ZZ instead of QQ and is thus quite cheap
    var_name = "q_" + str(n)
    L = PowerSeriesRing(RR,var_name)
    q_n = L.gen()
    tmp = L(j_inv.power_series().polynomial().subs(q=q_n**n)) #Because we currently cannot work with Puiseux series in Sage.
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
            print("We are using the index as max_extension_field_degree which is not correct in general!")
            max_extension_field_degree = G.index()
        
        (p3, p2, pc) = run_newton(S, starting_digit_prec, target_digit_prec, stop_when_coeffs_are_recognized=True)
        self.G = G
        self.p3, self.p2, self.pc = p3, p2, pc
        self.p3_constructed, self.p2_constructed, self.pc_constructed = p3.construct(), p2.construct(), pc.construct()
        self.princial_cusp_width = G.cusp_width(Cusp(1,0))

        self._e2_valuations = self._get_e2_fixed_point_valuations()
        self._e3_valuations = self._get_e3_fixed_point_valuations()
        self._cusp_valuations = self._get_cusp_valuations()
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
    
    def _get_cusp_valuations(self):
        """
        Returns algebraic valuations of hauptmodul at cusps.
        The cusp at infinity gets treated separately.
        """
        cusp_valuations = []
        for (p,multiplicity) in self.pc.factors:
            roots = p.roots(ring=QQbar)
            for (root, order) in roots:
                if order != 1:
                    raise ArithmeticError("Something went wrong, we need distinct roots here!")
                cusp_valuations.append(root)
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
        for cusp_valuation in cusp_valuations:
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

    def _get_p_list_cuspform(self, weight):
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
        B_factored = self._get_B_factored(weight) #!!! Maybe we want to precompute this
        B_factored_degree = get_B_factored_degree(B_factored)
        p_list = [[] for _ in range(cuspform_dim)]

        for n in range(1,cuspform_dim+1):
            for cusp_valuation in cusp_valuations:
                p_list[n-1].append([x-cusp_valuation,1]) #Prescribe zeros of order one at all cusps != infty
            #Now we need to get the correct order of vanishing at infinity
            p_degree = -weight_half+B_factored_degree-n
            power = p_degree-len(cusp_valuations) #Exponent of x s.t. cuspform vanishes to degree n at infty.
            if power < 0:
                raise ArithmeticError("This should not happen...")
            p_list[n-1].append([x,power])
        
        return p_list
    
    def _get_p_list_modform(self, weight):
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
        B_factored = self._get_B_factored(weight) #!!! Maybe we want to precompute this
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
    
    def get_cuspforms(self, weight, trunc_order):
        """
        Use Hauptmodul to return a basis of cuspforms with specified weight in reduced row-echelon form.
        """
        F = self._get_regularized_modular_form_q_expansion(weight, trunc_order)
        p_list = self._get_p_list_cuspform(weight)
        j_G = self.get_hauptmodul_q_expansion(trunc_order) #Again, we could precompute this
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
    
    def get_modforms(self, weight, trunc_order):
        """
        Use Hauptmodul to return a basis of modforms with specified weight in reduced row-echelon form.
        """
        F = self._get_regularized_modular_form_q_expansion(weight, trunc_order)
        p_list = self._get_p_list_modform(weight)
        j_G = self.get_hauptmodul_q_expansion(trunc_order) #Again, we could precompute this
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
        parent = self.p2_constructed[0].parent()
        L = LaurentSeriesRing(parent,"x")
        x = L.gen()
        princial_cusp_width = self.princial_cusp_width
        s = (L(self.pc_constructed).subs(x=1/x).O(trunc_order)/L(self.p3_constructed).subs(x=1/x).O(trunc_order)).power_series().nth_root(self.princial_cusp_width)
        r = s.reverse().inverse()
        n_sqrt_j_inverse = get_n_th_root_of_1_over_j(trunc_order,princial_cusp_width)
        j_G = r.subs(x=n_sqrt_j_inverse)
        return j_G
    
    def get_hauptmodul_approx_q_expansion(self, trunc_order, digit_prec):
        print("This function is not working correctly yet")
        #The .subs is not working correctly yet!
    
    def get_hauptmodul_q_expansion_derivative(self, trunc_order):
        """
        Returns 1/(2*pi*i) * d/dtau j_Gamma(tau).
        """
        j_G = self.get_hauptmodul_q_expansion(trunc_order)
        q = j_G.parent().gen()
        j_G_prime = -1/q #Derivative of q^-1
        for n in range(1,j_G.degree()+1): #The derivative of the q^0 term is zero
            j_G_prime += n*j_G[n]*q**n
        return j_G_prime.O(j_G.degree()+1)
    
    def _get_regularized_modular_form_q_expansion(self, weight, trunc_order):
        """
        Returns a (non-holomorphic!) modular form of G that has no poles outside infinity and no zeros.
        This form is given by (j'_Gamma)^weight_half/B.
        """
        #Note that many of these expressions could be precomputed
        weight_half = weight//2
        j_Gamma = self.get_hauptmodul_q_expansion(trunc_order)
        num = self.get_hauptmodul_q_expansion_derivative(trunc_order)**weight_half
        B_factored = self._get_B_factored(weight)
        #It is probably best to avoid divisions of PowerSeries so we first build the denominator through multiplication
        den = 1
        for (B_factor,order) in B_factored:
            x = B_factor.parent().gen()
            factor_q_expansion = B_factor.subs(x=j_Gamma)
            den *= factor_q_expansion**order #Working with powers of these factors should generally be faster than constructing p(x) and substituting with Horner
        return num/den
    
    def get_hauptmodul_approx_q_expansion(self, trunc_order, digit_prec):
        raise NotImplementedError("Not implemented yet")