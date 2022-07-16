from sage.libs.arb.acb cimport *
from sage.rings.complex_arb cimport *
from sage.libs.arb.acb_poly cimport *
from sage.rings.qqbar import QQbar, AlgebraicField
from sage.modular.cusps import Cusp
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.power_series_poly import PowerSeries_poly
from sage.modular.modform.j_invariant import j_invariant_qexp
from sage.rings.complex_field import ComplexField
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.real_mpfr import RealField
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.polynomial.polynomial_complex_arb import Polynomial_complex_arb
from sage.rings.complex_arb import ComplexBallField
from sage.rings.integer import Integer
from sage.misc.misc import newton_method_sizes
from sage.functions.other import ceil

from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from psage.modform.arithgroup.mysubgroup import MySubgroup

from arblib_helpers.acb_approx cimport *
from belyi.newton_genus_zero import run_newton
from belyi.number_fields import get_decimal_digit_prec
from point_matching.point_matching_arb_wrap import get_coefficients_haupt_ir_arb_wrap, digits_to_bits, bits_to_digits, get_pi_ball
from classes.fourier_expansion import FourierExpansion, to_reduced_row_echelon_form, get_trunc_orders_convergence
from classes.factored_polynomial import get_factored_polynomial_in_u_v

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
    tmp = L(j_inv.V(n).subs(q=q_n)) #Because we currently cannot work with Puiseux series in Sage.
    res = tmp.nth_root(n)
    return res

def my_n_th_root(p, n):
    """
    Given a Laurent series p defined over acbs for which the first n-1 coefficients are zero, compute the n-th root using division free newton iterations.
    The reason for implementing this ourselves is that the following example does not work in sage:
    sage: P.<x> = PowerSeriesRing(CBF)
    sage: p = P(2*x**3-16*x**4).O(10)
    sage: p.nth_root(3)
    TypeError: Cannot convert int to sage.rings.integer.Integer
    """
    cdef ComplexBall c
    if n == 1:
        return p
    elif n == 2: #Use arbs implementation which is more optimized
        return p.parent()(p.power_series().sqrt())
    else:
        CBF = p.parent().base_ring()
        c = CBF(p[n]) #Make sure we do not modify p[n]
        acb_root_ui(c.value,c.value,n,CBF.precision()) #because CBF.nth_root is not implemented...
        prec = p.prec()-n    
        r = ((1/c)*p.parent().gen()**(-1)).O(0)
        for i in newton_method_sizes(prec)[1:]:
            r = r.lift_to_precision(i-1)
            r = (r*((n+1)-p*r**n))/n
        return r.inverse()

def get_all_nth_roots_cbf(x, n):
    """
    Given a CBF x, compute all the nth roots of unity with rigorous error bounds.
    """
    cdef ComplexBall r
    CBF = x.parent()
    bit_prec = CBF.precision()
    r = CBF(x.abs())
    acb_root_ui(r.value,r.value,n,bit_prec) #because CBF.nth_root is not implemented...
    phi = x.arg()
    pi = get_pi_ball(bit_prec)
    nth_roots = [r*(((phi+2*pi*l)*CBF(0,1)/n).exp()) for l in range(n)]
    return nth_roots

def my_n_th_root_with_correct_embedding(p, n, q_coefficient):
    """
    Given a Laurent series p for which the first n-1 coefficients are zero, compute the n-th root using division free newton iterations.
    We also make sure that we choose the correct root by comparing the roots to the numerical value of the q^1 coefficient of the hauptmodul obtained from Hejhal's method.
    """
    if n == 1:
        return p
    else:
        correct_embedding = 1/q_coefficient
        if isinstance(p.parent().base_ring(),ComplexBallField):
            CBF = p.parent().base_ring()
            bit_prec = CBF.precision()
            nth_roots = get_all_nth_roots_cbf(p[n],n)
        elif isinstance(p.parent().base_ring(),AlgebraicField):
            bit_prec = 1024 #This should be enough to identify the correct root...
            nth_roots = p[n].nth_root(n,all=True)
        else:
            raise NotImplementedError("We have not considered this case yet!")
        CC = ComplexField(bit_prec)
        diffs = [(CC(nth_root)-correct_embedding).abs() for nth_root in nth_roots]
        c = nth_roots[diffs.index(min(diffs))]

        prec = p.prec()-n    
        r = ((1/c)*p.parent().gen()**(-1)).O(0)
        for i in newton_method_sizes(prec)[1:]:
            r = r.lift_to_precision(i-1)
            r = (r*((n+1)-p*r**n))/n
        return r.inverse()

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

def get_u_str(u_interior_Kv, principal_cusp_width):
    """
    Returns string of symbolic expression of u.
    """
    if principal_cusp_width == 1:
        return u_interior_Kv.__str__() #u is just a rational
    else:
        return "(" + u_interior_Kv.__str__() + ")^(1/" + str(principal_cusp_width) + ")"

class BelyiMap():
    """
    Class for working with Belyi maps that are expressed in terms of Factored_Polynomials with algebraic coefficients.
    """
    def __init__(self, G, newton_res=None, starting_digit_prec=42, target_digit_prec=50000, max_extension_field_degree=None):
        """
        Compute Belyi map from scratch. This is done by first using Hejhal's method to get approximate values at the elliptic points and cusps.
        Afterwards, Newton iterations are used to refine the solution. In the last step, the result is verified.
        """
        if G.genus() != 0:
            raise ArithmeticError("This function only works for genus zero subgroups!")
        if max_extension_field_degree == None:
            print("We are using the index as max_extension_field_degree which is not correct in general!")
            max_extension_field_degree = G.index()
        
        S = AutomorphicFormSpace(G,0)
        G = S.group() #In case the input has been a sage-subgroup, we now make sure to get a psage-subgroup
        if newton_res == None:
            (p3, p2, pc), cusp_rep_values, j_G_hejhal, v_Kw, u_interior_Kv = run_newton(S,starting_digit_prec,target_digit_prec,stop_when_coeffs_are_recognized=True,return_cusp_rep_values=True,max_extension_field_degree=max_extension_field_degree)
        else: #We have loaded the result and use this to reconstruct the belyi map
            (p3, p2, pc), cusp_rep_values, j_G_hejhal, v_Kw, u_interior_Kv = newton_res
        if len(cusp_rep_values) != 0:
            CIF = ComplexIntervalField(cusp_rep_values[0][1].parent().precision())
            cusp_rep_values = [(cusp,CIF(cusp_evaluation)) for (cusp,cusp_evaluation) in cusp_rep_values]
        self.G = G
        self.principal_cusp_width = G.cusp_width(Cusp(1,0))
        self.p3, self.p2, self.pc = p3, p2, pc
        self.p3_constructed, self.p2_constructed, self.pc_constructed = p3.construct(), p2.construct(), pc.construct()
        self.p3_u_v, self.p2_u_v, self.pc_u_v = get_factored_polynomial_in_u_v(p3,u_interior_Kv,self.principal_cusp_width), get_factored_polynomial_in_u_v(p2,u_interior_Kv,self.principal_cusp_width), get_factored_polynomial_in_u_v(pc,u_interior_Kv,self.principal_cusp_width)
        self._v_Kw, self._u_interior_Kv = v_Kw, u_interior_Kv
        self._Kv, self._Kw = u_interior_Kv.parent(), v_Kw.parent()
        self.v_QQbar = QQbar(self._Kv.gen())
        if self.principal_cusp_width == 1:
            self.u_QQbar = QQbar(self._u_interior_Kv)
        else:
            self.u_QQbar = QQbar(self._Kw.gen())
        self._j_G_hejhal = j_G_hejhal._convert_to_CC() #Note that CBFs currently cannot be stored
        self._cusp_rep_values = cusp_rep_values

        self._p2_fixed = self._get_e2_fixed_point_polynomial()
        self._p3_fixed = self._get_e3_fixed_point_polynomial()
        self._p_cusp_evaluations, self._cusp_evaluations = self._get_cusp_polynomial(cusp_rep_values)
        self.verify_polynomial_equation()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return self.p3_u_v.__str__() + " / " + self.pc_u_v.__str__()
    
    def _return_res_as_dict(self):
        """
        Returns most important results as a dictionary.
        The requirements on this function are that the output should be readable from sage only without requiring psage to be installed.
        """
        res = dict()
        res["p2_factored_raw"] = self.p2.factors
        res["p3_factored_raw"] = self.p3.factors
        res["pc_factored_raw"] = self.pc.factors
        res["p2_factored_pretty"] = self.p2_u_v.factors
        res["p3_factored_pretty"] = self.p3_u_v.factors
        res["pc_factored_pretty"] = self.pc_u_v.factors
        return res
    
    def _return_newton_res(self):
        """
        Returns parameters of this instance from which the class can be initialized again.
        This is useful because it offers a way to load belyi maps even after the code of this class has changed.
        """
        if len(self._cusp_rep_values) != 0:
            CIF = ComplexIntervalField(self._cusp_rep_values[0][1].parent().precision()) #CBFs currently cannot be stored
            cusp_rep_values_cif = [(cusp,CIF(cusp_evaluation)) for (cusp,cusp_evaluation) in self._cusp_rep_values]
        else:
            cusp_rep_values_cif = []
        return (self.p3, self.p2, self.pc), cusp_rep_values_cif, self._j_G_hejhal, self._v_Kw, self._u_interior_Kv

    def print_u(self):
        """
        Prints symbolic expression and embedding of u.
        """
        if self.principal_cusp_width == 1:
            print("u = " + self.get_u_str())
        else:
            CC = ComplexField(53)
            print("u = " + self.get_u_str())
            print("with embedding: " + CC(self.u_QQbar).__str__())

    def get_u_str(self):
        """
        Returns string of symbolic expression of u.
        """
        return get_u_str(self._u_interior_Kv,self.principal_cusp_width)
    
    def print_v(self):
        """
        Prints numberfield and embedding of v.
        """
        CC = ComplexField(53)
        print("v is a root of " + self._Kv.polynomial().__str__() + " with embedding: " + CC(self.v_QQbar).__str__())
    
    def _get_e2_fixed_point_polynomial(self):
        """
        Returns a polynomial whose roots are the elliptic fixed points of order two.
        """
        for (p,multiplicity) in self.p2.factors:
            if multiplicity == 1: #We are only interested in the fixed points
                return p
        return self._Kw(1)

    def _get_e3_fixed_point_polynomial(self):
        """
        Returns a polynomial whose roots are the elliptic fixed points of order three.
        """
        for (p,multiplicity) in self.p3.factors:
            if multiplicity == 1: #We are only interested in the fixed points
                return p
        return self._Kw(1)
    
    def _get_cusp_polynomial(self, cusp_rep_values):
        """
        Returns a polynomial whose roots are the evaluations at the cusps (with multiplicity one).
        We also return the cusp representatives and the algebraic evaluations of the hauptmodul attached to these.
        The cusp at infinity gets treated separately.
        The parameter "cusp_rep_values" contains contains floating point approximations for each cusp which we try to use
        to attach each cusp to a root of pc.
        """
        G = self.G
        Kw = self._Kw
        cusp_evaluations = dict()
        p_cusp = PolynomialRing(Kw,"x").one() #Polynomial with multiplicity one roots at the cusp evaluations
        for (p,multiplicity) in self.pc.factors:
            roots = p.roots(ring=QQbar)
            if len(roots) != p.degree():
                raise ArithmeticError("We are missing a root here!")
            for (root, order) in roots:
                if order != 1:
                    raise ArithmeticError("Something went wrong, we need distinct roots here!")
                cusp = identify_cusp_from_pc_root(root,cusp_rep_values)
                if G.cusp_width(cusp) != multiplicity:
                    raise ArithmeticError("This should not happen!")
                cusp_evaluations[cusp] = root
            p_cusp *= p
        if len(cusp_evaluations) != G.ncusps()-1: #Some cusps have multiple values attached to them which is invalid
            raise ArithmeticError("This should not happen!")

        return p_cusp, cusp_evaluations

    def verify_polynomial_equation(self):
        p3_constructed, p2_constructed, pc_constructed = self.p3_constructed, self.p2_constructed, self.pc_constructed
        if p3_constructed-p2_constructed-1728*pc_constructed != 0: #Verify result
            raise ArithmeticError("Verification of polynomial_equation failed!")
        if pc_constructed.degree()+self.principal_cusp_width != self.G.index():
            raise ArithmeticError("Wrong behavior at infinity!")
        #We also have to verify that the branching behavior is correct (i.e. that the roots of the polynomials are of correct multiplicity)
        #We have already done this for pc inside _get_cusp_polynomial but should do it for p2 and p3 here
    
    def _get_B_elliptic(self, weight):
        """
        Returns product of factors (x-elliptic_evaluation)^beta_e for all elliptic fixed points.
        beta_e is chosen to cancel the zeros of j_G'^weight_half.
        """
        weight_half = weight//2
        p2_fixed, p3_fixed = self._p2_fixed, self._p3_fixed
        beta_e2 = weight_half*(2-1)//2 #We need to divide because (x-j_G(e_2)) has a double zero
        beta_e3 = weight_half*(3-1)//3 #We need to divide because (x-j_G(e_3)) has a triple zero
        return (p2_fixed**beta_e2)*(p3_fixed**beta_e3)

    def _get_B_cusp(self, weight):
        """
        Returns product of factors (x-cusp_evaluation)^alpha_c for all cusps not equal infinity.
        alpha_c is chosen to cancel the zeros of j_G'^weight_half.
        """
        weight_half = weight//2
        p_cusp = self._p_cusp_evaluations
        alpha_c = weight_half
        return p_cusp**alpha_c
   
    def _get_B(self, weight):
        """
        Returns B_cusp*B_elliptic.
        """
        B_cusp = self._get_B_cusp(weight)
        B_elliptic = self._get_B_elliptic(weight)
        return B_cusp*B_elliptic

    def _get_p_list_cuspform(self, weight, B):
        """
        This function returns a list of polynomials p such that for each p in list, (j_G'^weight_half)*p/B corresponds to a cuspform.
        p is written as a product of factors of (x-cusp_evaluation) which prescribe zeros of given order at the cusps.
        We assume that p vanishes to degree 1 at all cusps not equal infinity.
        For the cusp at infinity we assume orders of vanishing up to the dimension of the space.
        """
        weight_half = weight//2
        p_cusp = self._p_cusp_evaluations
        x = B.parent().gen()
        cuspform_dim = self.G.dimension_cusp_forms(weight)
        if cuspform_dim == 0:
            raise ArithmeticError("The dimension of cuspforms is zero for this space!")
        B_degree = B.degree()
        p_list = [p_cusp for _ in range(cuspform_dim)]

        for n in range(1,cuspform_dim+1):
            #Now we need to get the correct order of vanishing at infinity
            p_degree = -weight_half+B_degree-n
            power = p_degree-(self.G.ncusps()-1) #Exponent of x s.t. cuspform vanishes to degree n at infty.
            if power < 0:
                raise ArithmeticError("This should not happen...")
            p_list[n-1] *= x**power
        
        return p_list
    
    def _get_p_list_modform(self, weight, B):
        """
        This function returns a list of polynomials p such that for each p in list, (j_G'^weight_half)*p/B corresponds to a holomorphic modform.
        p is written as a product of factors of (x-cusp_evaluation) which prescribe zeros of given order at the cusps.
        We assume that p is constant and non-zero at all cusps not equal infinity.
        For the cusp at infinity we assume orders of vanishing up to the dimension of the space.
        """
        weight_half = weight//2
        p_cusp = self._p_cusp_evaluations
        x = B.parent().gen()
        modform_dim = self.G.dimension_modular_forms(weight)
        if modform_dim == 0:
            raise ArithmeticError("The dimension of modforms is zero for this space!")
        B_degree = B.degree()
        p_list = [x**0 for _ in range(modform_dim)]

        for n in range(1,modform_dim+1):
            #Now we need to get the correct order of vanishing at infinity
            p_degree = -weight_half+B_degree-n+1
            power = p_degree #Exponent of x s.t. modform vanishes to degree n-1 at infty.
            if power < 0:
                raise ArithmeticError("This should not happen...")
            p_list[n-1] *= x**power
        
        return p_list

    def get_cuspforms(self, weight, trunc_orders, digit_prec=None, j_G=None):
        """
        Use Hauptmodul to return a basis of cuspforms with specified weight in reduced row-echelon form.
        If "digit_prec" is given, use approximate ball arithmetic with rigorous error bounds.

        trunc_orders can either be an int in which case only the expansion at infinity will be returned or a dict
        containing the cusps as well as the expansion orders at these cusps.
        """
        if j_G == None:
            if digit_prec == None:
                j_G = self.get_hauptmodul_q_expansion(trunc_orders)
            else:
                j_G = self.get_hauptmodul_q_expansion_approx(trunc_orders,digit_prec)
        B = self._get_B(weight)
        p_list = self._get_p_list_cuspform(weight,B)
        F = self._get_regularized_modular_form_q_expansion(weight,j_G,B) #We could re-use this for modforms of the same weight

        cuspforms = []
        base_ring = j_G.cusp_expansions[Cusp(1,0)].base_ring()
        for p in p_list:
            coeffs = list(p)
            prefactor = coeffs[-1]
            for i in range(len(coeffs)-2,-1,-1): #Horner's method
                prefactor = j_G*prefactor+coeffs[i]
            cuspform = F*prefactor
            cuspform.modform_type = "CuspForm"
            cuspforms.append(cuspform)
        cuspforms_reduced_row_echelon_form = to_reduced_row_echelon_form(cuspforms)
        if digit_prec != None:
            for i in range(len(cuspforms_reduced_row_echelon_form)):
                cuspforms_reduced_row_echelon_form[i]._set_constant_coefficients_to_zero_inplace()
        return cuspforms_reduced_row_echelon_form
    
    def get_modforms(self, weight, trunc_orders, digit_prec=None, j_G=None):
        """
        Use Hauptmodul to return a basis of modforms with specified weight in reduced row-echelon form.
        If "digit_prec" is given, use approximate ball arithmetic with rigorous error bounds.

        trunc_orders can either be an int in which case only the expansion at infinity will be returned or a dict
        containing the cusps as well as the expansion orders at these cusps.
        """
        if j_G == None:
            if digit_prec == None:
                j_G = self.get_hauptmodul_q_expansion(trunc_orders)
            else:
                j_G = self.get_hauptmodul_q_expansion_approx(trunc_orders,digit_prec)
        B = self._get_B(weight)
        p_list = self._get_p_list_modform(weight,B)
        F = self._get_regularized_modular_form_q_expansion(weight,j_G,B) #We could re-use this for cuspforms of the same weight

        modforms = []
        base_ring = j_G.cusp_expansions[Cusp(1,0)].base_ring()
        for p in p_list:
            coeffs = list(p)
            prefactor = coeffs[-1]
            for i in range(len(coeffs)-2,-1,-1): #Horner's method
                prefactor = j_G*prefactor+coeffs[i]
            modform = F*prefactor
            modform.modform_type = "ModForm"
            modforms.append(modform)
        return to_reduced_row_echelon_form(modforms)

    def get_hauptmodul_q_expansion(self, trunc_orders):
        """
        Return q-expansion of hauptmodul at all cusps using rigorous arithmetic.
        We return the result as an instance of FourierExpansion.

        trunc_orders can either be an int in which case only the expansion at infinity will be returned or a dict
        containing the cusps as well as the expansion orders at these cusps.
        """
        cusp_expansions = dict()
        if isinstance(trunc_orders,int) == True or isinstance(trunc_orders,Integer): #We are only considering the expansion at infinity
            only_principal_cusp_expansion = True
            cusp_expansions[Cusp(1,0)] = self._get_hauptmodul_q_expansion_infinity(trunc_orders)
        else:
            only_principal_cusp_expansion = False
            for cusp in self.G.cusps():
                if cusp == Cusp(1,0):
                    cusp_expansion = self._get_hauptmodul_q_expansion_infinity(trunc_orders[cusp]).change_ring(QQbar) #We need to work with QQbar because of the other cusps
                else:
                    cusp_expansion = self._get_hauptmodul_q_expansion_non_infinity(cusp,trunc_orders[cusp])
                cusp_expansions[cusp] = cusp_expansion
        return FourierExpansion(self.G,0,cusp_expansions,"Hauptmodul",self._Kw,only_principal_cusp_expansion=only_principal_cusp_expansion,Kw=self._Kw,Kv=self._Kv,u_interior_Kv=self._u_interior_Kv)
    
    def get_hauptmodul_q_expansion_approx(self, trunc_orders, digit_prec, min_digits_last_coeff=1, only_principal_cusp_expansion=True):
        """
        Return q-expansion of hauptmodul at all cusps using ball arithmetic with rigorous error bounds.
        We return the result as an instance of FourierExpansion.

        trunc_orders can either be an int in which case only the expansion at infinity will be returned or a dict
        containing the cusps as well as the expansion orders at these cusps.
        """
        cusp_expansions = dict()
        if isinstance(trunc_orders,int) == True or isinstance(trunc_orders,Integer): #We are only considering the expansion at infinity
            only_principal_cusp_expansion = True
            cusp_expansions[Cusp(1,0)] = self._get_hauptmodul_q_expansion_infinity_approx(trunc_orders,digit_prec,min_digits_last_coeff=min_digits_last_coeff)
        else:
            for cusp in self.G.cusps():
                if cusp == Cusp(1,0):
                    cusp_expansion = self._get_hauptmodul_q_expansion_infinity_approx(trunc_orders[cusp],digit_prec,min_digits_last_coeff=min_digits_last_coeff)
                else:
                    cusp_expansion = self._get_hauptmodul_q_expansion_non_infinity_approx(cusp,trunc_orders[cusp],digit_prec,min_digits_last_coeff=min_digits_last_coeff)
                cusp_expansions[cusp] = cusp_expansion
        base_ring = ComplexBallField(digits_to_bits(digit_prec))
        return FourierExpansion(self.G,0,cusp_expansions,"Hauptmodul",base_ring,only_principal_cusp_expansion=only_principal_cusp_expansion)

    def _get_hauptmodul_q_expansion_infinity(self, trunc_order):
        """
        Computes the q-expansion at the principal cusp which is normalized to j_Gamma = 1/q_N + 0 + c_1*q_N + ...
        """
        parent = self.p2_constructed[0].parent()
        L = LaurentSeriesRing(parent,"x")
        x = L.gen()
        principal_cusp_width = self.principal_cusp_width
        s = (L(self.pc_constructed).V(-1).O(trunc_order)/L(self.p3_constructed).V(-1).O(trunc_order)).power_series().nth_root(self.principal_cusp_width)
        r = s.reverse().inverse()
        n_sqrt_j_inverse = get_n_th_root_of_1_over_j(trunc_order,principal_cusp_width)
        j_G = r.subs(x=n_sqrt_j_inverse)
        return j_G

    def _get_hauptmodul_q_expansion_non_infinity(self, cusp, trunc_order):
        """
        Computes the q-expansion at a non-principal cusp which is normalized to j_Gamma = c_0 + c_1*q_N + ...
        Because these q-expansions are usually defined over a field Kw times another root of an element in Kv, which seems very tedious
        to implement, we work over QQbar which can become very slow for large examples.
        """
        q_coefficient = self._j_G_hejhal.get_cusp_expansion(cusp)[1] #Coefficient of the q^1 term of the hauptmodul which we will need to identify the correct nth root
        cusp_evaluation = QQbar(self._cusp_evaluations[cusp])
        cusp_width = self.G.cusp_width(cusp)
        working_trunc_order = trunc_order+cusp_width #Precision with which the arithmetic needs to be performed
        L = LaurentSeriesRing(QQbar,"x") #The expansions at other cusps can be defined over a different numberfield than Kw, so we have to use QQbar...
        x = L.gen()
        s_no_nth_root = (L(self.pc_constructed).subs(x=x+cusp_evaluation).O(working_trunc_order)/L(self.p3_constructed).subs(x=x+cusp_evaluation).O(working_trunc_order))
        s = my_n_th_root_with_correct_embedding(s_no_nth_root,cusp_width,q_coefficient).power_series()
        print("Warning, computing q-expansions at other cusps explicitly can be very slow because we use the QQbar type!")
        r = s.reverse()
        n_sqrt_j_inverse = get_n_th_root_of_1_over_j(working_trunc_order,cusp_width)
        j_G = QQbar(cusp_evaluation) + r.subs(x=n_sqrt_j_inverse)
        return j_G.O(trunc_order)

    # def _get_hauptmodul_q_expansion_infinity_approx_sage(self, trunc_order, digit_prec):
    #     """
    #     Computes q-expansion of hauptmodul using Sage implementations.
    #     These are significantly slower than the ones provided by Arb and it is hence advised to use "get_hauptmodul_approx_q_expansion" instead.
    #     """
    #     bit_prec = digits_to_bits(digit_prec)
    #     CC = ComplexField(bit_prec)
    #     L = LaurentSeriesRing(CC,"x")
    #     x = L.gen()
    #     principal_cusp_width = self.principal_cusp_width
    #     s = (L(self.pc_constructed).subs({x:1/x}).O(trunc_order)/L(self.p3_constructed).subs({x:1/x}).O(trunc_order)).power_series().nth_root(self.principal_cusp_width)
    #     r = s.reverse().inverse()
    #     n_sqrt_j_inverse = get_n_th_root_of_1_over_j(trunc_order,principal_cusp_width)
    #     j_G = r.subs({x:n_sqrt_j_inverse})
    #     return j_G #Note that some of the higher coefficients might be wrong due to rounding errors

    def _get_r_for_laurent_expansion(self, trunc_order, bit_prec):
        """
        Get the reversed series of the reciprocal of the Belyi map expanded in 1/x.
        We need to compute this for "get_hauptmodul_q_expansion_infinity_approx".
        """
        CBF = ComplexBallField(bit_prec)
        L = LaurentSeriesRing(CBF,"x")
        x = L.gen()
        s = my_n_th_root(L(self.pc_constructed).V(-1).O(trunc_order)/L(self.p3_constructed).V(-1).O(trunc_order),self.principal_cusp_width)
        s_prec = s.prec() #Exponent of the O-term
        s_arb_reverted = s.power_series().polynomial().revert_series(s_prec) #Perform the reversion in arb because it is expensive
        r = L(s_arb_reverted).O(s_prec).inverse()
        return r

    def _get_hauptmodul_q_expansion_infinity_approx(self, trunc_order, digit_prec, min_digits_last_coeff=1, guess_minimal_working_precision=True):
        """
        Compute approximation of the coefficients of the q-expansion of the hauptmodul truncated to "trunc_order" and with "digit_prec"
        working precision. We make use of fast Arb implementations for performance and rigorous error bounds.
        Because of the large coefficients involved, the arithmetic might become ill-conditioned. 
        We therefore try to guess the minimal working precision required and return a result for which the last coefficient
        has at least min_digits_last_coeff digits precision.
        """
        principal_cusp_width = self.principal_cusp_width
        working_trunc_order = trunc_order
        working_prec = digit_prec
        if guess_minimal_working_precision == True:
            #Now guess the minimal precision required to get the correct order of magnitude of the last coefficient
            #We do this by constructing "r" to low precision to get the size of its largest exponent
            r_low_prec = self._get_r_for_laurent_expansion(working_trunc_order,64)
            CC = ComplexField(64)
            required_prec = int((trunc_order/4+CC(r_low_prec[r_low_prec.degree()]).abs().log10()).round()) #Because log10 is not defined for arb...
            working_prec = max(digit_prec,required_prec)
            if working_prec > digit_prec:
                print("Used higher digit precision during Hauptmodul q-expansion computation: ", working_prec)
            if r_low_prec.prec() < trunc_order: #It somehow sometimes happens for large orders that r is only computed to fewer terms
                working_trunc_order += trunc_order-r_low_prec.prec()+5
        n_sqrt_j_inverse = get_n_th_root_of_1_over_j(working_trunc_order,principal_cusp_width)
        working_bit_prec = digits_to_bits(working_prec)
        while True:
            CBF = ComplexBallField(working_bit_prec)
            n_sqrt_j_inverse_CBF = n_sqrt_j_inverse.polynomial().change_ring(CBF)
            r = self._get_r_for_laurent_expansion(working_trunc_order,working_bit_prec)
            r_pos_degree = r[0:r.degree()+1].power_series().polynomial() #We treat the 1/x term later because arb only supports polynomials
            tmp = r_pos_degree.compose_trunc(n_sqrt_j_inverse_CBF,working_trunc_order) #Perform composition in arb because it is expensive
            P = PowerSeriesRing(CBF,n_sqrt_j_inverse_CBF.variable_name())
            one_over_x_term = n_sqrt_j_inverse.inverse() #Treat 1/x term separately because it is not supported by arb
            j_G = (one_over_x_term + P(tmp)).O(trunc_order)
            if get_decimal_digit_prec(j_G[j_G.degree()].rad()) < min_digits_last_coeff: #Result is not correct up to desired accuracy, therefore we need redo the computations with higher precision
                working_bit_prec = int(round(1.2*working_bit_prec)) #Obviously we could try to come up with a more sophisticated choice
            else:
                return j_G

    def _get_r_for_taylor_expansion(self, cusp, trunc_order, q_coefficient, bit_prec):
        """
        Get the reversed series of the reciprocal of the Belyi map expanded in x.
        We need to compute this for "get_hauptmodul_q_expansion_non_infinity_approx".
        """
        cusp_evaluation = self._cusp_evaluations[cusp]
        cusp_width = self.G.cusp_width(cusp)
        CBF = ComplexBallField(bit_prec)
        L = LaurentSeriesRing(CBF,"x")
        pc_shifted_QQbar, p3_shifted_QQbar = self.pc_constructed.change_ring(QQbar), self.p3_constructed.change_ring(QQbar)
        x = pc_shifted_QQbar.parent().gen()
        #We need to perform this shift in QQbar because with CBF we might get empty error balls
        #In this case, performing the shift in QQbar also seems feasible performance wise
        pc_shifted_QQbar, p3_shifted_QQbar = pc_shifted_QQbar.subs({x:x+cusp_evaluation}), p3_shifted_QQbar.subs({x:x+cusp_evaluation})
        
        #It is very important that the leading order terms of pc_shifted are truely zero, otherwise we can get very large rounding errors
        #We therefore set the coefficients that are effectively zero to true zeros
        pc_shifted_coeffs = []
        for i in range(pc_shifted_QQbar.degree()+1):
            if i < cusp_width: #It is an easy calculation to check that these are always zero
                pc_shifted_coeffs.append(CBF(0))
            else:
                if i == cusp_width: #If imag or real are empty then this can cause problems too...
                    if pc_shifted_QQbar[i].real() == 0:
                        pc_shifted_coeffs.append(CBF(0,pc_shifted_QQbar[i].imag()))
                    elif pc_shifted_QQbar[i].imag() == 0:
                        pc_shifted_coeffs.append(CBF(pc_shifted_QQbar[i].real(),0))
                    else:
                        pc_shifted_coeffs.append(CBF(pc_shifted_QQbar[i]))
                else:
                    pc_shifted_coeffs.append(CBF(pc_shifted_QQbar[i]))
        pc_shifted = L(pc_shifted_coeffs).O(trunc_order)
        p3_shifted_coeffs = []
        for i in range(p3_shifted_QQbar.degree()+1):
            if i == 0: #Again, zeros can cause problems...
                if p3_shifted_QQbar[i].real() == 0:
                        p3_shifted_coeffs.append(CBF(0,p3_shifted_QQbar[i].imag()))
                elif p3_shifted_QQbar[i].imag() == 0:
                    p3_shifted_coeffs.append(CBF(p3_shifted_QQbar[i].real(),0))
                else:
                    p3_shifted_coeffs.append(CBF(p3_shifted_QQbar[i]))
            else:
                p3_shifted_coeffs.append(CBF(p3_shifted_QQbar[i]))
        p3_shifted = L(p3_shifted_coeffs).O(trunc_order)

        s_no_nth_root = pc_shifted/p3_shifted
        if s_no_nth_root[s_no_nth_root.valuation()].real().accuracy() <= 0: #Again, we have to look for empty balls...
            s_no_nth_root_list = list(s_no_nth_root.power_series().polynomial())
            s_no_nth_root_list[s_no_nth_root.valuation()] = CBF(0,s_no_nth_root[s_no_nth_root.valuation()].imag())
            s_no_nth_root = L(s_no_nth_root_list).O(s_no_nth_root.prec())
        elif s_no_nth_root[s_no_nth_root.valuation()].imag().accuracy() <= 0:
            s_no_nth_root_list = list(s_no_nth_root.power_series().polynomial())
            s_no_nth_root_list[s_no_nth_root.valuation()] = CBF(s_no_nth_root[s_no_nth_root.valuation()].real(),0)
            s_no_nth_root = L(s_no_nth_root_list).O(s_no_nth_root.prec())

        s = my_n_th_root_with_correct_embedding(s_no_nth_root,cusp_width,q_coefficient)
        s_prec = s.prec() #Exponent of the O-term
        s_arb_reverted = s.power_series().polynomial().revert_series(s_prec) #Perform the reversion in arb because it is expensive
        s_arb_reverted_prec = s_arb_reverted.prec() #Exponent of the O-term
        r = L(s_arb_reverted).O(s_arb_reverted_prec)
        return r

    def _get_hauptmodul_q_expansion_non_infinity_approx(self, cusp, trunc_order, digit_prec, min_digits_last_coeff=1, guess_minimal_working_precision=True):
        """
        Compute approximation of the coefficients of the q-expansion of the hauptmodul at a non-principal cusp truncated to "trunc_order" 
        and with "digit_prec" working precision. We make use of fast Arb implementations for performance and rigorous error bounds.
        Because of the large coefficients involved, the arithmetic might become ill-conditioned. 
        We therefore try to guess the minimal working precision required and return a result for which the last coefficient
        has at least min_digits_last_coeff digits precision.
        """
        q_coefficient = self._j_G_hejhal.get_cusp_expansion(cusp)[1] #Coefficient of the q^1 term of the hauptmodul which we will need to identify the correct nth root
        cusp_width = self.G.cusp_width(cusp)
        working_trunc_order = trunc_order+cusp_width #Precision with which the arithmetic needs to be performed
        working_prec = digit_prec
        if guess_minimal_working_precision == True:
            #Now guess the minimal precision required to get the correct order of magnitude of the last coefficient
            #We do this by constructing "r" to low precision to get the size of its largest exponent
            r_low_prec = self._get_r_for_taylor_expansion(cusp,working_trunc_order,q_coefficient,64)
            CC = ComplexField(64)
            required_prec = int((trunc_order/4+CC(r_low_prec[r_low_prec.degree()]).abs().log10()).round()) #Because log10 is not defined for arb...
            working_prec = max(digit_prec,required_prec)
            if working_prec > digit_prec:
                print("Used higher digit precision during Hauptmodul q-expansion computation: ", working_prec)
            if r_low_prec.prec() < trunc_order: #It somehow sometimes happens for large orders that r is only computed to fewer terms
                working_trunc_order += trunc_order-r_low_prec.prec()+5
        n_sqrt_j_inverse = get_n_th_root_of_1_over_j(working_trunc_order,cusp_width)
        working_bit_prec = digits_to_bits(working_prec)
        cusp_evaluation = self._cusp_evaluations[cusp]
        while True:
            CBF = ComplexBallField(working_bit_prec)
            n_sqrt_j_inverse_CBF = n_sqrt_j_inverse.polynomial().change_ring(CBF)
            r = self._get_r_for_taylor_expansion(cusp,working_trunc_order,q_coefficient,working_bit_prec)
            r_pos_degree = r.power_series().polynomial()
            tmp = r_pos_degree.compose_trunc(n_sqrt_j_inverse_CBF,r.degree()+1) #Perform composition in arb because it is expensive
            P = PowerSeriesRing(CBF,n_sqrt_j_inverse_CBF.variable_name())
            j_G = (CBF(cusp_evaluation) + P(tmp)).O(trunc_order)
            if get_decimal_digit_prec(j_G[j_G.degree()].rad()) < min_digits_last_coeff: #Result is not correct up to desired accuracy, therefore we need redo the computations with higher precision
                working_bit_prec = int(round(1.2*working_bit_prec)) #Obviously we could try to come up with a more sophisticated choice
            else:
                return j_G

    def _get_hauptmodul_q_expansion_derivative(self, j_G):
        """
        Returns 1/(2*pi*i) * d/dtau j_G(tau) where j_G is an instance of "FourierExpansion".
        This function works for both rigorous and floating arithmetic.
        """
        cusp_expansions = dict()
        for cusp in j_G.cusp_expansions.keys():
            cusp_expansion = j_G.get_cusp_expansion(cusp)
            q = cusp_expansion.parent().gen()
            cusp_expansion_prime = 0
            for n in range(1,cusp_expansion.prec()): #The derivative of the q^0 term is zero
                cusp_expansion_prime += n*cusp_expansion[n]*q**n
            if cusp_expansion[-1] != 0: #cusp_expansion starts with 1/q instead of a constant term
                cusp_expansion_prime += -1/q #Derivative of q^-1
            cusp_expansions[cusp] = cusp_expansion_prime.O(cusp_expansion.prec())
        return FourierExpansion(j_G.G,2,cusp_expansions,"ModForm",j_G.base_ring,
                only_principal_cusp_expansion=j_G.only_principal_cusp_expansion,Kw=j_G._Kw,Kv=j_G._Kv,u_interior_Kv=j_G._u_interior_Kv)

    def _get_regularized_modular_form_q_expansion(self, weight, j_G, B):
        """
        Returns a (non-holomorphic!) modular form of G that has no poles outside infinity and no zeros.
        This form is given by (j'_Gamma)^weight_half/B.

        Note:
        Inside this function we have to compute quotients of power series which can be ill-conditioned.
        If the coefficients have too low precision, we therefore mutate j_G to a higher precision and redo the computations until
        a sufficient precision is reached.
        """
        weight_half = weight//2
        restart_computation = True
        while restart_computation == True:
            restart_computation = False
            j_G_prime = self._get_hauptmodul_q_expansion_derivative(j_G)
            num = j_G_prime**weight_half
            base_ring = j_G.base_ring
            coeffs = list(B)
            den = coeffs[-1]
            for i in range(len(coeffs)-2,-1,-1): #Horner's method
                den = j_G*den+coeffs[i]
            if isinstance(base_ring,ComplexBallField) == True:
                #This part is really ugly:
                #Because arb does not detect the true zeros and returns empty error bounds, we would get NaN's when dividing with "den"
                #We therefore need to set the first weight_half coefficients at the cusps outside infinity to zero
                #(It is obvious that these need to be zero because their constant term is equivalent to one of the cusp evaluations)
                for c in den.cusp_expansions.keys():
                    if c != Cusp(1,0):
                        P = den.cusp_expansions[c].parent()
                        den_cusp_expansions_list = list(den.cusp_expansions[c]) #We need to work with lists because power series are immutable
                        for i in range(weight_half):
                            den_cusp_expansions_list[i] = 0 #Set these coefficients to be truely zero
                        den.cusp_expansions[c] = P(den_cusp_expansions_list).O(den.cusp_expansions[c].prec()) #update expression in class instance
            F = num/den

            if isinstance(base_ring,ComplexBallField) == True: #We need to check that the last coefficient has enough precision for the modforms to be valid
                for c in F.cusp_expansions.keys():
                    F_cusp_expansion = F.cusp_expansions[c]
                    if get_decimal_digit_prec(F_cusp_expansion[F_cusp_expansion.degree()].rad()) < 0.05*F_cusp_expansion.degree(): #The construction of modforms or cuspforms might produce empty error balls
                        restart_computation = True
                        current_digit_prec = bits_to_digits(F_cusp_expansion.base_ring().precision())
                        j_G_cusp_expansion = j_G.cusp_expansions[c]
                        new_digit_prec = int(round(1.2*current_digit_prec)) #Obviously we could try to come up with a more sophisticated choice
                        if c == Cusp(1,0):
                            j_G.cusp_expansions[c] = self._get_hauptmodul_q_expansion_infinity_approx(j_G_cusp_expansion.prec(),new_digit_prec,guess_minimal_working_precision=False)
                        else:
                            j_G.cusp_expansions[c] = self._get_hauptmodul_q_expansion_non_infinity_approx(c,j_G_cusp_expansion.prec(),new_digit_prec,guess_minimal_working_precision=False)
        return F
    
    def _get_trunc_orders_equal_sized(self, trunc_order):
        """
        Returns a dictionary with cusps of self.G as keys and trunc_order as values.
        """
        trunc_orders = dict()
        for c in self.G.cusps():
            trunc_orders[c] = trunc_order
        return trunc_orders
    
    def _get_trunc_orders_convergence(self, weight, digit_prec):
        """
        Returns a dictionary with cusps of self.G as keys and trunc_orders as values.
        We choose the trunc_order for each cusp in a way that the expansion of f(z) (approximately) converges inside the fundamental domain.
        Where f(z) is a modular form of specified weight k for which we assume that the coefficients grow like O(n^k)
        """
        return get_trunc_orders_convergence(self.G,weight,digit_prec)
