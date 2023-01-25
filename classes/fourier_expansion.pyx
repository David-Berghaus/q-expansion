# ****************************************************************************
#       Copyright (C) 2022 David Berghaus <berghaus@th.physik.uni-bonn.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import copy, deepcopy

from sage.misc.banner import version_dict
vers_dict = version_dict()
from sage.modular.cusps import Cusp
if vers_dict['major'] == 9 and vers_dict['minor'] == 2: #We still need to support sage 9.2 for now
    from sage.rings.complex_field import ComplexField
else:
    from sage.rings.complex_mpfr import ComplexField
from sage.rings.complex_arb import ComplexBallField
from sage.rings.qqbar import AlgebraicField
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.number_field.number_field import NumberField_generic
from sage.functions.other import ceil

from psage.modform.maass.automorphic_forms_alg import get_M_for_holom

from point_matching.point_matching_arb_wrap import (
    get_coefficients_cuspform_ir_arb_wrap, get_coefficients_modform_ir_arb_wrap, get_coefficients_haupt_ir_arb_wrap,
    digits_to_bits, bits_to_digits, _get_normalization_cuspforms, 
    _get_normalization_modforms, get_cuspform_basis_ir_arb_wrap, get_modform_basis_ir_arb_wrap
)
from belyi.expression_in_u_and_v import factor_q_expansion_into_u_v, convert_from_Kv_to_Kw
from belyi.number_fields import get_decimal_digit_prec, to_K, is_effectively_zero

def get_cuspform_q_expansion_approx(S, digit_prec, Y=0, M_0=0, label=0, c_vec=None, prec_loss=None, use_FFT=True, use_splitting=True, use_scipy_lu=True):
    """
    Computes q-expansion of cuspform numerically and returns result as instance of "FourierExpansion".
    """
    starting_order = 1
    normalization = _get_normalization_cuspforms(S,label=label)
    if c_vec == None: #We compute c_vec from scratch
        c_vec, M_0 = get_coefficients_cuspform_ir_arb_wrap(S,digit_prec,Y=Y,M_0=M_0,return_M=True,use_FFT=use_FFT,use_splitting=use_splitting,label=label,prec_loss=prec_loss,use_scipy_lu=use_scipy_lu)
    else: #We construct ApproxModForm from previously computed solution
        if M_0 == 0:
            raise ArithmeticError("Cannot construct FourierExpansion from c_vec without specifying M_0!")
    bit_prec = digits_to_bits(digit_prec)
    c_vec_mcbd = c_vec._get_mcbd(bit_prec)
    cusp_expansions = c_vec_to_cusp_expansions(c_vec_mcbd,S,starting_order,normalization,M_0)
    base_ring = cusp_expansions[Cusp(1,0)].base_ring()
    return FourierExpansion(S.group(),S.weight(),cusp_expansions,"CuspForm",base_ring)

def get_modform_q_expansion_approx(S, digit_prec, Y=0, M_0=0, label=0, c_vec=None, prec_loss=None, use_FFT=True, use_splitting=True, use_scipy_lu=True):
    """
    Computes q-expansion of modform numerically and returns result as instance of "FourierExpansion".
    """
    starting_order = 0
    normalization = _get_normalization_modforms(S,label=label)
    if c_vec == None: #We compute c_vec from scratch
        c_vec, M_0 = get_coefficients_modform_ir_arb_wrap(S,digit_prec,Y=Y,M_0=M_0,return_M=True,use_FFT=use_FFT,use_splitting=use_splitting,label=label,prec_loss=prec_loss,use_scipy_lu=use_scipy_lu)
    else: #We construct ApproxModForm from previously computed solution
        if M_0 == 0:
            raise ArithmeticError("Cannot construct FourierExpansion from c_vec without specifying M_0!")
    bit_prec = digits_to_bits(digit_prec)
    c_vec_mcbd = c_vec._get_mcbd(bit_prec)
    cusp_expansions = c_vec_to_cusp_expansions(c_vec_mcbd,S,starting_order,normalization,M_0)
    base_ring = cusp_expansions[Cusp(1,0)].base_ring()
    return FourierExpansion(S.group(),S.weight(),cusp_expansions,"ModForm",base_ring)

def get_hauptmodul_q_expansion_approx(S, digit_prec, Y=0, M_0=0, c_vec=None, prec_loss=None, use_FFT=True, use_splitting=True, use_scipy_lu=True):
    """
    Computes q-expansion of hauptmodul numerically and returns result as instance of "FourierExpansion".
    """
    if c_vec == None: #We compute c_vec from scratch
        c_vec, M_0 = get_coefficients_haupt_ir_arb_wrap(S,digit_prec,only_principal_expansion=False,Y=Y,M_0=M_0,return_M=True,use_FFT=use_FFT,use_splitting=use_splitting,prec_loss=prec_loss,use_scipy_lu=use_scipy_lu)
    else: #We construct ApproxModForm from previously computed solution
        if M_0 == 0:
            raise ArithmeticError("Cannot construct FourierExpansion from c_vec without specifying M_0!")
    bit_prec = digits_to_bits(digit_prec)
    c_vec_mcbd = c_vec._get_mcbd(bit_prec)
    cusp_expansions = c_vec_to_cusp_expansions_hauptmodul(c_vec_mcbd,S,M_0)
    base_ring = cusp_expansions[Cusp(1,0)].base_ring()
    return FourierExpansion(S.group(),S.weight(),cusp_expansions,"Hauptmodul",base_ring)

def get_cuspform_basis_approx(S,digit_prec,Y=0,M_0=0,labels=None,prec_loss=None):
    """
    Computes a basis of cuspforms numerically (in reduced row echelon form) and returns results as instances of "FourierExpansion".
    This function is more efficient than iterating over "get_cuspform_q_expansion_approx".
    """
    c_vecs, M_0, labels = get_cuspform_basis_ir_arb_wrap(S,digit_prec,Y=Y,M_0=M_0,return_M_and_labels=True,labels=labels,prec_loss=prec_loss)
    starting_order = 1
    bit_prec = digits_to_bits(digit_prec)
    basis = []
    for i in range(len(c_vecs)):
        normalization = _get_normalization_cuspforms(S,label=labels[i])
        c_vec_mcbd = c_vecs[i]._get_mcbd(bit_prec)
        cusp_expansions = c_vec_to_cusp_expansions(c_vec_mcbd,S,starting_order,normalization,M_0)
        base_ring = cusp_expansions[Cusp(1,0)].base_ring()
        basis.append(FourierExpansion(S.group(),S.weight(),cusp_expansions,"CuspForm",base_ring))
    return basis

def get_modform_basis_approx(S,digit_prec,Y=0,M_0=0,labels=None,prec_loss=None):
    """
    Computes a basis of modforms numerically (in reduced row echelon form) and returns results as instances of "FourierExpansion".
    This function is more efficient than iterating over "get_modform_q_expansion_approx".
    """
    c_vecs, M_0, labels = get_modform_basis_ir_arb_wrap(S,digit_prec,Y=Y,M_0=M_0,return_M_and_labels=True,labels=labels,prec_loss=prec_loss)
    starting_order = 0
    bit_prec = digits_to_bits(digit_prec)
    basis = []
    for i in range(len(c_vecs)):
        normalization = _get_normalization_modforms(S,label=labels[i])
        c_vec_mcbd = c_vecs[i]._get_mcbd(bit_prec)
        cusp_expansions = c_vec_to_cusp_expansions(c_vec_mcbd,S,starting_order,normalization,M_0)
        base_ring = cusp_expansions[Cusp(1,0)].base_ring()
        basis.append(FourierExpansion(S.group(),S.weight(),cusp_expansions,"ModForm",base_ring))
    return basis

def get_trunc_orders_convergence(G, weight, digit_prec):
        """
        Returns a dictionary with cusps of self.G as keys and trunc_orders as values.
        We choose the trunc_order for each cusp in a way that the expansion of f(z) (approximately) converges inside the fundamental domain.
        Where f(z) is a modular form of specified weight k for which we assume that the coefficients grow like O(n^k)
        """
        trunc_orders = dict()
        Y_0 = 0.866025403784439 #height of fundamental domain of SL2Z
        for c in G.cusps():
            Y = Y_0/G.cusp_width(c)
            trunc_order = ceil(get_M_for_holom(Y,2*weight,digit_prec))
            trunc_orders[c] = trunc_order
        return trunc_orders

def c_vec_to_cusp_expansions(c_vec_mcbd, S, starting_order, normalization, M_0):
    """
    Given a c_vec of coefficients as a acb-matrix in Sage, convert expression to a dictionary with cusps as keys and Laurent series as values.
    """
    cusp_expansions = dict()
    weight = S.weight()
    G = S.group()
    bit_prec = c_vec_mcbd[0][0].parent().precision()
    CF = ComplexField(bit_prec)
    for ci in G._cusps:
        cii = G._cusps.index(ci)
        cusp_width = G.cusp_width(ci)
        var_string = 'q_' + str(cusp_width)
        R = PowerSeriesRing(CF, var_string)
        q = R.gen()
        f = 0
        normalization_len = len(normalization[cii])
        for i in range(normalization_len): #Set normalized coefficients if there are any
            if normalization[cii][i] != 0:
                f += normalization[cii][i]*q**(starting_order+i)
        for i in range(M_0):
            f += CF(c_vec_mcbd[i+cii*M_0][0])*q**(i+starting_order+normalization_len)
        cusp_expansions[Cusp(ci)] = f.O(M_0)
    return cusp_expansions

def c_vec_to_cusp_expansions_hauptmodul(c_vec_mcbd, S, M_0):
    """
    Given a c_vec of coefficients as a acb-matrix in Sage, convert expression to a dictionary with cusps as keys and Laurent series as values.
    This function is different to "c_vec_to_cusp_expansions" because the q-expansion of the hauptmodul starts with 1/q at the principal cusp.
    """
    cusp_expansions = dict()
    weight = S.weight()
    G = S.group()
    bit_prec = c_vec_mcbd[0][0].parent().precision()
    CF = ComplexField(bit_prec)
    for ci in G._cusps:
        cii = G._cusps.index(ci)
        cusp_width = G.cusp_width(ci)
        var_string = 'q_' + str(cusp_width)
        R = LaurentSeriesRing(CF, var_string)
        q = R.gen()
        f = 0
        if cii == 0: #The normalization here is different
            f += 1/q
            for i in range(1,M_0):
                f += CF(c_vec_mcbd[i+cii*M_0-1][0])*q**(i)
        else:
            for i in range(M_0):
                f += CF(c_vec_mcbd[i+cii*M_0][0])*q**(i)
        cusp_expansions[Cusp(ci)] = f.O(M_0)
    return cusp_expansions

def recognize_cusp_expansion_using_u(cusp_expansion, weight, G, max_rigorous_trunc_order, modform_type, label, Kv, u, v_Kw, u_interior_Kv, estimated_bit_prec=None):
    """
    Given a cusp expansion as a (floating point) power series, use u and the LLL algorithm to try to recognize coefficients.
    This function returns FourierExpansion instance of the recognized cusp_expansion truncated to the order up to which the coefficients have been recognized.
    """
    if estimated_bit_prec == None: #We have not specified the precision so we use the full working precision
        bit_prec = cusp_expansion[0].parent().precision()
    else:
        bit_prec = estimated_bit_prec
    if modform_type == "CuspForm":
        starting_order = 1
        dim = G.dimension_cusp_forms(weight)
    elif modform_type == "ModForm":
        starting_order = 0
        dim = G.dimension_modular_forms(weight)
    else:
        raise NotImplementedError("We have only implemented this functionalty for CuspForm and ModForm yet!")
    CC = ComplexField(bit_prec)
    Kw = v_Kw.parent()
    coeff_list_recognized = list()
    for i in range(min(max_rigorous_trunc_order+1,cusp_expansion.degree()+1)):
        if i <= label+starting_order:
            u_pow = 0
        else:
            u_pow = i-(label+starting_order)
        expression_to_recognize = cusp_expansion[i]/(CC(u)**u_pow)
        if expression_to_recognize.is_one() == True:
            recognized_expression = Kv(1)
        elif is_effectively_zero(expression_to_recognize,int(0.8*bits_to_digits(bit_prec))) == True: #This should only detect true zeros while also accounting for precision loss
            recognized_expression = Kv(0)
        else:
            recognized_expression = to_K(expression_to_recognize,Kv)
        if recognized_expression == None: #Stop because we have been unable to recognize coeff
            break
        recognized_expression_Kw = convert_from_Kv_to_Kw(recognized_expression, v_Kw)
        algebraic_expression = recognized_expression_Kw*(u**u_pow)
        coeff_list_recognized.append(algebraic_expression)
    P = PowerSeriesRing(Kw,cusp_expansion.variable())
    cusp_expansion_rig = P(coeff_list_recognized).O(len(coeff_list_recognized))
    cusp_expansions = dict()
    cusp_expansions[Cusp(1,0)] = cusp_expansion_rig
    return FourierExpansion(G,weight,cusp_expansions,modform_type,Kv,only_principal_cusp_expansion=True,Kw=Kw,Kv=Kv,u_interior_Kv=u_interior_Kv)

def to_reduced_row_echelon_form(fourier_expansions):
    """
    Given a list of instances of "FourierExpansion" return a list of "FourierExpansion" normalized to reduced row echelon form (at infinity).
    """
    rowCount = len(fourier_expansions)
    expansion_degree = min([fourier_expansion.cusp_expansions[Cusp(1,0)].prec()-1 for fourier_expansion in fourier_expansions])
    if expansion_degree < len(fourier_expansions):
        raise ArithmeticError("Not enough terms to transform basis into reduced row echelon form!")
    columnCount = expansion_degree
    fourier_expansions_copy = copy(fourier_expansions)
    #We use the algorithm from https://en.wikipedia.org/wiki/Row_echelon_form
    lead = 0
    for r in range(rowCount):
        if columnCount <= lead:
            break
        i = r
        while fourier_expansions_copy[i].cusp_expansions[Cusp(1,0)][lead] == 0:
            i += 1
            if rowCount == i:
                i = r
                lead += 1
                if columnCount == lead:
                    break
        if i != r:
            fourier_expansions_copy[i], fourier_expansions_copy[r] = fourier_expansions_copy[r], fourier_expansions_copy[i]
        tmp = 1/(fourier_expansions_copy[r].cusp_expansions[Cusp(1,0)][lead])
        fourier_expansions_copy[r] *= tmp
        for i in range(rowCount):
            if i != r:
                tmp = fourier_expansions_copy[i].cusp_expansions[Cusp(1,0)][lead]
                fourier_expansions_copy[i] -= fourier_expansions_copy[r]*tmp
        lead += 1
    return fourier_expansions_copy

class FourierExpansion():
    """
    Class for storing Fourier expansions (q-expansions) of modular forms over general numberfields.
    We use the notation q_N = exp(2*pi*I*z/N).
    """
    def __init__(self, G, weight, cusp_expansions, modform_type, base_ring, only_principal_cusp_expansion=False, Kw=None, Kv=None, u_interior_Kv=None):
        """
        Parameters
        ----------
        G : Group
        weight : weight
        cusp_expansions : dictionary with cusp representatives as keys and Laurent series as values
        modform_type: One of "CuspForm", "ModForm", "Hauptmodul"
        base_ring: Specifies the base_ring over which the q-expansions are defined. Note that for CBFs the precision might be higher for different cusps.
        only_principal_cusp_expansion: Boolean that lets one decide if only the expansion at infinity should be considered
        """
        self.G = G
        self.weight = weight
        self.cusp_expansions = cusp_expansions
        if modform_type not in ("CuspForm","ModForm","Hauptmodul"):
            raise ArithmeticError("Invalid modform_type!")
        self.modform_type = modform_type
        self.base_ring = base_ring
        self.only_principal_cusp_expansion = only_principal_cusp_expansion
        self._Kw = Kw
        self._Kv = Kv
        self._u_interior_Kv = u_interior_Kv

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        trunc_order = min(10,self.cusp_expansions[Cusp(1,0)].prec())
        base_ring = self.cusp_expansions[Cusp(1,0)].base_ring()
        if isinstance(base_ring,NumberField_generic) == True: #q-expansion is defined over numberfield so we can factor it
            diplay_u = self._u_interior_Kv != 1
            c_str = self.get_cusp_expansion(Cusp(1,0),trunc_order=trunc_order,factor_into_u_v=diplay_u).__str__()
        else: #We only know the numerical values of the q-expansions
            c_str = self.get_cusp_expansion(Cusp(1,0),trunc_order=trunc_order,factor_into_u_v=False).__str__() #We only know the numerical values of the q-expansions
        return self.modform_type + " of weight " + str(self.weight) + " with leading order expansion at infinity given by:\n" + c_str

    def _get_convergence_truncated_instance(self):
        """
        Return a new instance for which the expansions at each cusp are not larger than what is needed
        to obtain convergence up to 10^(-digit_prec).
        This function can be used to truncate results that have been obtained from Hejhal's method before computing Petersson products.
        """
        new_self = deepcopy(self)
        digit_prec = bits_to_digits(new_self.base_ring.precision())
        trunc_orders_convergence = get_trunc_orders_convergence(new_self.G,new_self.weight,digit_prec)
        for c in new_self.G.cusps():
            if new_self.cusp_expansions[c].prec() > trunc_orders_convergence[c]+8:
                new_self.cusp_expansions[c] = new_self.cusp_expansions[c].O(trunc_orders_convergence[c]+8)
        return new_self

    def get_cusp_expansion(self, c, trunc_order=None, digit_prec=None, factor_into_u_v=False):
        """
        Returns Fourier expansion at cusp representative 'c' truncated to 'trunc_order' terms and with 'digit_prec' digits precision.
        """
        cusp_expansion = self.cusp_expansions[c]
        if trunc_order == None:
            trunc_order = cusp_expansion.prec()
        elif trunc_order > cusp_expansion.prec():
            raise ArithmeticError("Truncation order is beyond the amount of computed coefficients!")
        
        if digit_prec == None:
            if factor_into_u_v == True:
                base_ring = cusp_expansion.base_ring()
                if isinstance(base_ring,NumberField_generic) == False:
                    raise ArithmeticError("We can only factor expressions into u and v if they are defined over a numberfield!")
                if self._u_interior_Kv == None:
                    raise ArithmeticError("Invalid construction")
                cusp_width = self.G.cusp_width(c)
                return factor_q_expansion_into_u_v(cusp_expansion,self._u_interior_Kv,cusp_width,trunc_order)
            return cusp_expansion.O(trunc_order)
        else:
            if isinstance(cusp_expansion.base_ring(),ComplexField): #Results are non-rigorous
                CF = ComplexField(digits_to_bits(digit_prec))
                return cusp_expansion.O(trunc_order).change_ring(CF)
            #Return rigorous error bounds
            CBF = ComplexBallField(digits_to_bits(digit_prec))
            return cusp_expansion.O(trunc_order).change_ring(CBF)
    
    def _convert_to_CC(self, bit_prec=None):
        """
        Return new instance of FourierExpansion with coefficients converted to a ComplexField.
        Note that this function also sets the precision of the expansions at all cusps to be the same.
        """
        base_ring = self.base_ring
        if bit_prec == None:
            bit_prec = base_ring.precision()
        CC = ComplexField(bit_prec)
        cusp_expansions_new = dict()
        for c in self.cusp_expansions.keys():
            trunc_order = self.cusp_expansions[c].prec()
            new_trunc_order = trunc_order #new_trunc_order denotes the amount of terms that have non-empty error balls
            if isinstance(base_ring,ComplexBallField) and get_decimal_digit_prec(self.cusp_expansions[c][trunc_order-1].rad()) < 0:
                if abs(get_decimal_digit_prec(abs(self.cusp_expansions[c][trunc_order-1]))) < abs(get_decimal_digit_prec(self.cusp_expansions[c][trunc_order-1].rad())):
                    raise ArithmeticError("The q-expansion contains empty error balls which cannot be rounded to ComplexFields. This might effect future results so we decide to stop here.")
            cusp_expansions_new[c] = self.cusp_expansions[c].change_ring(CC).O(new_trunc_order)
        return FourierExpansion(self.G,self.weight,cusp_expansions_new,self.modform_type,CC,only_principal_cusp_expansion=self.only_principal_cusp_expansion)

    def _set_constant_coefficients_to_zero_inplace(self):
        """
        When working with numerical expressions of cuspforms it might happen that the constant terms at the cusps outside infinity are
        only effectively but not truely zero. This function sets these coefficients to zero as well (INPLACE!).
        """
        if self.modform_type != "CuspForm":
            raise ArithmeticError("This function should only be used for cuspforms!")
        for c in self.cusp_expansions.keys():
            if c != Cusp(1,0) and self.cusp_expansions[c][0].is_zero() == False: #For some weird reason, self.cusp_expansions[c][0] != 0 does not work
                P = self.cusp_expansions[c].parent()
                cusp_expansions_list = list(self.cusp_expansions[c]) #We need to work with lists because power series are immutable
                cusp_expansions_list[0] = 0 #Set this coefficient to be truely zero
                self.cusp_expansions[c] = P(cusp_expansions_list).O(self.cusp_expansions[c].prec()) #update expression in class instance

    def __add__(self, a):
        """
        Return self+a where "a" is a constant factor or another instance of "FourierExpansion".
        """
        cusp_expansions = dict()
        cusp_expansions_self = self.cusp_expansions
        weight = self.weight
        if isinstance(a,FourierExpansion):
            for c in cusp_expansions_self.keys():
                cusp_expansions[c] = cusp_expansions_self[c]+a.get_cusp_expansion(c)
            if a.weight != weight:
                raise ArithmeticError("Please only add forms of equal weight!")
        else: #Scalar addition
            for c in cusp_expansions_self.keys():
                cusp_expansion_base_ring = cusp_expansions_self[c].base_ring()
                a_converted = cusp_expansion_base_ring(a) #This is useful if a is a numberfield element and cusp_expansions are defined over CBFs of different prec
                cusp_expansions[c] = cusp_expansions_self[c]+a_converted
        return FourierExpansion(self.G,weight,cusp_expansions,self.modform_type,self.base_ring,
                only_principal_cusp_expansion=self.only_principal_cusp_expansion,Kw=self._Kw,Kv=self._Kv,u_interior_Kv=self._u_interior_Kv)
    
    def __sub__(self, a):
        """
        Return self-a where "a" is a constant factor or another instance of "FourierExpansion".
        """
        cusp_expansions = dict()
        cusp_expansions_self = self.cusp_expansions
        weight = self.weight
        if isinstance(a,FourierExpansion):
            for c in cusp_expansions_self.keys():
                cusp_expansions[c] = cusp_expansions_self[c]-a.get_cusp_expansion(c)
            if a.weight != weight:
                raise ArithmeticError("Please only subtract forms of equal weight!")
        else: #Scalar subtraction
            for c in cusp_expansions_self.keys():
                cusp_expansion_base_ring = cusp_expansions_self[c].base_ring()
                a_converted = cusp_expansion_base_ring(a) #This is useful if a is a numberfield element and cusp_expansions are defined over CBFs of different prec
                cusp_expansions[c] = cusp_expansions_self[c]-a_converted
        return FourierExpansion(self.G,weight,cusp_expansions,self.modform_type,self.base_ring,
                only_principal_cusp_expansion=self.only_principal_cusp_expansion,Kw=self._Kw,Kv=self._Kv,u_interior_Kv=self._u_interior_Kv)
    
    def __mul__(self, a):
        """
        Return self*a where "a" is a constant factor or another instance of "FourierExpansion".
        """
        cusp_expansions = dict()
        cusp_expansions_self = self.cusp_expansions
        weight = self.weight
        modform_type = self.modform_type
        if isinstance(a,FourierExpansion):
            for c in cusp_expansions_self.keys():
                cusp_expansions[c] = cusp_expansions_self[c]*a.get_cusp_expansion(c)
            weight += a.weight
            if a.modform_type == "CuspForm":
                self.modform_type = "CuspForm" #Even if self is a modform, multiplying a modform and a cuspform results in a modform
        else: #Scalar multiplication
            for c in cusp_expansions_self.keys():
                cusp_expansion_base_ring = cusp_expansions_self[c].base_ring()
                a_converted = cusp_expansion_base_ring(a) #This is useful if a is a numberfield element and cusp_expansions are defined over CBFs of different prec
                cusp_expansions[c] = cusp_expansions_self[c]*a_converted
        return FourierExpansion(self.G,weight,cusp_expansions,modform_type,self.base_ring,
                only_principal_cusp_expansion=self.only_principal_cusp_expansion,Kw=self._Kw,Kv=self._Kv,u_interior_Kv=self._u_interior_Kv)
    
    def __truediv__(self, a):
        """
        Return self/a where "a" is a constant factor or another instance of "FourierExpansion".
        """
        cusp_expansions = dict()
        cusp_expansions_self = self.cusp_expansions
        weight = self.weight
        if isinstance(a,FourierExpansion):
            for c in cusp_expansions_self.keys():
                cusp_expansions[c] = cusp_expansions_self[c]/a.get_cusp_expansion(c)
            weight -= a.weight
        else: #Scalar division
            for c in cusp_expansions_self.keys():
                cusp_expansion_base_ring = cusp_expansions_self[c].base_ring()
                a_converted = cusp_expansion_base_ring(a) #This is useful if a is a numberfield element and cusp_expansions are defined over CBFs of different prec
                cusp_expansions[c] = cusp_expansions_self[c]/a_converted
        return FourierExpansion(self.G,weight,cusp_expansions,self.modform_type,self.base_ring,
                only_principal_cusp_expansion=self.only_principal_cusp_expansion,Kw=self._Kw,Kv=self._Kv,u_interior_Kv=self._u_interior_Kv)
    
    def __pow__(self, n):
        """
        Return self**n where "n" is an integer
        """
        cusp_expansions = dict()
        cusp_expansions_self = self.cusp_expansions
        weight = n*self.weight
        for c in cusp_expansions_self.keys():
            cusp_expansions[c] = cusp_expansions_self[c]**n
        return FourierExpansion(self.G,weight,cusp_expansions,self.modform_type,self.base_ring,
                only_principal_cusp_expansion=self.only_principal_cusp_expansion,Kw=self._Kw,Kv=self._Kv,u_interior_Kv=self._u_interior_Kv)
    
    def __one__(self):
        """
        Return one-element of self.
        """
        cusp_expansions = dict()
        cusp_expansions_self = self.cusp_expansions
        for c in cusp_expansions_self.keys():
            cusp_expansions[c] = cusp_expansions_self[c].parent().one()
        return FourierExpansion(self.G,0,cusp_expansions,self.modform_type,self.base_ring,
                only_principal_cusp_expansion=self.only_principal_cusp_expansion,Kw=self._Kw,Kv=self._Kv,u_interior_Kv=self._u_interior_Kv)