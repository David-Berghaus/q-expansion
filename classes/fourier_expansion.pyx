from copy import copy, deepcopy

from sage.modular.cusps import Cusp
from sage.rings.complex_field import ComplexField
from sage.rings.complex_arb import ComplexBallField
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.laurent_series_ring import LaurentSeriesRing

from point_matching.point_matching_arb_wrap import (
    get_coefficients_cuspform_ir_arb_wrap, get_coefficients_modform_ir_arb_wrap, get_coefficients_haupt_ir_arb_wrap,
    digits_to_bits, _get_normalization_cuspforms, 
    _get_normalization_modforms, get_cuspform_basis_ir_arb_wrap, get_modform_basis_ir_arb_wrap
)

def get_cuspform_q_expansion_approx(S, digit_prec, Y=0, M_0=0, label=0, c_vec=None, prec_loss=None, use_FFT=True, use_Horner=False):
    """
    Computes q-expansion of cuspform numerically and returns result as instance of "FourierExpansion".
    """
    starting_order = 1
    normalization = _get_normalization_cuspforms(S,label=label)
    if c_vec == None: #We compute c_vec from scratch
        c_vec, M_0 = get_coefficients_cuspform_ir_arb_wrap(S,digit_prec,Y=Y,M_0=M_0,return_M=True,use_FFT=use_FFT,use_Horner=use_Horner,label=label,prec_loss=prec_loss)
    else: #We construct ApproxModForm from previously computed solution
        if M_0 == 0:
            raise ArithmeticError("Cannot construct FourierExpansion from c_vec without specifying M_0!")
    bit_prec = digits_to_bits(digit_prec)
    c_vec_mcbd = c_vec._get_mcbd(bit_prec)
    cusp_expansions = c_vec_to_cusp_expansions(c_vec_mcbd,S,starting_order,normalization,M_0,True)
    return FourierExpansion(S.group(),S.weight(),cusp_expansions,"CuspForm")

def get_modform_q_expansion_approx(S, digit_prec, Y=0, M_0=0, label=0, c_vec=None, prec_loss=None, use_FFT=True, use_Horner=False):
    """
    Computes q-expansion of modform numerically and returns result as instance of "FourierExpansion".
    """
    starting_order = 0
    normalization = _get_normalization_modforms(S,label=label)
    if c_vec == None: #We compute c_vec from scratch
        c_vec, M_0 = get_coefficients_modform_ir_arb_wrap(S,digit_prec,Y=Y,M_0=M_0,return_M=True,use_FFT=use_FFT,use_Horner=use_Horner,label=label,prec_loss=prec_loss)
    else: #We construct ApproxModForm from previously computed solution
        if M_0 == 0:
            raise ArithmeticError("Cannot construct FourierExpansion from c_vec without specifying M_0!")
    bit_prec = digits_to_bits(digit_prec)
    c_vec_mcbd = c_vec._get_mcbd(bit_prec)
    cusp_expansions = c_vec_to_cusp_expansions(c_vec_mcbd,S,starting_order,normalization,M_0,True)
    return FourierExpansion(S.group(),S.weight(),cusp_expansions,"ModForm")

def get_hauptmodul_q_expansion_approx(S, digit_prec, Y=0, M_0=0, c_vec=None, prec_loss=None, use_FFT=True, use_Horner=False):
    """
    Computes q-expansion of hauptmodul numerically and returns result as instance of "FourierExpansion".
    """
    if c_vec == None: #We compute c_vec from scratch
        c_vec, M_0 = get_coefficients_haupt_ir_arb_wrap(S,digit_prec,only_principal_expansion=False,Y=Y,M_0=M_0,return_M=True,use_FFT=use_FFT,use_Horner=use_Horner,prec_loss=prec_loss)
    else: #We construct ApproxModForm from previously computed solution
        if M_0 == 0:
            raise ArithmeticError("Cannot construct FourierExpansion from c_vec without specifying M_0!")
    bit_prec = digits_to_bits(digit_prec)
    c_vec_mcbd = c_vec._get_mcbd(bit_prec)
    cusp_expansions = c_vec_to_cusp_expansions_hauptmodul(c_vec_mcbd,S,M_0)
    return FourierExpansion(S.group(),S.weight(),cusp_expansions,"Hauptmodul")

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
        cusp_expansions = c_vec_to_cusp_expansions(c_vec_mcbd,S,starting_order,normalization,M_0,True)
        basis.append(FourierExpansion(S.group(),S.weight(),cusp_expansions,"CuspForm"))
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
        cusp_expansions = c_vec_to_cusp_expansions(c_vec_mcbd,S,starting_order,normalization,M_0,True)
        basis.append(FourierExpansion(S.group(),S.weight(),cusp_expansions,"ModForm"))
    return basis

def c_vec_to_cusp_expansions(c_vec_mcbd, S, starting_order, normalization, M_0, rescale_hejhal_coefficients):
    """
    Given a c_vec of coefficients as a acb-matrix in Sage, convert expression to a dictionary with cusps as keys and Laurent series as values.
    If rescale_hejhal_coefficients == True, we multiply expansions at non-inf cusps by cusp_width**(-weight_half) because
    our implementation of Hejhal's method uses special cusp-normalizers that absorb the cusp width.
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
        if rescale_hejhal_coefficients == True and Cusp(ci) != Cusp(1,0):
            f /= cusp_width**(weight//2)
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

def to_reduced_row_echelon_form(fourier_expansions):
    """
    Given a list of instances of "FourierExpansion" return a list of "FourierExpansion" normalized to reduced row echelon form (at infinity).
    """
    rowCount = len(fourier_expansions)
    expansion_degree = min([fourier_expansion.cusp_expansions[Cusp(1,0)].prec()-1 for fourier_expansion in fourier_expansions])
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
    def __init__(self, G, weight, cusp_expansions, modform_type):
        """
        Parameters
        ----------
        G : Group
        weight : weight
        cusp_expansions : dictionary with cusp representatives as keys and Laurent series as values
        modform_type: One of "CuspForm", "ModForm", "Hauptmodul"
        """
        self.G = G
        self.weight = weight
        self.cusp_expansions = cusp_expansions
        if modform_type not in ("CuspForm","ModForm","Hauptmodul"):
            raise ArithmeticError("Invalid modform_type!")
        self.modform_type = modform_type

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        trunc_order = min(10,self.cusp_expansions[Cusp(1,0)].prec())
        c_str = self.get_cusp_expansion(Cusp(1,0),trunc_order=trunc_order).__str__()
        return self.modform_type + " of weight " + str(self.weight) + " with leading order expansion at infinity given by:\n" + c_str 

    def get_cusp_expansion(self, c, trunc_order=None, digit_prec=None):
        """
        Returns Fourier expansion at cusp representative 'c' truncated to 'trunc_order' terms and with 'digit_prec' digits precision.
        """
        cusp_expansion = self.cusp_expansions[c]
        if trunc_order == None:
            trunc_order = cusp_expansion.prec()
        elif trunc_order > cusp_expansion.prec():
            raise ArithmeticError("Truncation order is beyond the amount of computed coefficients!")
        
        if digit_prec == None:
            return cusp_expansion.O(trunc_order)
        else:
            if isinstance(cusp_expansion.base_ring(),ComplexField): #Results are non-rigorous
                CF = ComplexField(digits_to_bits(digit_prec))
                return cusp_expansion.O(trunc_order).change_ring(CF)
            #Return rigorous error bounds
            CBF = ComplexBallField(digits_to_bits(digit_prec))
            return cusp_expansion.O(trunc_order).change_ring(CBF)
    
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
            weight += a.weight
        else: #Scalar addition
            for c in cusp_expansions_self.keys():
                cusp_expansions[c] = cusp_expansions_self[c]+a
        return FourierExpansion(self.G,weight,cusp_expansions,self.modform_type)
    
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
            weight += a.weight
        else: #Scalar subtraction
            for c in cusp_expansions_self.keys():
                cusp_expansions[c] = cusp_expansions_self[c]-a
        return FourierExpansion(self.G,weight,cusp_expansions,self.modform_type)
    
    def __mul__(self, a):
        """
        Return self*a where "a" is a constant factor or another instance of "FourierExpansion".
        """
        cusp_expansions = dict()
        cusp_expansions_self = self.cusp_expansions
        weight = self.weight
        if isinstance(a,FourierExpansion):
            for c in cusp_expansions_self.keys():
                cusp_expansions[c] = cusp_expansions_self[c]*a.get_cusp_expansion(c)
            weight += a.weight
        else: #Scalar multiplication
            for c in cusp_expansions_self.keys():
                cusp_expansions[c] = cusp_expansions_self[c]*a
        return FourierExpansion(self.G,weight,cusp_expansions,self.modform_type)
    
    def __div__(self, a):
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
                cusp_expansions[c] = cusp_expansions_self[c]/a
        return FourierExpansion(self.G,weight,cusp_expansions,self.modform_type)
    
    def __pow__(self, n):
        """
        Return self**n where "n" is an integer
        """
        cusp_expansions = dict()
        cusp_expansions_self = self.cusp_expansions
        weight = n*self.weight
        for c in cusp_expansions_self.keys():
            cusp_expansions[c] = cusp_expansions_self[c]**n
        return FourierExpansion(self.G,weight,cusp_expansions,self.modform_type)
    
    def __one__(self):
        """
        Return one-element of self.
        """
        cusp_expansions = dict()
        cusp_expansions_self = self.cusp_expansions
        for c in cusp_expansions_self.keys():
            cusp_expansions[c] = cusp_expansions_self[c].parent().one()
        return FourierExpansion(self.G,0,cusp_expansions,"ModForm")