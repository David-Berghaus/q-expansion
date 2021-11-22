from copy import copy, deepcopy

from sage.modular.cusps import Cusp
from sage.rings.complex_field import ComplexField

from point_matching.point_matching_arb_wrap import digits_to_bits

# def to_reduced_row_echelon_form(fourier_expansions):
#     """
#     Given a list of instances of "FourierExpansion" return a list of "FourierExpansion" normalized to reduced row echelon form.
#     """
#     M = [fourier_expansion.get_cusp_expansion(Cusp(1,0)).list() for fourier_expansion in fourier_expansions]
#     rowCount = len(M)
#     columnCount = len(M[0])
#     for i in range(rowCount):
#         if len(M[i]) != columnCount:
#             raise ArithmeticError("We need equal amounts of terms for all Fourier expansions!")
#     fourier_expansions_copy = copy(fourier_expansions)
#     #We use the algorithm from https://en.wikipedia.org/wiki/Row_echelon_form
#     lead = 0
#     for r in range(rowCount):
#         if columnCount <= lead:
#             break
#         i = r
#         while M[i][lead] == 0:
#             i += 1
#             if rowCount == i:
#                 i = r
#                 lead += 1
#                 if columnCount == lead:
#                     break
#         if i != r:
#             M[i], M[r] = M[r], M[i]
#             fourier_expansions_copy[i], fourier_expansions_copy[r] = fourier_expansions_copy[r], fourier_expansions_copy[i]
#         tmp = 1/(M[r][lead])
#         M[r] *= tmp
#         fourier_expansions_copy[r] *= tmp
#         for i in range(rowCount):
#             if i != r:
#                 tmp = M[i][lead]
#                 M[i] -= tmp*M[r]
#                 fourier_expansions_copy[i] -= tmp*fourier_expansions_copy[r]
#         lead += 1
#     return fourier_expansions_copy

def to_reduced_row_echelon_form(fourier_expansions):
    """
    Given a list of instances of "FourierExpansion" return a list of "FourierExpansion" normalized to reduced row echelon form.
    """
    rowCount = len(fourier_expansions)
    expansion_degree = fourier_expansions[0].get_cusp_expansion(Cusp(1,0)).prec()-1
    for i in range(rowCount):
        if fourier_expansions[i].get_cusp_expansion(Cusp(1,0)).prec()-1 != expansion_degree:
            raise ArithmeticError("We need equal amounts of terms for all Fourier expansions!")
    columnCount = expansion_degree
    fourier_expansions_copy = copy(fourier_expansions)
    #We use the algorithm from https://en.wikipedia.org/wiki/Row_echelon_form
    lead = 0
    for r in range(rowCount):
        if columnCount <= lead:
            break
        i = r
        while fourier_expansions_copy[i].get_cusp_expansion(Cusp(1,0))[lead] == 0:
            i += 1
            if rowCount == i:
                i = r
                lead += 1
                if columnCount == lead:
                    break
        if i != r:
            fourier_expansions_copy[i], fourier_expansions_copy[r] = fourier_expansions_copy[r], fourier_expansions_copy[i]
        tmp = 1/(fourier_expansions_copy[r].get_cusp_expansion(Cusp(1,0))[lead])
        fourier_expansions_copy[r] *= tmp
        for i in range(rowCount):
            if i != r:
                tmp = fourier_expansions_copy[i].get_cusp_expansion(Cusp(1,0))[lead]
                fourier_expansions_copy[i] -= fourier_expansions_copy[r]*tmp
        lead += 1
    return fourier_expansions_copy

class FourierExpansion():
    """
    Class for storing Fourier expansions (q-expansions) of modular forms over general rings.
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
            CBF = ComplexField(digits_to_bits(digit_prec))
            return cusp_expansion.O(trunc_order).change_ring(CBF)
    
    def __mul__(self, a):
        """
        Return a*self where "a" is a constant factor or another instance of "FourierExpansion".
        """
        if isinstance(a,FourierExpansion):
            raise NotImplementedError("We have not implemented this yet!")
        res = deepcopy(self)
        cusp_expansions = res.cusp_expansions
        for c in cusp_expansions.keys():
            cusp_expansions[c] *= a
        return res
    
    def __add__(self, g):
        """
        Return self+g where "g" is another instance of "FourierExpansion".
        """
        res = deepcopy(self)
        cusp_expansions = res.cusp_expansions
        cusp_expansions_g = g.cusp_expansions
        for c in cusp_expansions.keys():
            cusp_expansions[c] += cusp_expansions_g[c]
        return res
    
    def __sub__(self, g):
        """
        Return self-g where "g" is another instance of "FourierExpansion".
        """
        res = deepcopy(self)
        cusp_expansions = res.cusp_expansions
        cusp_expansions_g = g.cusp_expansions
        for c in cusp_expansions.keys():
            cusp_expansions[c] -= cusp_expansions_g[c]
        return res