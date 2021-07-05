from copy import copy, deepcopy

from sage.rings.real_mpfr import RealField
from sage.rings.complex_field import ComplexField
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.power_series_poly import PowerSeries_poly
from sage.modular.cusps import Cusp

from psage.modform.arithgroup.mysubgroups_alg import SL2Z_elt

from classes.acb_mat_class import Acb_Mat
from point_matching.point_matching_arb_wrap import (
    get_coefficients_cuspform_ir_arb_wrap, get_coefficients_modform_ir_arb_wrap, digits_to_bits, _get_normalization_cuspforms, 
    _get_normalization_modforms, get_pi_ball, get_cuspform_basis_ir_arb_wrap, get_modform_basis_ir_arb_wrap
)
from pullback.my_pullback import my_pullback_general_group_dp, simple_two_by_two_matmul, apply_moebius_transformation_arb_wrap

class ApproxModForm():
    """
    This class contains an approximation of a modular form with coefficients given by approximately 'digit_prec' digits precision.
    """
    def __init__(self,S,digit_prec,modform_type="CuspForm",Y=0,M_0=0,label=0,c_vec=None,prec_loss=None,use_FFT=False,use_Horner=False):
        if modform_type == "CuspForm":
            starting_order = 1
            normalization = _get_normalization_cuspforms(S,label=label)
        elif modform_type == "ModForm":
            starting_order = 0
            normalization = _get_normalization_modforms(S,label=label)
        else:
            raise ArithmeticError("Specified modform_type is not supported yet. Please choose between 'CuspForm' & 'ModForm'")
        G = S.group()
        bit_prec = digits_to_bits(digit_prec)
        RF = RealField(bit_prec)
        CF = ComplexField(bit_prec)
        cusp_expansions = dict()
        if c_vec == None: #We compute c_vec from scratch
            if modform_type == "CuspForm":
                c_vec, M_0 = get_coefficients_cuspform_ir_arb_wrap(S,digit_prec,Y=Y,M_0=M_0,return_M=True,use_FFT=use_FFT,use_Horner=use_Horner,label=label,prec_loss=prec_loss)
            elif modform_type == "ModForm":
                c_vec, M_0 = get_coefficients_modform_ir_arb_wrap(S,digit_prec,Y=Y,M_0=M_0,return_M=True,use_FFT=use_FFT,use_Horner=use_Horner,label=label,prec_loss=prec_loss)
        else: #We construct ApproxModForm from previously computed solution
            if M_0 == 0:
                raise ArithmeticError("Cannot construct ApproxModForm from c_vec without specifying M_0!")

        c_vec_mcbd = c_vec._get_mcbd(bit_prec)

        for ci in G._cusps:
            cii = G._cusps.index(ci)
            cusp_width = G.cusp_width(ci)
            var_string = 'q_' + str(cusp_width)
            R = PowerSeriesRing(CF, var_string)
            q = R.gen()
            f = PowerSeries_poly(R,prec=M_0) #This corresponds to O(q**M_0) in sage syntax
            normalization_len = len(normalization[cii])
            for i in range(normalization_len): #Set normalized coefficients if there are any
                if normalization[cii][i] != 0:
                    f += normalization[cii][i]*q**(starting_order+i)
            for i in range(M_0):
                f += CF(c_vec_mcbd[i+cii*M_0][0])*q**(i+starting_order+normalization_len)
            cusp_expansions[Cusp(ci)] = f


        #Now to some group variables that we store
        self.M_0 = M_0
        self.G = G
        self.weight = S.weight()
        self.digit_prec = digit_prec
        self.modform_type = modform_type
        self.cusp_expansions = cusp_expansions
        self._pi = RF(get_pi_ball(bit_prec))
        self.CF = CF
        self.bit_prec = bit_prec
    
    def get_cusp_expansion(self, c, trunc_order=None, digit_prec=None):
        """
        Returns expansion of modform at cusp representative 'c' truncated to 'trunc_order' terms and with 'digit_prec' digits precision.
        """
        if trunc_order != None:
            if trunc_order > self.M_0:
                raise ArithmeticError("Truncation order is beyond the amount of computed coefficients!")
        else:
            trunc_order = self.M_0
        if digit_prec != None:
            bit_prec = digits_to_bits(digit_prec)
            CF = ComplexField(bit_prec)
            cusp_expansion = self.cusp_expansions[c]
            var_string = cusp_expansion.variable()
            R = PowerSeriesRing(CF, var_string)
            q = R.gen()
            f = PowerSeries_poly(R,prec=trunc_order)
            for i in range(trunc_order):
                f += CF(cusp_expansion[i])*q**i
            return f
        else:
            return self.cusp_expansions[c].add_bigoh(trunc_order)

    def get_cusp_expansions(self, trunc_order=None, digit_prec=None):
        cusp_expansions = dict()
        for c in self.G.cusps():
            cusp_expansions[c] = self.get_cusp_expansion(c,trunc_order=trunc_order,digit_prec=digit_prec)
        return cusp_expansions
    
    def _get_cusp_expansions_dict(self, trunc_order=None, digit_prec=None):
        """
        Return cusp-expansions truncated to 'trunc_order' and 'digit_prec' digits precision as a dictionary.
        The output contains additional parameters such as the weight and is designed to be usable
        through native sage, i.e., without custom code available.
        """
        cusp_expansions = dict()
        for c in self.G.cusps():
            cusp_expansions[c] = self.get_cusp_expansion(c,trunc_order=trunc_order,digit_prec=digit_prec)
        tmp = cusp_expansions.values()[0] #This is the power series at the first cusp
        trunc_order = tmp.prec()
        bit_prec = tmp[0].prec()
        weight = self.weight

        #Now prepare object that we return
        cusp_expansions_dict = dict()
        cusp_expansions_dict['cusp_expansions'] = cusp_expansions
        cusp_expansions_dict['trunc_order'] = trunc_order
        cusp_expansions_dict['bit_prec'] = bit_prec
        cusp_expansions_dict['weight'] = weight
        cusp_expansions_dict['index'] = self.G.index()
        return cusp_expansions_dict
    
    def evaluate(self, z):
        """
        Evaluate Modform at complex point 'z'. In order to achieve ideal convergence, we pullback the point into the fundamental domain
        and afterwards choose the best cusp expansion for evaluation.
        (This function is not working correctly yet)
        """
        print("This function is not working correctly yet!")
        G = self.G
        x,y = float(z.real()),float(z.imag()) #We don't need to perform this computation with arbs...
        x1,y1,T_a,T_b,T_c,T_d = my_pullback_general_group_dp(G,x,y,ret_mat=1)
        vjj = G.closest_vertex(x1,y1,as_integers=1)
        cjj = G._vertex_data[vjj]['cusp'] #closest cusp
        cj = Cusp(G._cusps[cjj])
        U_w = G._vertex_data[vjj]['cusp_map']
        tmp = G.cusp_normalizer(cj).inverse()*U_w*SL2Z_elt(T_a,T_b,T_c,T_d) #This computation involves only integers

        a,b,c,d = tmp[0],tmp[1],tmp[2],tmp[3]
        z_fund = apply_moebius_transformation_arb_wrap(z,a,b,c,d)
        f = self.cusp_expansions[cj].polynomial()
        width = G.cusp_width(cj)
        q_fund = (self.CF(0,2*self._pi)*z_fund/width).exp()

        f_tmp = self.cusp_expansions[Cusp(1,0)].polynomial()
        q_tmp = (self.CF(0,2*self._pi)*z/G.cusp_width(Cusp(1,0))).exp()
        print(cj)
        print(f_tmp(q_tmp))

        automorphy_fact = (c*z+d)**(-self.weight)
        return automorphy_fact*f(q_fund)
    
    def _scal_mul(self, a):
        """
        Return a*Modform where 'a' is a constant factor.
        """
        res = deepcopy(self)
        cusp_expansions = res.cusp_expansions
        for c in cusp_expansions.keys():
            cusp_expansions[c] *= a
        return res
    
    def __add__(self, g):
        """
        Return self+g.
        """
        res = deepcopy(self)
        cusp_expansions = res.cusp_expansions
        cusp_expansions_g = g.cusp_expansions
        for c in cusp_expansions.keys():
            cusp_expansions[c] += cusp_expansions_g[c]
        return res
    
    def __sub__(self, g):
        """
        Return self-g.
        """
        res = deepcopy(self)
        cusp_expansions = res.cusp_expansions
        cusp_expansions_g = g.cusp_expansions
        for c in cusp_expansions.keys():
            cusp_expansions[c] -= cusp_expansions_g[c]
        return res

class CuspExpansions():
    """
    This class is designed to store the cusp expansions of a modular form. It allows convenient storing/loading and does not require any
    custom routines to work.
    """
    def __init__(self,cusp_expansions_dict,trunc_order,bit_prec,weight):
        self.cusp_expansions = cusp_expansions_dict
        self.trunc_order = trunc_order
        self.bit_prec = bit_prec
        self.weight = weight
    
    def get_cusp_expansion(self, c):
        """
        Returns expansion of modform at cusp representative 'c'.
        """
        return self.cusp_expansions[c]

def get_approxmodform_basis(S,digit_prec,modform_type="CuspForm",Y=0,M_0=0,labels=None,prec_loss=None):
    if modform_type == "CuspForm":
        c_vecs, M_0, labels = get_cuspform_basis_ir_arb_wrap(S,digit_prec,Y=Y,M_0=M_0,return_M_and_labels=True,labels=labels,prec_loss=prec_loss)
    elif modform_type == "ModForm":
        c_vecs, M_0, labels = get_modform_basis_ir_arb_wrap(S,digit_prec,Y=Y,M_0=M_0,return_M_and_labels=True,labels=labels,prec_loss=prec_loss)
    else:
        raise ArithmeticError("Specified modform_type is not supported yet. Please choose between 'CuspForm' & 'ModForm'")
    basis = []
    for i in range(len(c_vecs)):
        basis.append(ApproxModForm(S,digit_prec,modform_type=modform_type,Y=Y,M_0=M_0,label=labels[i],c_vec=c_vecs[i]))
    
    return basis
