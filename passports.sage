# ****************************************************************************
#       Copyright (C) 2022 David Berghaus <berghaus@th.physik.uni-bonn.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import deepcopy
from sys import exit
import os
import json
from textwrap import indent

from psage.modform.maass.automorphic_forms import AutomorphicFormSpace

from belyi.number_fields import is_effectively_zero, get_numberfield_of_coeff, to_K, get_v_Kw, get_u_minpoly
from belyi.expression_in_u_and_v import convert_from_Kv_to_Kw
from belyi.elliptic_curve import get_elliptic_curve
from eisenstein.eisenstein_computation import compute_eisenstein_series, echelon_basis_to_eisenstein_basis
from point_matching.point_matching_arb_wrap import _get_echelon_normalization_from_label, digits_to_bits, get_M_0, get_horo_height_arb_wrap
from classes.fourier_expansion import get_hauptmodul_q_expansion_approx, get_cuspform_basis_approx, get_cuspform_q_expansion_approx, get_modform_basis_approx, recognize_cusp_expansion_using_u, to_reduced_row_echelon_form
from classes.belyi_map import BelyiMap, get_u_str
from classes.factored_polynomial import get_updated_Kw_v_Kw

def compute_passport_data_genus_zero(passport, rigorous_trunc_order, eisenstein_digit_prec, max_weight, compare_result_to_numerics=True, compute_embeddings=True, return_floating_expansions=True, numerics_digit_prec=50, tol=1e-10, state_file_path=None):
    """
    Compute relevant data for a specified passport.
    Input:
    ------
    passport: List of (Galois conjugate) subgroups that are part of the passport.
    rigorous_trunc_order: Amount of terms in the rigorous computation of the q-expansion at infinity
    eisenstein_digit_prec: Digit precision of the numerical approximation of the Eisenstein series
    max_weight: The maximum weight of the computed modular forms
    compare_result_to_numerics: Boolean that decides if the results that have been computed through the Belyi map should be compared to the numerical values
    numerics_digit_prec: The precision at which the numerical computation that we use to compare the results is performed
    tol: Maximum difference between q-expansion coefficients compared to the numerical results
    state_file_path: If specified, restart computation before numerical results are being computed.
    Output:
    ------
    A dictionary with the computed information that can be read using Sage only.
    """
    max_extension_field_degree = get_max_extension_field_degree(passport)
    G = passport[0]
    CC = ComplexField(digits_to_bits(eisenstein_digit_prec))

    if state_file_path != None and os.path.exists(state_file_path) == True:
        B, j_G_rig, cuspforms_rig, modforms_rig = dict_to_state_genus_zero(load(state_file_path))
    else:
        B = BelyiMap(G,max_extension_field_degree=max_extension_field_degree)
        if B._Kv.degree() != max_extension_field_degree:
            if has_equal_list_entry(G.cusp_widths(),0) == False: #If two cusps are identical it sometimes happens that they are in the same numberfield which we do not need to investigate further
                raise ArithmeticError("We have not considered the case of decaying Galois orbits yet!")

        #First do the rigorous computation of the q-expansions
        j_G_rig = B.get_hauptmodul_q_expansion(rigorous_trunc_order)
        cuspforms_rig, modforms_rig = dict(), dict()
        for weight in range(2,max_weight+1,2): #We only consider even weights
            if G.dimension_modular_forms(weight) != 0:
                modforms_rig[weight] = B.get_modforms(weight,rigorous_trunc_order,j_G=j_G_rig)
            if G.dimension_cusp_forms(weight) != 0:
                cuspforms_rig[weight] = B.get_cuspforms(weight,rigorous_trunc_order,j_G=j_G_rig)
    if state_file_path != None and os.path.exists(state_file_path) == False:
        curr_state = state_to_dict_genus_zero(B, j_G_rig, cuspforms_rig, modforms_rig)
        save(curr_state,state_file_path)
        exit(1)

    #Now to the numerical stuff
    eisenstein_trunc_orders = B._get_trunc_orders_convergence(max_weight,eisenstein_digit_prec)
    j_G_fl = B.get_hauptmodul_q_expansion_approx(eisenstein_trunc_orders,eisenstein_digit_prec)
    cuspforms_fl, modforms_fl = dict(), dict()
    eis_scaling_constants = dict()
    eis_scaling_constants_canonical = dict()
    for weight in range(2,max_weight+1,2): #We only consider even weights
        if G.dimension_modular_forms(weight) != 0:
            modforms_fl[weight] = B.get_modforms(weight,eisenstein_trunc_orders,digit_prec=eisenstein_digit_prec,j_G=j_G_fl)
            if G.dimension_cusp_forms(weight) != 0:
                cuspforms_fl[weight] = B.get_cuspforms(weight,eisenstein_trunc_orders,digit_prec=eisenstein_digit_prec,j_G=j_G_fl)
                eisforms_fl, eis_scaling_constant_list = compute_eisenstein_series(cuspforms_fl[weight],modforms_fl[weight],return_scaling_constants=True)
                for i in range(len(eis_scaling_constant_list)):
                    for j in range(len(eis_scaling_constant_list[i])):
                        if eis_scaling_constant_list[i][j] != 0 and eis_scaling_constant_list[i][j] != 1 and is_effectively_zero(eis_scaling_constant_list[i][j],int(round(0.99*eisenstein_digit_prec))) == True:
                            eis_scaling_constant_list[i][j] = 0 #We have a numerical zero which we now set to a true zero
            else:
                eisforms_fl = modforms_fl[weight] #In this case the eisforms are equivalent to modforms
                eis_scaling_constant_list = [_get_echelon_normalization_from_label(i,len(eisforms_fl)) for i in range(len(eisforms_fl))]
            eis_scaling_constants[weight] = eis_scaling_constant_list
            _, eis_scaling_constants_canonical[weight] = echelon_basis_to_eisenstein_basis(eisforms_fl,return_scaling_constants=True)

    if compare_result_to_numerics == True:
        compare_results_to_numerics(G,max_weight,modforms_rig,cuspforms_rig,eis_scaling_constants,B.u_QQbar,numerics_digit_prec,tol)#Verify results by comparing them to numerical values

    #Convert numerical expressions of Eisenstein series to CC because we cannot hash CBFs
    for weight in eis_scaling_constants:
        for i in range(len(eis_scaling_constants[weight])):
            for j in range(len(eis_scaling_constants[weight][i])):
                eis_scaling_constants[weight][i][j] = CC(eis_scaling_constants[weight][i][j])
        for i in range(len(eis_scaling_constants_canonical[weight])):
            for j in range(len(eis_scaling_constants_canonical[weight][i])):
                eis_scaling_constants_canonical[weight][i][j] = CC(eis_scaling_constants_canonical[weight][i][j])

    res = dict()
    res["G"] = B.G.as_permutation_group()
    res["is_congruence"] = G.is_congruence()
    res["monodromy_group"] = G.perm_group().structure_description()
    res["Kv"] = B._Kv
    res["v"] = B._Kv.gen()
    res["u"] = B.u_QQbar
    res["u_interior_Kv"] = B._u_interior_Kv
    res["v_Kw"] = B._v_Kw
    res["u_str"] = B.get_u_str()
    res["curve"] = B._return_res_as_dict()
    res["q_expansions"] = dict()
    display_u = B._u_interior_Kv != 1
    for weight in range(0,max_weight+1,2): #We only consider even weights
        res["q_expansions"][weight] = dict()
        if weight == 0:
            res["q_expansions"][weight]["hauptmodul_raw"] = j_G_rig.get_cusp_expansion(Cusp(1,0))
            res["q_expansions"][weight]["hauptmodul_pretty"] = j_G_rig.get_cusp_expansion(Cusp(1,0),factor_into_u_v=display_u)
        else:
            if G.dimension_modular_forms(weight) != 0:
                res["q_expansions"][weight]["modforms_raw"] = [modform.get_cusp_expansion(Cusp(1,0)) for modform in modforms_rig[weight]]
                res["q_expansions"][weight]["modforms_pretty"] = [modform.get_cusp_expansion(Cusp(1,0),factor_into_u_v=display_u) for modform in modforms_rig[weight]]
                res["q_expansions"][weight]["eisenstein_basis_factors"] = eis_scaling_constants[weight]
                res["q_expansions"][weight]["eisenstein_canonical_normalizations"] = eis_scaling_constants_canonical[weight]
                if G.dimension_cusp_forms(weight) != 0:
                    res["q_expansions"][weight]["cuspforms_raw"] = [cuspform.get_cusp_expansion(Cusp(1,0)) for cuspform in cuspforms_rig[weight]]
                    res["q_expansions"][weight]["cuspforms_pretty"] = [cuspform.get_cusp_expansion(Cusp(1,0),factor_into_u_v=display_u) for cuspform in cuspforms_rig[weight]]
    if compute_embeddings == True:
        res["embeddings"] = get_all_embeddings(passport,res["q_expansions"],res["Kv"],B._u_interior_Kv,G.cusp_width(Cusp(1,0)))
    if state_file_path != None:
        os.remove(state_file_path)
    if return_floating_expansions == True:
        floating_expansions = dict()
        for weight in range(0,max_weight+1,2): #We only consider even weights
            floating_expansions[weight] = dict()
            if weight == 0:
                floating_expansions[weight]["hauptmodul_float"] = j_G_fl._convert_to_CC() #Note that we cannot store arbs
            else:
                if weight in modforms_fl:
                    floating_expansions[weight]["modforms_float"] = [modform_fl._convert_to_CC() for modform_fl in modforms_fl[weight]] #Note that we cannot store arbs
                if weight in cuspforms_fl:
                    floating_expansions[weight]["cuspforms_float"] = [cuspform_fl._convert_to_CC() for cuspform_fl in cuspforms_fl[weight]] #Note that we cannot store arbs
        return res, floating_expansions
    return res

def state_to_dict_genus_zero(B, j_G_rig, cuspforms_rig, modforms_rig):
    """
    Save the current state into a dictionary which we can store to restart the computation.
    """
    curr_state = dict()
    curr_state['B'] = B
    curr_state['j_G_rig'] = j_G_rig
    curr_state['cuspforms_rig'] = cuspforms_rig
    curr_state['modforms_rig'] = modforms_rig
    return curr_state

def dict_to_state_genus_zero(curr_state):
    """
    Given a dictionary of the current state, unpack the variables.
    """
    return curr_state['B'], curr_state['j_G_rig'], curr_state['cuspforms_rig'], curr_state['modforms_rig']

def compute_passport_data_higher_genera(passport, max_closed_form_trunc_order, digit_prec, max_weight, construct_higher_weight_from_lower_weight_forms=False, compare_result_to_numerics=True, compute_embeddings=True, return_floating_expansions=True, numerics_digit_prec=40, tol=1e-10, state_file_path=None):
    """
    Compute database entry for given passport.
    If state_file_path != None we store intermediate steps and exit the computation.
    """
    if max_closed_form_trunc_order == None:
        max_closed_form_trunc_order = 10000000 #Recognize as many coeffs as possible
    max_extension_field_degree = get_max_extension_field_degree(passport)
    G = passport[0]
    CC = ComplexField(digits_to_bits(digit_prec))
    principal_cusp_width = G.cusp_width(Cusp(1,0))
    modforms_fl = dict()
    modforms_rig = dict()

    #First compute the lowest-weight cuspform space and recognize u from that
    if state_file_path != None and os.path.exists(state_file_path) == True:
        cuspforms_fl, cuspforms_rig, modforms_fl, modforms_rig, Kv, Kw, v_Kw, u_interior_Kv, u, lowest_non_zero_cuspform_weight = dict_to_state_higher_genera(load(state_file_path))
    else:
        cuspforms_fl, cuspforms_rig, Kv, Kw, v_Kw, u_interior_Kv, u, lowest_non_zero_cuspform_weight = compute_lowest_weight_cuspform_space_to_get_u(G,max_closed_form_trunc_order,digit_prec,max_extension_field_degree,principal_cusp_width)
    if state_file_path != None and os.path.exists(state_file_path) == False:
        curr_state = state_to_dict_higher_genera(cuspforms_fl, cuspforms_rig, modforms_fl, modforms_rig, Kv, Kw, v_Kw, u_interior_Kv, u, lowest_non_zero_cuspform_weight)
        save(curr_state,state_file_path)
        exit(1)

    #Then continue with computing modforms
    #We start with these before computing the remaining cuspforms because some of them can be used to construct cuspforms (while the reverse it not true)
    for weight in range(2,max_weight+1,2): #We only consider even weights
        if weight not in modforms_fl and G.dimension_modular_forms(weight) != 0 and G.dimension_eis(weight) != 0:
            compute_modforms_higher_genera(weight, G, digit_prec, Kv, u, v_Kw, u_interior_Kv, modforms_fl, modforms_rig, cuspforms_fl, cuspforms_rig, max_closed_form_trunc_order, construct_higher_weight_from_lower_weight_forms)
            if state_file_path != None:
                curr_state = state_to_dict_higher_genera(cuspforms_fl, cuspforms_rig, modforms_fl, modforms_rig, Kv, Kw, v_Kw, u_interior_Kv, u, lowest_non_zero_cuspform_weight)
                save(curr_state,state_file_path)
                exit(1)

    #Now compute remaining cuspforms
    for weight in range(lowest_non_zero_cuspform_weight+2,max_weight+1,2): #We only consider even weights
        if weight not in cuspforms_fl and G.dimension_cusp_forms(weight) != 0:
            compute_cuspforms_higher_genera(weight, G, digit_prec, Kv, u, v_Kw, u_interior_Kv, modforms_fl, modforms_rig, cuspforms_fl, cuspforms_rig, max_closed_form_trunc_order, construct_higher_weight_from_lower_weight_forms)
            if state_file_path != None:
                curr_state = state_to_dict_higher_genera(cuspforms_fl, cuspforms_rig, modforms_fl, modforms_rig, Kv, Kw, v_Kw, u_interior_Kv, u, lowest_non_zero_cuspform_weight)
                save(curr_state,state_file_path)
                exit(1)

    #Now to the Eisenstein series
    eis_scaling_constants = dict()
    eis_scaling_constants_canonical = dict()
    for weight in range(2,max_weight+1,2): #We only consider even weights
        if G.dimension_eis(weight) != 0:
            if G.dimension_cusp_forms(weight) != 0:
                eisforms_fl, eis_scaling_constant_list = compute_eisenstein_series(cuspforms_fl[weight],modforms_fl[weight],return_scaling_constants=True)
                for i in range(len(eis_scaling_constant_list)):
                    for j in range(len(eis_scaling_constant_list[i])):
                        if eis_scaling_constant_list[i][j] != 0 and eis_scaling_constant_list[i][j] != 1 and is_effectively_zero(eis_scaling_constant_list[i][j],int(round(0.99*digit_prec))) == True:
                            eis_scaling_constant_list[i][j] = 0 #We have a numerical zero which we now set to a true zero
            else:
                eisforms_fl = modforms_fl[weight] #In this case the eisforms are equivalent to modforms
                eis_scaling_constant_list = [_get_echelon_normalization_from_label(i,len(eisforms_fl)) for i in range(len(eisforms_fl))]
            eis_scaling_constants[weight] = eis_scaling_constant_list
            _, eis_scaling_constants_canonical[weight] = echelon_basis_to_eisenstein_basis(eisforms_fl,return_scaling_constants=True)

    u_QQbar = QQbar(u)
    if compare_result_to_numerics == True:
        compare_results_to_numerics(G,max_weight,modforms_rig,cuspforms_rig,eis_scaling_constants,u_QQbar,numerics_digit_prec,tol,Y_fact=0.9)
    
    #Convert numerical expressions of Eisenstein series to CC because we cannot hash CBFs
    for weight in eis_scaling_constants:
        for i in range(len(eis_scaling_constants[weight])):
            for j in range(len(eis_scaling_constants[weight][i])):
                eis_scaling_constants[weight][i][j] = CC(eis_scaling_constants[weight][i][j])
        for i in range(len(eis_scaling_constants_canonical[weight])):
            for j in range(len(eis_scaling_constants_canonical[weight][i])):
                eis_scaling_constants_canonical[weight][i][j] = CC(eis_scaling_constants_canonical[weight][i][j])
    
    res = dict()
    res["G"] = G.as_permutation_group()
    res["is_congruence"] = G.is_congruence()
    res["monodromy_group"] = G.perm_group().structure_description()
    res["Kv"] = Kv
    res["v"] = Kv.gen()
    res["u"] = u_QQbar
    res["u_str"] = get_u_str(u_interior_Kv,principal_cusp_width)
    res["u_interior_Kv"] = u_interior_Kv
    res["v_Kw"] = v_Kw
    if G.genus() == 1:
        res["curve"] = get_elliptic_curve(cuspforms_fl[2][0],Kv,digit_prec)
    else:
        res["curve"] = None
    res["q_expansions"] = dict()
    display_u = u_interior_Kv != 1
    for weight in range(2,max_weight+1,2): #We only consider even weights
        res["q_expansions"][weight] = dict()
        if G.dimension_modular_forms(weight) != 0:
            if G.dimension_eis(weight) == 0: #The cuspforms form the modform basis and we don't have any eisenstein_basis_factors
                res["q_expansions"][weight]["modforms_raw"] = [cuspform.get_cusp_expansion(Cusp(1,0)) for cuspform in cuspforms_rig[weight]]
                res["q_expansions"][weight]["modforms_pretty"] = [cuspform.get_cusp_expansion(Cusp(1,0),factor_into_u_v=display_u) for cuspform in cuspforms_rig[weight]]
            else:
                res["q_expansions"][weight]["modforms_raw"] = [modform.get_cusp_expansion(Cusp(1,0)) for modform in modforms_rig[weight]]
                res["q_expansions"][weight]["modforms_pretty"] = [modform.get_cusp_expansion(Cusp(1,0),factor_into_u_v=display_u) for modform in modforms_rig[weight]]
                res["q_expansions"][weight]["eisenstein_basis_factors"] = eis_scaling_constants[weight]
                res["q_expansions"][weight]["eisenstein_canonical_normalizations"] = eis_scaling_constants_canonical[weight]
            if G.dimension_cusp_forms(weight) != 0:
                res["q_expansions"][weight]["cuspforms_raw"] = [cuspform.get_cusp_expansion(Cusp(1,0)) for cuspform in cuspforms_rig[weight]]
                res["q_expansions"][weight]["cuspforms_pretty"] = [cuspform.get_cusp_expansion(Cusp(1,0),factor_into_u_v=display_u) for cuspform in cuspforms_rig[weight]]
    if compute_embeddings == True:
        res["embeddings"] = get_all_embeddings(passport,res["q_expansions"],res["Kv"],u_interior_Kv,principal_cusp_width)
    if state_file_path != None:
        os.remove(state_file_path)
    if return_floating_expansions == True:
        floating_expansions = dict()
        for weight in range(2,max_weight+1,2): #We only consider even weights
            floating_expansions[weight] = dict()
            if weight in modforms_fl:
                floating_expansions[weight]["modforms_float"] = [modform_fl._convert_to_CC() for modform_fl in modforms_fl[weight]] #Note that we cannot store arbs
            if weight in cuspforms_fl:
                floating_expansions[weight]["cuspforms_float"] = [cuspform_fl._convert_to_CC() for cuspform_fl in cuspforms_fl[weight]] #Note that we cannot store arbs
        return res, floating_expansions
    return res

def state_to_dict_higher_genera(cuspforms_fl, cuspforms_rig, modforms_fl, modforms_rig, Kv, Kw, v_Kw, u_interior_Kv, u, lowest_non_zero_cuspform_weight):
    """
    Save the current state into a dictionary which we can store to restart the computation.
    """
    curr_state = dict()
    curr_state['cuspforms_fl'] = cuspforms_fl
    curr_state['cuspforms_rig'] = cuspforms_rig
    curr_state['modforms_fl'] = modforms_fl
    curr_state['modforms_rig'] = modforms_rig
    curr_state['Kv'] = Kv
    curr_state['Kw'] = Kw
    curr_state['v_Kw'] = v_Kw
    curr_state['u_interior_Kv'] = u_interior_Kv
    curr_state['u'] = u
    curr_state['lowest_non_zero_cuspform_weight'] = lowest_non_zero_cuspform_weight
    return curr_state

def dict_to_state_higher_genera(curr_state):
    """
    Given a dictionary of the current state, unpack the variables.
    """
    return curr_state['cuspforms_fl'], curr_state['cuspforms_rig'], curr_state['modforms_fl'], curr_state['modforms_rig'], curr_state['Kv'], curr_state['Kw'], curr_state['v_Kw'], curr_state['u_interior_Kv'], curr_state['u'], curr_state['lowest_non_zero_cuspform_weight']

def compute_modforms_higher_genera(weight, G, digit_prec, Kv, u, v_Kw, u_interior_Kv, modforms_fl, modforms_rig, cuspforms_fl, cuspforms_rig, max_closed_form_trunc_order, construct_higher_weight_from_lower_weight_forms):
    """
    Compute the modular forms of specified weight and add them to modforms_fl and modforms_rig dicts.
    """
    dim_M = G.dimension_modular_forms(weight)
    if dim_M != 0 and G.dimension_eis(weight) != 0: #If the eisenstein space is zero-dimensional then we could only get cuspforms which have a different normalization...
        if construct_higher_weight_from_lower_weight_forms == False: #Compute everything from scratch
            modforms_fl[weight] = get_modform_basis_approx(AutomorphicFormSpace(G,weight),digit_prec)
            modforms_rig[weight] = [recognize_cusp_expansion_using_u(modforms_fl[weight][label].get_cusp_expansion(Cusp(1,0)),weight,G,max_closed_form_trunc_order,"ModForm",Kv,u,v_Kw,u_interior_Kv) for label in range(dim_M)]
        else:
            modforms_fl_weight, modforms_rig_weight = dict(), dict() #Temporary variables that we use to store forms for each label
            product_formulas = get_product_formulas_for_forms(G,weight,False)
            constructable_labels, computable_labels = get_constructable_and_computable_labels(product_formulas,dim_M)
            for constructable_label in constructable_labels:
                product_formula = product_formulas[constructable_label]
                modforms_fl_weight[constructable_label] = construct_form_from_product_formula(product_formula,modforms_fl,cuspforms_fl)
                modforms_rig_weight[constructable_label] = construct_form_from_product_formula(product_formula,modforms_rig,cuspforms_rig)
            modforms_fl_computed = get_modform_basis_approx(AutomorphicFormSpace(G,weight),digit_prec,labels=computable_labels)
            for i in range(len(modforms_fl_computed)):
                label = computable_labels[i]
                modforms_fl_weight[label] = modforms_fl_computed[i]
                modforms_rig_weight[label] = recognize_cusp_expansion_using_u(modforms_fl_computed[i].get_cusp_expansion(Cusp(1,0)),weight,G,max_closed_form_trunc_order,"ModForm",Kv,u,v_Kw,u_interior_Kv)
            #Now we got full bases for the spaces which we only need to transform into reduced-row-echelon-form
            modforms_fl[weight] = to_reduced_row_echelon_form([modforms_fl_weight[i] for i in range(dim_M)])
            modforms_rig[weight] = to_reduced_row_echelon_form([modforms_rig_weight[i] for i in range(dim_M)])

def compute_cuspforms_higher_genera(weight, G, digit_prec, Kv, u, v_Kw, u_interior_Kv, modforms_fl, modforms_rig, cuspforms_fl, cuspforms_rig, max_closed_form_trunc_order, construct_higher_weight_from_lower_weight_forms):
    """
    Compute the cusp forms of specified weight and add them to cuspforms_fl and cuspforms_rig dicts.
    """
    dim_S = G.dimension_cusp_forms(weight)
    if dim_S != 0:
        if construct_higher_weight_from_lower_weight_forms == False: #Compute everything from scratch
            cuspforms_fl[weight] = get_cuspform_basis_approx(AutomorphicFormSpace(G,weight),digit_prec)
            cuspforms_rig[weight] = [recognize_cusp_expansion_using_u(cuspforms_fl[weight][label].get_cusp_expansion(Cusp(1,0)),weight,G,max_closed_form_trunc_order,"CuspForm",Kv,u,v_Kw,u_interior_Kv) for label in range(dim_S)]
        else:
            cuspforms_fl_weight, cuspforms_rig_weight = dict(), dict() #Temporary variables that we use to store forms for each label
            product_formulas = get_product_formulas_for_forms(G,weight,True)
            constructable_labels, computable_labels = get_constructable_and_computable_labels(product_formulas,dim_S)
            for constructable_label in constructable_labels:
                product_formula = product_formulas[constructable_label]
                cuspforms_fl_weight[constructable_label] = construct_form_from_product_formula(product_formula,modforms_fl,cuspforms_fl)
                cuspforms_rig_weight[constructable_label] = construct_form_from_product_formula(product_formula,modforms_rig,cuspforms_rig)
            if len(computable_labels) != 0:
                cuspforms_fl_computed = get_cuspform_basis_approx(AutomorphicFormSpace(G,weight),digit_prec,labels=computable_labels)
            else:
                cuspforms_fl_computed = []
            for i in range(len(cuspforms_fl_computed)):
                label = computable_labels[i]
                cuspforms_fl_weight[label] = cuspforms_fl_computed[i]
                cuspforms_rig_weight[label] = recognize_cusp_expansion_using_u(cuspforms_fl_computed[i].get_cusp_expansion(Cusp(1,0)),weight,G,max_closed_form_trunc_order,"CuspForm",Kv,u,v_Kw,u_interior_Kv)
            #Now we got full bases for the spaces which we only need to transform into reduced-row-echelon-form
            cuspforms_fl[weight] = to_reduced_row_echelon_form([cuspforms_fl_weight[i] for i in range(dim_S)])
            cuspforms_rig[weight] = to_reduced_row_echelon_form([cuspforms_rig_weight[i] for i in range(dim_S)])

def compute_lowest_weight_cuspform_space_to_get_u(G, max_closed_form_trunc_order, digit_prec, max_extension_field_degree, principal_cusp_width):
    """
    Compute lowest weight non-empty space of cuspforms to determine u from one of the cuspforms.
    Note that one needs to be careful that the selected form is not an oldform.
    For our examples, choosing one of the lowest weight cuspforms worked but this might not work in general!
    """
    cuspforms_fl = dict()
    cuspforms_rig = dict()
    Kv, Kw, v_Kw, u_interior_Kv, u, lowest_non_zero_cuspform_weight = None, None, None, None, None, None
    for weight in range(2,100,2): #We only consider even weights
        dim_S = G.dimension_cusp_forms(weight)
        if dim_S != 0:
            lowest_non_zero_cuspform_weight = weight
            cuspforms_fl[weight] = get_cuspform_basis_approx(AutomorphicFormSpace(G,weight),digit_prec)
            if u == None:
                #Note that we need to be careful to select a cusp_expansion that is not an oldform!
                #For our examples we always used one of the lowest-weight cuspforms, which does however not always work in general!
                cuspform_index = -1 #The last cuspform in row-echelon form has linear term in u as first non-trivial coefficient
                Kv, Kw, v_Kw, u_interior_Kv, u = get_u_from_q_expansion(cuspforms_fl[weight][cuspform_index].get_cusp_expansion(Cusp(1,0)),dim_S+1,digit_prec,max_extension_field_degree,principal_cusp_width)
                #Now also try to recognize the second coefficient to see if we can factor out additional factors
                expression_to_recognize = cuspforms_fl[weight][cuspform_index].get_cusp_expansion(Cusp(1,0))/u**2

            cuspforms_rig[weight] = [recognize_cusp_expansion_using_u(cuspforms_fl[weight][label].get_cusp_expansion(Cusp(1,0)),weight,G,max_closed_form_trunc_order,"CuspForm",Kv,u,v_Kw,u_interior_Kv) for label in range(dim_S)]
            break
    return cuspforms_fl, cuspforms_rig, Kv, Kw, v_Kw, u_interior_Kv, u, lowest_non_zero_cuspform_weight

def get_u_from_q_expansion(cusp_expansion, coeff_index, digit_prec, max_extension_field_degree, principal_cusp_width):
    """
    Given a cusp_expansion defined over CC, try to determine u by recognizing a coefficient that is linear in u.
    """
    expression_linear_in_u = cusp_expansion[coeff_index]
    if is_effectively_zero(expression_linear_in_u,digit_prec-10) == True:
        first_non_zero_coeff_index = coeff_index+1
        while is_effectively_zero(cusp_expansion[first_non_zero_coeff_index],digit_prec-10) == True:
            first_non_zero_coeff_index += 1
        u_pow = first_non_zero_coeff_index - coeff_index + 1
        coeff_index = first_non_zero_coeff_index #We need this later in "try_to_improve_choice_of_u"
        if principal_cusp_width%u_pow == 0: #Try to recognize u from u_interior^(u_pow//principal_cusp_width)
            #We now re-use the function to try to get u_interior and Kv
            tmp = get_numberfield_of_coeff(cusp_expansion[first_non_zero_coeff_index],max_extension_field_degree,principal_cusp_width//u_pow)
            if tmp == None:
                raise ArithmeticError("Not enough precision to identify numberfield!")
            Kv, _, _, u_interior_Kv = tmp
            CC = ComplexField(digits_to_bits(digit_prec))
            minpoly = get_u_minpoly(QQbar(u_interior_Kv).nth_root(principal_cusp_width),principal_cusp_width,Kv.degree(),CC.prec())
            #We now need to find a working embedding of u
            for u in minpoly.roots(ring=QQbar,multiplicities=False):
                c = cusp_expansion[first_non_zero_coeff_index]/u**u_pow
                if to_K(c,Kv) == None: #Check if we can recognize the coefficient which indicates that u embedding is valid
                    continue
                Kw.<w> = NumberField(minpoly,embedding=u)
                v_Kw = get_v_Kw(Kv,Kw,principal_cusp_width,CC.prec())
                break
            try:
                Kw
            except NameError:
                raise ArithmeticError("Unable to find u that works!")
        else:
            raise NotImplementedError("Please only use cuspforms with non-zero coefficients to recognize u for now!")
    else:
        tmp = get_numberfield_of_coeff(expression_linear_in_u,max_extension_field_degree,principal_cusp_width)
        if tmp == None:
            raise ArithmeticError("Not enough precision to identify numberfield!")
        Kv, Kw, v_Kw, u_interior_Kv = tmp
    if principal_cusp_width == 1:
        u = convert_from_Kv_to_Kw(u_interior_Kv,v_Kw)
    else:
        u = Kw.gen()
    #Try to recognize second coeff to see if one can find a better denominator, otherwise leave parameters unchanged
    Kw, v_Kw, u_interior_Kv, u = try_to_improve_choice_of_u(cusp_expansion,coeff_index+1,Kv,Kw,v_Kw,u_interior_Kv,u,principal_cusp_width,digit_prec)
    print("u_interior_Kv: ", u_interior_Kv)
    return Kv, Kw, v_Kw, u_interior_Kv, u

def try_to_improve_choice_of_u(cusp_expansion, coeff_index, Kv, Kw, v_Kw, u_interior_Kv, u, principal_cusp_width, digit_prec):
    """
    Given an initial choice of u, try to improve u by recognizing second coefficient (i.e., a term that is quadratic in u)
    and test if this has a non-trivial denominator.
    """
    expression_to_recognize = cusp_expansion[coeff_index]/u**2
    if is_effectively_zero(expression_to_recognize,digit_prec-10) == True:
        return Kw, v_Kw, u_interior_Kv, u #Unable to improve choice
    tmp = to_K(expression_to_recognize,Kv)
    if tmp == None:
        return Kw, v_Kw, u_interior_Kv, u #Unable to improve choice
    largest_denominator = max([c.denominator() for c in list(tmp)])
    u_interior_Kv_updated = u_interior_Kv
    if largest_denominator != 1:
        largest_denominator_factored = list(factor(largest_denominator))
        for (fact,power) in largest_denominator_factored:
            if power >= 2: #Recall that our expression is quadratic in u
                u_interior_Kv_updated /= fact**principal_cusp_width
        if principal_cusp_width == 1:
            u_interior_Kv = u_interior_Kv_updated
            u = convert_from_Kv_to_Kw(u_interior_Kv,v_Kw)
        else:
            Kw, v_Kw = get_updated_Kw_v_Kw(u_interior_Kv_updated,u_interior_Kv,Kw,v_Kw,principal_cusp_width)
            u_interior_Kv = u_interior_Kv_updated
            u = Kw.gen()
    return Kw, v_Kw, u_interior_Kv, u #Return updated parameters

def get_product_formulas_for_forms(G, weight, is_cuspform):
    """
    Return products that can be used to generate forms of specified weight by multiplying forms of lower weight.
    To represent a modular form, we use tuples (in fact lists for mutation) of the form "[is_cuspform,weight,label]".
    """
    weight_combinations = get_potential_combinations_to_construct_weight(weight)
    cuspform_dims, modform_dims = get_mod_and_cusp_form_dimensions_up_to_weight(G,weight)
    max_valuation_list = [get_max_valuation(weight_combination,modform_dims,cuspform_dims,is_cuspform) for weight_combination in weight_combinations]
    max_valuation, max_valuation_product = None, None #Try to get the product of largest valuation
    for tmp in max_valuation_list:
        if tmp != None:
            if tmp[0] > max_valuation:
                max_valuation, max_valuation_product = tmp[0], tmp[1]
    if max_valuation == None:
        return None #Impossible to use any lower weight forms
    valuation_products = dict() #Stores products of forms for each valuation
    valuation_products[max_valuation] = max_valuation_product
    if is_cuspform == True:
        min_valuation = 1
    else:
        min_valuation = 0
    for valuation in range(max_valuation-1,min_valuation-1,-1): #Use lower labels to construct products of lower valuation
        new_valuation_product = deepcopy(valuation_products[valuation+1])
        lowered_valuation = False
        for i in range(len(new_valuation_product)):
            label = new_valuation_product[i][2]
            if label > 0:
                new_valuation_product[i][2] -= 1
                lowered_valuation = True
                break
        if lowered_valuation == False:
            break
        valuation_products[valuation] = new_valuation_product
    label_products = dict() #Stores product for each label
    if is_cuspform == False:
        label_products = valuation_products #For modforms, valuations at infinity are equal to labels
    else:
        for valuation in valuation_products.keys():
            label_products[valuation-1] = valuation_products[valuation]
    return label_products

def get_potential_combinations_to_construct_weight(weight):
    """
    This function returns possible combinations to construct a form of specified weight using lower (even) weights.
    Note that the weight_combinations are always sorted with decreasing order.
    """
    if weight == 2:
        return []
    if weight == 4:
        return [(2,2)]
    if weight == 6:
        return [(4,2),(2,2,2)]
    if weight == 8:
        return [(6,2),(4,4),(4,2,2),(2,2,2,2)]
    #One could write a recursive function for general weights but we have been too lazy to implement it because it was not needed yet :P
    raise NotImplementedError("We have not implemented this functions for weights larger than 8 yet!")

def get_mod_and_cusp_form_dimensions_up_to_weight(G, max_weight):
    """
    Returns dictionaries that store dimension of cuspforms and modforms up to (and including) max_weight.
    """
    cuspform_dims = dict()
    modform_dims = dict()
    for weight in range(2,max_weight+1,2):
        cuspform_dims[weight] = G.dimension_cusp_forms(weight)
        if G.dimension_eis(weight) == 0:
            modform_dims[weight] = 0 #We cannot use these for construction of modforms because they have a wrong valuation
        else:
            modform_dims[weight] = G.dimension_modular_forms(weight)
    return cuspform_dims, modform_dims

def get_max_valuation(weight_combination, modform_dims, cuspform_dims, is_cuspform):
    """
    Given a weight_combination, return the largest possible valuation of a modular (or cusp) form that
    can be achieved by multiplying forms inside weight_combination.
    This functions returns the maximal valuation as well as a list of forms that are used to generate this valuation.
    If weight_combination is invalid (meaning that there are empty dimensional spaces inside it), this function returns None.
    """
    if is_cuspform == False: #To construct modforms, we are only allowed to multiply modular forms and not cuspforms
        max_valuation = 0
        product_list = []
        for weight in weight_combination:
            if modform_dims[weight] == 0:
                return None #This weight_combination is invalid so we cannot construct any forms from it
            max_valuation += (modform_dims[weight]-1)
            product_list.append([False,weight,modform_dims[weight]-1])
    else: #To construct cuspforms we need at least one cuspform in order to get vanishing at the other cusps
        max_valuation = 0
        product_list = []
        if cuspform_dims[weight_combination[0]] == 0: #Recall that we need at least one cuspform
            return None
        max_valuation += cuspform_dims[weight_combination[0]] #Choose a cuspform for the first for in weight_combinations
        product_list.append([True,weight_combination[0],cuspform_dims[weight_combination[0]]-1])
        for i in range(1,len(weight_combination)):
            weight = weight_combination[i]
            if modform_dims[weight] == 0:
                return None #This weight_combination is invalid so we cannot construct any forms from it
            if cuspform_dims[weight] > modform_dims[weight]-1: #For some cases it would still be preferable to use modforms to get to low valuations but we ignore this here
                max_valuation += cuspform_dims[weight]
                product_list.append([True,weight,cuspform_dims[weight]-1])
            else:
                max_valuation += modform_dims[weight]-1
                product_list.append([False,weight,modform_dims[weight]-1])
    return max_valuation, product_list

def get_constructable_and_computable_labels(product_formulas, dim):
    """
    Returns labels that are constructible (i.e. that can be constructed from lower-weight forms) as well as labels that need to be computed froms scratch.
    """
    if product_formulas == None:
        constructable_labels = [] #Labels that can be constructed from lower weight forms
        computable_labels = list(range(dim)) #Labels that we need to compute from scratch
    else:
        constructable_labels, computable_labels = [], []
        for label in range(dim):
            if label in product_formulas:
                constructable_labels.append(label)
            else:
                computable_labels.append(label)
    return constructable_labels, computable_labels

def construct_form_from_product_formula(product_formula, modforms, cuspforms):
    """
    Given a product formula, construct a form.
    Note that this function can take both rigorous and floating-point results as input.
    """
    res = 1
    for (is_cuspform,weight,label) in product_formula:
        if is_cuspform == True:
            res = cuspforms[weight][label]*res
        else:
            res = modforms[weight][label]*res
    return res

def compare_results_to_numerics(G, max_weight, modforms_rig, cuspforms_rig, eis_scaling_constants, u_QQbar, numerics_digit_prec, tol, Y_fact=1):
    """
    Compare the results that have been computed through the Belyi map with values obtained from a numerical method.
    More specifically, we test the rigorous expansions (and their u-v-factorization) at infinity of all modular objects
    by comparing them to the numerical values.
    We do the same for the Eisenstein scaling coefficients which tests the floating-point expansions at all cusps.
    """
    prec_loss = 10
    M_modform_rig_max, M_cuspform_rig_max = 0, 0
    for weight in range(2,max_weight+1,2):
        if weight in modforms_rig:
            M_modform_rig_max = max([modform_rig.get_cusp_expansion(Cusp(1,0)).degree() for modform_rig in modforms_rig[weight]]+[M_modform_rig_max])
        if weight in cuspforms_rig:
            M_cuspform_rig_max = max([cuspform_rig.get_cusp_expansion(Cusp(1,0)).degree() for cuspform_rig in cuspforms_rig[weight]]+[M_cuspform_rig_max])
    M_rig_max = max(M_modform_rig_max,M_cuspform_rig_max)
    S = AutomorphicFormSpace(G,max_weight)
    M_0_numerics = get_M_0(S,numerics_digit_prec,is_cuspform=False)
    if M_0_numerics-8 < M_rig_max: #We need to increase the size of M_numerics in order to get enough values to compare
        M_0_numerics = M_rig_max+8
    RBF = RealBallField(digits_to_bits(numerics_digit_prec))
    Y = get_horo_height_arb_wrap(S,RBF,M_0_numerics,is_cuspform=False,prec_loss=prec_loss)*RBF(Y_fact)

    for weight in range(2,max_weight+1,2): #We only consider even weights
        if G.dimension_modular_forms(weight) != 0:
            if G.dimension_eis(weight) != 0: #If the eisenstein space is zero-dimensional then we could only get cuspforms which have a different normalization...
                modforms_num = get_modform_basis_approx(AutomorphicFormSpace(G,weight),numerics_digit_prec,prec_loss=prec_loss,Y=Y,M_0=M_0_numerics)
                for i in range(len(modforms_num)):
                    if do_coefficients_match_the_numerics(modforms_rig[weight][i],modforms_num[i],tol,u_QQbar) == False:
                        print("weight: ", weight)
                        print("label: ", i)
                        raise ArithmeticError("We detected a modform coefficient that does not match the numerical values!")
            if G.dimension_cusp_forms(weight) != 0:
                cuspforms_num = get_cuspform_basis_approx(AutomorphicFormSpace(G,weight),numerics_digit_prec,prec_loss=prec_loss,Y=Y,M_0=M_0_numerics)
                for i in range(len(cuspforms_num)):
                    if do_coefficients_match_the_numerics(cuspforms_rig[weight][i],cuspforms_num[i],tol,u_QQbar) == False:
                        print("weight: ", weight)
                        print("label: ", i)
                        raise ArithmeticError("We detected a cuspform coefficient that does not match the numerical values!")
                if G.dimension_eis(weight) != 0:
                    eisforms_num, eis_scaling_constant_list_num = compute_eisenstein_series(cuspforms_num,modforms_num,return_scaling_constants=True)
                    for i in range(len(eis_scaling_constant_list_num)):
                        for j in range(len(eis_scaling_constant_list_num[i])):
                            if does_result_match_numerics(eis_scaling_constants[weight][i][j],eis_scaling_constant_list_num[i][j],tol) == False:
                                print("We detected a eis_scaling_constants that does not match the numerical values!")
                                print("weight: ", weight)
                                print("i, j: ", i, j)
                                print("eis_scaling_constants[weight][i][j]: ", eis_scaling_constants[weight][i][j])
                                print("eis_scaling_constant_list_num[i][j]: ", eis_scaling_constant_list_num[i][j])
                                print("diff: ", abs(eis_scaling_constants[weight][i][j]-eis_scaling_constant_list_num[i][j]))
                                eisforms_num, eis_scaling_constant_list_num = None, None #Delete the Eisenstein results becaues they seem to be wrong
                                break
                        if eis_scaling_constant_list_num == None:
                            break

def do_coefficients_match_the_numerics(f, f_numerics, tol, u_QQbar):
    """
    Compare coefficients of rigorously computed q-expansion (at infinity) to values obtained from a numerical method.
    """
    f_expansion = f.get_cusp_expansion(Cusp(1,0))
    f_expansion_factored = f.get_cusp_expansion(Cusp(1,0),factor_into_u_v=True)
    f_numerics_expansion = f_numerics.get_cusp_expansion(Cusp(1,0))
    CC = ComplexField(2048) #It is important not to use too little precision here because the substitution of the rigorous expressions can be ill-conditioned...
    u_CC = CC(u_QQbar)
    for i in range(f_expansion.degree()+1): #We could also get a 1/q part for j_G which is however always correct
        if does_result_match_numerics(CC(f_expansion[i]),f_numerics_expansion[i],tol) == False:
            print("coeff_index: ", i)
            print("f_expansion[i]: ", f_expansion[i])
            print("CC(f_expansion[i]): ", CC(f_expansion[i]))
            print("f_numerics_expansion[i]: ", f_numerics_expansion[i])
            print("diff: ", abs(CC(f_expansion[i])-f_numerics_expansion[i]))
            return False
        #Also test that the u-v-factoriazation works correctly
        if does_result_match_numerics(f_expansion_factored[i].change_ring(CC).subs(u=u_CC),f_numerics_expansion[i],tol) == False:
            print("coeff_index: ", i)
            print("f_expansion_factored[i]: ", f_expansion_factored[i])
            print("f_expansion_factored[i].change_ring(CC).subs(u=u_CC): ", f_expansion_factored[i].change_ring(CC).subs(u=u_CC))
            print("f_numerics_expansion[i]: ", f_numerics_expansion[i])
            print("diff: ", abs(f_expansion_factored[i].change_ring(CC).subs(u=u_CC)-f_numerics_expansion[i]))
            return False

    return True

def does_result_match_numerics(res, res_num, tol):
    """
    Check if result agrees to numerical result up to tolerated precision.
    This is first done by computing abs(res-res_num).
    If this fails, we also check for the relative gap abs(res-res_num)/abs(res_num) for the case of large exponents.
    If this fails as well, we return False.
    """
    if abs(res-res_num) < tol:
        return True
    if abs(res-res_num)/abs(res_num) < tol:
        return True
    return False

def get_max_extension_field_degree(passport):
    """
    Returns the maximal degree of the extension field. This is usually given by the amount of elements of the given passport.
    The only exception for this is when there are several cusps that have equal cusp-width to the cusp at infinity.
    """
    principal_cusp_width = passport[0].cusp_width(Cusp(1,0))
    cusp_widths = passport[0].cusp_widths()
    amount_of_equal_cusp_widths = 0 #Amount of cusp-widths that are equal to principal_cusp_width
    for cusp_width in cusp_widths:
        if cusp_width == principal_cusp_width:
            amount_of_equal_cusp_widths += 1
    return amount_of_equal_cusp_widths*len(passport)

def has_equal_list_entry(list, index):
    """
    If there exists an element in list (outside index) that is equal to list[index], return True, otherwise return False.
    """
    for i in range(len(list)):
        if i != index and list[i] == list[index]:
            return True
    return False

def get_lin_u_v_term(q_expansions):
    """
    Choose one coefficient of a cuspform that is linear in u.
    """
    for weight in range(2,10,2):
        try:
            cuspform = q_expansions[weight]["cuspforms_pretty"][-1] #We are not guaranteed that this is a newform though
        except KeyError:
            continue
        i = 1
        while cuspform[i][1] == 0 and i < cuspform.degree(): 
            i += 1
        if cuspform[i][1] == 0:
            continue
        label = len(q_expansions[weight]["cuspforms_pretty"])-1
        return cuspform[i], weight, i, label
    raise ArithmeticError("Did not find suitable term linear in u.")

def get_expr_for_other_embeddings(passport, weight, i, label, digit_prec=30):
    """
    Returns floating-point approximations for expression that is linear in u-v for Galois conjugate passport elements.
    """
    expr_for_other_embeddings = []
    for G in passport[1:]:
        expr_for_other_embedding = get_cuspform_q_expansion_approx(AutomorphicFormSpace(G,weight),digit_prec,label=label).get_cusp_expansion(Cusp(1,0))[i]
        expr_for_other_embeddings.append(expr_for_other_embedding)
    return expr_for_other_embeddings

def identify_other_embeddings_and_diffs(non_zero_lin_u_expr, expr_for_other_embeddings, Kv, u_interior_Kv, principal_cusp_width, digit_prec=30):
    """
    Return different embeddings of v in the same order as the elements of 'expr_for_other_embeddings'
    as well as the corresponding differences of the floating expressions to the embedded expressions.
    """
    CC = ComplexField(digits_to_bits(digit_prec))
    current_embedding = QQbar(Kv.gen())
    remaining_embeddings = []
    for embedding in Kv.polynomial().roots(ring=QQbar,multiplicities=False):
        if embedding != current_embedding:
            remaining_embeddings.append(embedding)
    nth_power_of_embedding_expr_list = [] #List of expressions for each embedding which we use to idenfity embeddings
    for embedding in remaining_embeddings:
        nth_power_of_embedding_expr = non_zero_lin_u_expr[1].polynomial().subs(x=embedding)**principal_cusp_width * u_interior_Kv.polynomial().subs(x=embedding)
        nth_power_of_embedding_expr_list.append(nth_power_of_embedding_expr)
    ordered_embeddings = [] #List of embeddings ordered in the same way as 'expr_for_other_embeddings'
    for expr in expr_for_other_embeddings:
        nth_power_of_expr = expr**principal_cusp_width
        diffs = [abs(nth_power_of_expr-nth_power_of_embedding_expr) for nth_power_of_embedding_expr in nth_power_of_embedding_expr_list]
        min_diff = min(diffs)
        embedding_index = diffs.index(min_diff)
        ordered_embeddings.append((remaining_embeddings[embedding_index],min_diff))
    return ordered_embeddings
    
def get_all_embeddings(passport, q_expansions, Kv, u_interior_Kv, principal_cusp_width):
    """
    Identify the corresponding passport elements for each root of Kv.
    Passport elements that do not correspond to roots of Kv (in case of decaying passports) are ignored.
    """
    G = passport[0]
    res = {(str(G.permS),str(G.permR),str(G.permT)): QQbar(Kv.gen())}
    if len(passport) == 1:
        return res
    lin_u_v_term, weight, i, label = get_lin_u_v_term(q_expansions)
    expr_for_other_embeddings = get_expr_for_other_embeddings(passport,weight,i,label)
    other_embeddings_and_diffs = identify_other_embeddings_and_diffs(lin_u_v_term,expr_for_other_embeddings,Kv,u_interior_Kv,principal_cusp_width)
    
    for i in range(1,len(passport)):
        G = passport[i]
        perms = (str(G.permS),str(G.permR),str(G.permT))
        embedding, diff = other_embeddings_and_diffs[i-1]
        if diff < 1e-12 and embedding.abs() > 1e-6: #This seems to be a true embedding
            res[perms] = embedding
    embeddings = list(res.values())
    if len(embeddings) != len(set(embeddings)):
        raise ArithmeticError("We found duplicate embeddings which means that the embeddings could not be uniquely specified!")
    if len(embeddings) != Kv.degree():
        raise ArithmeticError("The amount of embeddings does not match the degree of Kv!")
    return res

def get_unresolved_passport_elements(passport, embeddings):
    """
    Get a list of passport elements whose embedding has not been specified yet.
    """
    unresolved_embeddings = []
    print("embeddings: ", embeddings)
    for G in passport:
        perms = (str(G.permS),str(G.permR),str(G.permT))
        if perms not in embeddings:
            print("perms: ", perms)
            unresolved_embeddings.append(G)
    return unresolved_embeddings

def to_JSON(passport_data, filename=None):
    """
    Convert passport to JSON format.
    If filename is specified, store the JSON data in the file.
    """
    res = {}
    indent = 4
    for k in passport_data.keys():
        if isinstance(passport_data[k], dict):
            res[str(k)] = {}
            for k2 in passport_data[k].keys():
                if isinstance(passport_data[k][k2], dict):
                    res[str(k)][str(k2)] = {}
                    for k3 in passport_data[k][k2].keys():
                        res[str(k)][str(k2)][str(k3)] = str(passport_data[k][k2][k3])
                else:
                    res[str(k)][str(k2)] = str(passport_data[k][k2])
        else:
            res[str(k)] = str(passport_data[k])
    res["is_congruence"] = passport_data["is_congruence"] #Store is_congruence as bool instead of string because it is supported by JSON
    L = passport_data["q_expansions"][4]["modforms_raw"][0].base_ring()
    res["L"] = str(L)
    for weight in res["q_expansions"]: #Remove data defined over Kw because it creates very large JSON files
        if 'hauptmodul_raw' in res["q_expansions"][weight]:
            del res["q_expansions"][weight]['hauptmodul_raw']
        if 'modforms_raw' in res["q_expansions"][weight]:
            del res["q_expansions"][weight]['modforms_raw']
        if 'cuspforms_raw' in res["q_expansions"][weight]:
            del res["q_expansions"][weight]['cuspforms_raw']
    del res["v"] #This just returns v = v and is hence obsolete
    if passport_data["G"].genus() == 0: #Store Belyi map in a prettier form
        del res["curve"]
        res["curve"] = {}
        pc = 1
        for (p,n) in passport_data["curve"]["pc_factored_raw"]:
            pc *= p**n
        p2 = 1
        for (p,n) in passport_data["curve"]["p2_factored_raw"]:
            p2 *= p**n
        p3 = 1
        for (p,n) in passport_data["curve"]["p3_factored_raw"]:
            p3 *= p**n
        res["curve"]["belyi_map_raw"] = "(" + str(p3) + ")" + "/" + "(" + str(pc) + ")" + " = 1728 + " + "(" + str(p2) + ")" + "/" + "(" + str(pc) + ")"
        pc_factored_pretty = ""
        for (p,n) in passport_data["curve"]["pc_factored_pretty"]:
            pc_factored_pretty += "(" + str(p) + ")"
            if n != 1:
                pc_factored_pretty += "^" + str(n)
            pc_factored_pretty += "*"
        pc_factored_pretty = pc_factored_pretty[:-1] #Remove last "*"
        if pc_factored_pretty == "":
            pc_factored_pretty = "1"
        p3_factored_pretty = ""
        for (p,n) in passport_data["curve"]["p3_factored_pretty"]:
            p3_factored_pretty += "(" + str(p) + ")"
            if n != 1:
                p3_factored_pretty += "^" + str(n)
            p3_factored_pretty += "*"
        p3_factored_pretty = p3_factored_pretty[:-1] #Remove last "*"
        if p3_factored_pretty == "":
            p3_factored_pretty = "1"
        p2_factored_pretty = ""
        for (p,n) in passport_data["curve"]["p2_factored_pretty"]:
            p2_factored_pretty += "(" + str(p) + ")"
            if n != 1:
                p2_factored_pretty += "^" + str(n)
            p2_factored_pretty += "*"
        p2_factored_pretty = p2_factored_pretty[:-1] #Remove last "*"
        if p2_factored_pretty == "":
            p2_factored_pretty = "1"
        res["curve"]["belyi_map_pretty"] = "(" + p3_factored_pretty + ")" + "/" + "(" + pc_factored_pretty + ")" + " = 1728 + " + "(" + p2_factored_pretty + ")" + "/" + "(" + pc_factored_pretty + ")"

    if filename is not None:
        with open(filename, 'w') as outfile:
            json.dump(res, outfile, indent=indent)
    return json.dumps(res, indent=indent)

#Inspired by: https://stackoverflow.com/questions/1364640/python-generating-python
class CodeGenerator:
    def __init__(self, tab='    '):
        self.code = ''
        self.tab = tab
        self.indent_level = 0

    def print_code(self, file=None):
        if file is None:
            print(self.code)
        else:
            with open(file, 'w') as outfile:
                outfile.write(self.code)

    def write_line(self, string):
        self.code += self.tab * self.indent_level + string + '\n'

    def indent(self):
        self.indent_level = self.indent_level + 1

    def dedent(self):
        if self.indent_level == 0:
            raise SyntaxError("internal error in code generator")
        self.indent_level = self.indent_level - 1

def coefficients_to_list(polynomial):
    """
    Convert a polynomial to a list of coefficients.
    Somehow Sage currently does not support this because ".coefficients" neglects zero coefficients,
    while ".list()" does not start at valuation 0.
    """
    if isinstance(polynomial,sage.rings.integer.Integer): #Sometimes the polynomial is an integer
        return [polynomial]
    return [polynomial[i] for i in range(polynomial.degree()+1)]

def to_sage_script(passport_data, file_path=None):
    """
    Convert passport to a sage script.
    """
    cg = CodeGenerator()
    cg.write_line("# Sage script to reproduce the database data")
    cg.write_line("")
    cg.write_line("# Define some helper functions that are needed to transform u-v-factored expansions to L")
    cg.write_line("def convert_from_K_to_L(expression_in_K, v_L):")
    cg.indent()
    cg.write_line("\"\"\"")
    cg.write_line("Given an expression in K, convert the expression efficiently to L by plugging in v(w).")
    cg.write_line("\"\"\"")
    cg.write_line("if expression_in_K == 0:")
    cg.indent()
    cg.write_line("return 0")
    cg.dedent()
    cg.write_line("coeffs = list(expression_in_K.polynomial())")
    cg.write_line("res = coeffs[-1]")
    cg.write_line("for i in range(len(coeffs)-2,-1,-1): #Horner's method")
    cg.indent()
    cg.write_line("res = res*v_L+coeffs[i]")
    cg.dedent()
    cg.write_line("return res")
    cg.dedent()
    cg.write_line("")
    cg.write_line("def transform_u_v_factored_q_expansion_to_L(q_expansion, L, v_L, u_interior_K, principal_cusp_width):")
    cg.indent()
    cg.write_line("\"\"\"")
    cg.write_line("Given a q_expansion that has coefficients of the form (expression_in_K)*u**u_pow, convert the coefficients to L.")
    cg.write_line("\"\"\"")
    cg.write_line("if principal_cusp_width == 1:")
    cg.indent()
    cg.write_line("u = L(u_interior_K)")
    cg.dedent()
    cg.write_line("else:")
    cg.indent()
    cg.write_line("u = L.gen()")
    cg.dedent()
    cg.write_line("leading_order_exponent = q_expansion.valuation()")
    cg.write_line("coeffs = list(q_expansion)")
    cg.write_line("coeffs_L = []")
    cg.write_line("for coeff in coeffs:")
    cg.indent()
    cg.write_line("if coeff == 0 or coeff == 1:")
    cg.indent()
    cg.write_line("u_pow = 0")
    cg.dedent()
    cg.write_line("else:")
    cg.indent()
    cg.write_line("u_pow = coeff.degree()")
    cg.dedent()
    cg.write_line("coeffs_L.append(convert_from_K_to_L(coeff[u_pow],v_L)*u**u_pow)")
    cg.dedent()
    cg.write_line("P = LaurentSeriesRing(L,q_expansion.variable())")
    cg.write_line("return P(coeffs_L).shift(leading_order_exponent).O(q_expansion.prec())")
    cg.dedent()
    cg.write_line("")

    cg.write_line("res = {} #The dictionary in which we store the results")
    cg.write_line("res[\"G\"] = ArithmeticSubgroup_Permutation(S2=\"{}\",S3=\"{}\")".format(passport_data["G"].S2(),passport_data["G"].S3()))
    cg.write_line("res[\"is_congruence\"] = {}".format(passport_data["is_congruence"]))
    cg.write_line("principal_cusp_width = {}".format(passport_data["G"].cusp_width(Cusp(1,0))))
    cg.write_line("res[\"monodromy_group\"] = \"{}\"".format(passport_data["monodromy_group"]))
    cg.write_line("P.<T> = PolynomialRing(QQ)")
    K = passport_data["Kv"]
    cg.write_line("K.<v> = NumberField(P({}),embedding={}+{}*1j)".format(coefficients_to_list(K.polynomial()),CC(K.gen()).real(),CC(K.gen()).imag())) #Somehow sage does not support x+y*I embedding syntax...
    cg.write_line("res[\"K\"] = K")
    cg.write_line("res[\"v\"] = K.gen()")
    cg.write_line("u_interior_K = {}".format(passport_data["u_interior_Kv"]))
    L = passport_data["q_expansions"][4]["modforms_raw"][0].base_ring()
    cg.write_line("L.<w> = NumberField(P({}),embedding={}+{}*1j) #Base ring of the q-expansions in raw format".format(coefficients_to_list(L.polynomial()),CC(L.gen()).real(),CC(L.gen()).imag())) #Somehow sage does not support x+y*I embedding syntax...
    cg.write_line("v_L = {}".format(passport_data["v_Kw"]))
    cg.write_line("")
    cg.write_line("# The pretty form of q-expansions involves powers of u, where u is given by")
    cg.write_line("u_str = \"{}\"".format(passport_data["u_str"]))
    cg.write_line("res[\"u_str\"] = u_str")
    cg.write_line("# With embedding")
    if passport_data["G"].cusp_width(Cusp(1,0)) == 1:
        cg.write_line("res[\"u\"] = QQbar({})".format(passport_data["u_str"]))
    else:
        cg.write_line("res[\"u\"] = QQbar(L.gen())")
    cg.write_line("Pu.<u> = PolynomialRing(K)")
    cg.write_line("Pq.<{}> = LaurentSeriesRing(Pu)".format(passport_data["q_expansions"][4]["modforms_raw"][0].parent().variable_name()))

    cg.write_line("")
    cg.write_line("# Store the q-expansions in pretty form into res")
    cg.write_line("res[\"q_expansions\"] = {}")
    weights = sorted(list(passport_data["q_expansions"].keys()))
    cg.write_line("weights = {}".format(weights))
    cg.write_line("for weight in weights:")
    cg.indent()
    cg.write_line("res[\"q_expansions\"][weight] = {}")
    cg.dedent()
    for weight in weights:
        if "hauptmodul_pretty" in passport_data["q_expansions"][weight]:
            cg.write_line("res[\"q_expansions\"][{}][\"hauptmodul_pretty\"] = Pq({}).O({})".format(weight, coefficients_to_list(passport_data["q_expansions"][weight]["hauptmodul_pretty"]), passport_data["q_expansions"][weight]["hauptmodul_pretty"].prec()))
            cg.write_line("res[\"q_expansions\"][{}][\"hauptmodul_pretty\"] += Pq.gen()**(-1) #Don't forget the first coeff, which we cannot construct from a list".format(weight))
        if "cuspforms_pretty" in passport_data["q_expansions"][weight]:
            cg.write_line("res[\"q_expansions\"][{}][\"cuspforms_pretty\"] = []".format(weight))
            for cuspform in passport_data["q_expansions"][weight]["cuspforms_pretty"]:
                cg.write_line("res[\"q_expansions\"][{}][\"cuspforms_pretty\"].append(Pq({}).O({}))".format(weight, coefficients_to_list(cuspform), cuspform.prec()))
        if "modforms_pretty" in passport_data["q_expansions"][weight]:
            cg.write_line("res[\"q_expansions\"][{}][\"modforms_pretty\"] = []".format(weight))
            for modform in passport_data["q_expansions"][weight]["modforms_pretty"]:
                cg.write_line("res[\"q_expansions\"][{}][\"modforms_pretty\"].append(Pq({}).O({}))".format(weight, coefficients_to_list(modform), modform.prec()))
    
    cg.write_line("")
    cg.write_line("# Now consider the q-expansions in raw format defined over L")
    cg.write_line("Pq.<{}> = LaurentSeriesRing(L)".format(passport_data["q_expansions"][4]["modforms_raw"][0].parent().variable_name()))
    cg.write_line("for weight in weights:")
    cg.indent()
    cg.write_line("if \"hauptmodul_pretty\" in res[\"q_expansions\"][weight]:")
    cg.indent()
    cg.write_line("res[\"q_expansions\"][weight][\"hauptmodul_raw\"] = transform_u_v_factored_q_expansion_to_L(res[\"q_expansions\"][weight][\"hauptmodul_pretty\"],L,v_L,u_interior_K,principal_cusp_width)")
    cg.dedent()
    cg.write_line("if \"cuspforms_pretty\" in res[\"q_expansions\"][weight]:")
    cg.indent()
    cg.write_line("res[\"q_expansions\"][weight][\"cuspforms_raw\"] = [transform_u_v_factored_q_expansion_to_L(cuspform,L,v_L,u_interior_K,principal_cusp_width) for cuspform in res[\"q_expansions\"][weight][\"cuspforms_pretty\"]]")
    cg.dedent()
    cg.write_line("if \"modforms_pretty\" in res[\"q_expansions\"][weight]:")
    cg.indent()
    cg.write_line("res[\"q_expansions\"][weight][\"modforms_raw\"] = [transform_u_v_factored_q_expansion_to_L(modform,L,v_L,u_interior_K,principal_cusp_width) for modform in res[\"q_expansions\"][weight][\"modforms_pretty\"]]")
    cg.dedent()
    cg.dedent()

    cg.write_line("")
    cg.write_line("# Add the Eisenstein scaling constants as well")
    for weight in weights:
        if "eisenstein_basis_factors" in passport_data["q_expansions"][weight]:
            cg.write_line("res[\"q_expansions\"][{}][\"eisenstein_basis_factors\"] = {}".format(weight, passport_data["q_expansions"][weight]["eisenstein_basis_factors"]))
        if "eisenstein_canonical_normalizations" in passport_data["q_expansions"][weight]:
            cg.write_line("res[\"q_expansions\"][{}][\"eisenstein_canonical_normalizations\"] = {}".format(weight, passport_data["q_expansions"][weight]["eisenstein_canonical_normalizations"]))
    cg.write_line("")
    cg.write_line("# Now add the curve")
    if passport_data["G"].genus() == 0:
        pc = 1
        for (p,n) in passport_data["curve"]["pc_factored_raw"]:
            pc *= p**n
        p3 = 1
        for (p,n) in passport_data["curve"]["p3_factored_raw"]:
            p3 *= p**n
        cg.write_line("F.<x> = FunctionField(L)")
        cg.write_line("res[\"curve\"] = (F({}))/(F({}))".format(coefficients_to_list(p3),coefficients_to_list(pc)))
    elif passport_data["G"].genus() == 1:
        if passport_data["curve"]:
            cg.write_line("res[\"curve\"] = EllipticCurve({})".format(passport_data["curve"].a_invariants()))
        else:
            cg.write_line("res[\"curve\"] = None")
    else:
        raise NotImplementedError("Not implemented for genus > 1")
    cg.write_line("")
    cg.write_line("# Now add the embeddings")
    cg.write_line("res[\"embeddings\"] = {}".format(str(passport_data["embeddings"]).replace("?","")))
    cg.print_code(file=file_path)