from psage.modform.maass.automorphic_forms import AutomorphicFormSpace

from belyi.number_fields import is_effectively_zero, get_numberfield_of_coeff, to_K
from belyi.expression_in_u_and_v import convert_from_Kv_to_Kw
from eisenstein.eisenstein_computation import compute_eisenstein_series
from point_matching.point_matching_arb_wrap import _get_echelon_normalization_from_label, digits_to_bits
from classes.fourier_expansion import get_hauptmodul_q_expansion_approx, get_cuspform_basis_approx, get_modform_basis_approx, recognize_cusp_expansion_using_u
from classes.belyi_map import BelyiMap

def compute_passport_data_genus_zero(passport, rigorous_trunc_order, eisenstein_digit_prec, max_weight, return_newton_res=False, compare_result_to_numerics=True, numerics_digit_prec=30, tol=1e-10):
    """
    Compute relevant data for a specified passport.
    Input:
    ------
    passport: List of (Galois conjugate) subgroups that are part of the passport.
    rigorous_trunc_order: Amount of terms in the rigorous computation of the q-expansion at infinity
    eisenstein_digit_prec: Digit precision of the numerical approximation of the Eisenstein series
    max_weight: The maximum weight of the computed modular forms
    return_newton_res: If true, return object that allows for an efficient reconstruction of a new Belyi map instance
    compare_result_to_numerics: Boolean that decides if the results that have been computed through the Belyi map should be compared to the numerical values
    numerics_digit_prec: The precision at which the numerical computation that we use to compare the results is performed
    tol: Maximum difference between q-expansion coefficients compared to the numerical results
    Output:
    ------
    A dictionary with the computed information that can be read using Sage only.
    """
    max_extension_field_degree = get_max_extension_field_degree(passport)
    B = BelyiMap(passport[0],max_extension_field_degree=max_extension_field_degree)
    G = B.G
    if B._Kv.degree() != max_extension_field_degree:
        if has_equal_list_entry(G.cusp_widths(),0) == False: #If two cusps are identical it sometimes happens that they are in the same numberfield which we do not need to investigate further
            raise ArithmeticError("We have not considered the case of decaying numberfields yet!")

    #First do the rigorous computation of the q-expansions
    j_G_rig = B.get_hauptmodul_q_expansion(rigorous_trunc_order)
    cuspforms_rig, modforms_rig = dict(), dict()
    for weight in range(2,max_weight+1,2): #We only consider even weights
        if G.dimension_modular_forms(weight) != 0:
            modforms_rig[weight] = B.get_modforms(weight,rigorous_trunc_order,j_G=j_G_rig)
        if G.dimension_cusp_forms(weight) != 0:
            cuspforms_rig[weight] = B.get_cuspforms(weight,rigorous_trunc_order,j_G=j_G_rig)

    #Now to the numerical stuff
    eisenstein_trunc_orders = B._get_trunc_orders_convergence(max_weight,eisenstein_digit_prec)
    j_G_fl = B.get_hauptmodul_q_expansion_approx(eisenstein_trunc_orders,eisenstein_digit_prec)
    cuspforms_fl, modforms_fl = dict(), dict()
    eis_scaling_constants = dict()
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
                eisforms_fl = modforms_fl #In this case the eisforms are equivalent to modforms
                eis_scaling_constant_list = [_get_echelon_normalization_from_label(i,len(eisforms_fl)) for i in range(len(eisforms_fl))]
            eis_scaling_constants[weight] = eis_scaling_constant_list

    if compare_result_to_numerics == True:
        compare_results_to_numerics(G,max_weight,modforms_rig,cuspforms_rig,eis_scaling_constants,B.u_QQbar,numerics_digit_prec,tol)#Verify results by comparing them to numerical values

    #Should we also specify the embeddings of v into CC for different passport elements?

    res = dict()
    res["G"] = B.G.as_permutation_group()
    res["Kv"] = B._Kv
    res["v"] = B._Kv.gen()
    res["u"] = B.u_QQbar
    res["u_str"] = B.get_u_str()
    res["curve"] = B._return_res_as_dict()
    res["q_expansions"] = dict()
    CIF = ComplexIntervalField(digits_to_bits(101)) #Unfortunately arbs currently cannot be stored, see: https://trac.sagemath.org/ticket/33310#ticket
    #We limit the floating-point precision to 100 digits for storage-space reasons
    for weight in range(0,max_weight+1,2): #We only consider even weights
        res["q_expansions"][weight] = dict()
        if weight == 0:
            res["q_expansions"][weight]["hauptmodul_raw"] = j_G_rig.get_cusp_expansion(Cusp(1,0))
            res["q_expansions"][weight]["hauptmodul_pretty"] = j_G_rig.get_cusp_expansion(Cusp(1,0),factor_into_u_v=True)
            res["q_expansions"][weight]["hauptmodul_float"] = j_G_fl.get_cusp_expansion(Cusp(1,0)).change_ring(CIF)
        else:
            if G.dimension_modular_forms(weight) != 0:
                res["q_expansions"][weight]["modforms_raw"] = [modform.get_cusp_expansion(Cusp(1,0)) for modform in modforms_rig[weight]]
                res["q_expansions"][weight]["modforms_pretty"] = [modform.get_cusp_expansion(Cusp(1,0),factor_into_u_v=True) for modform in modforms_rig[weight]]
                res["q_expansions"][weight]["modforms_float"] = [modform.get_cusp_expansion(Cusp(1,0)).change_ring(CIF) for modform in modforms_fl[weight]]
                res["q_expansions"][weight]["eisenstein_basis_factors"] = eis_scaling_constants[weight]
                if G.dimension_cusp_forms(weight) != 0:
                    res["q_expansions"][weight]["cuspforms_raw"] = [cuspform.get_cusp_expansion(Cusp(1,0)) for cuspform in cuspforms_rig[weight]]
                    res["q_expansions"][weight]["cuspforms_pretty"] = [cuspform.get_cusp_expansion(Cusp(1,0),factor_into_u_v=True) for cuspform in cuspforms_rig[weight]]
                    res["q_expansions"][weight]["cuspforms_float"] = [cuspform.get_cusp_expansion(Cusp(1,0)).change_ring(CIF) for cuspform in cuspforms_fl[weight]]
    if return_newton_res == True:
        return res, B._return_newton_res()
    return res

def compute_passport_data_higher_genera(passport, max_rigorous_trunc_order, eisenstein_digit_prec, max_weight):
    #To Do:
    #-Generate (some of) higher forms by multiplying lower ones
    max_extension_field_degree = get_max_extension_field_degree(passport)
    G = passport[0]
    CC = ComplexField(digits_to_bits(eisenstein_digit_prec))
    principal_cusp_width = G.cusp_width(Cusp(1,0))

    cuspforms_fl = dict()
    cuspforms_rig = dict()
    Kv, Kw, v_Kw, u_interior_Kv, u = None, None, None, None, None
    for weight in range(2,max_weight+1,2): #We only consider even weights
        dim_S = G.dimension_cusp_forms(weight)
        if dim_S != 0:
            cuspforms_fl[weight] = get_cuspform_basis_approx(AutomorphicFormSpace(G,weight),eisenstein_digit_prec)
            if u == None:
                #Note that we need to be careful to select a cusp_expansion that is not an oldform!
                #For our examples we always used one of the lowest-weight cuspforms, which does however not always work in general!
                cuspform_index = -1 #The last cuspform in row-echelon form has linear term in u as first non-trivial coefficient
                Kv, Kw, v_Kw, u_interior_Kv, u = get_u_from_q_expansion(cuspforms_fl[weight][cuspform_index].get_cusp_expansion(Cusp(1,0)),dim_S+1,eisenstein_digit_prec,max_extension_field_degree,principal_cusp_width)
                if Kv.degree() != max_extension_field_degree:
                    if has_equal_list_entry(G.cusp_widths(),0) == False: #If two cusps are identical it sometimes happens that they are in the same numberfield which we do not need to investigate further
                        raise ArithmeticError("We have not considered the case of decaying numberfields yet! Please also make sure that the selected cusp_expansion is not an oldform!")
            cuspforms_rig[weight] = [recognize_cusp_expansion_using_u(cuspforms_fl[weight][label].get_cusp_expansion(Cusp(1,0)),weight,G,max_rigorous_trunc_order,"CuspForm",label,Kv,u,v_Kw,u_interior_Kv) for label in range(dim_S)]

    modforms_fl = dict()
    modforms_rig = dict()
    for weight in range(2,max_weight+1,2): #We only consider even weights
        dim_M = G.dimension_modular_forms(weight)
        if dim_M != 0:
            modforms_fl[weight] = get_modform_basis_approx(AutomorphicFormSpace(G,weight),eisenstein_digit_prec)
            modforms_rig[weight] = [recognize_cusp_expansion_using_u(modforms_fl[weight][label].get_cusp_expansion(Cusp(1,0)),weight,G,max_rigorous_trunc_order,"ModForm",label,Kv,u,v_Kw,u_interior_Kv) for label in range(dim_M)]

    return cuspforms_rig, modforms_rig

def get_u_from_q_expansion(cusp_expansion, coeff_index, digit_prec, max_extension_field_degree, principal_cusp_width):
    """
    Given a cusp_expansion defined over CC, try to determine u by recognizing a coefficient that is linear in u.
    """
    expression_linear_in_u = cusp_expansion[coeff_index]
    if is_effectively_zero(expression_linear_in_u,digit_prec-5) == True:
        raise NotImplementedError("Please only use cuspforms with non-zero coefficients to recognize u for now!")
    tmp = get_numberfield_of_coeff(expression_linear_in_u,max_extension_field_degree,principal_cusp_width)
    if tmp == None:
        raise ArithmeticError("Not enough precision to identify numberfield!")
    Kv, Kw, v_Kw, u_interior_Kv = tmp
    print("u_interior_Kv: ", u_interior_Kv)
    if principal_cusp_width == 1:
        u = convert_from_Kv_to_Kw(u_interior_Kv,v_Kw)
    else:
        u = Kw.gen()
    return Kv, Kw, v_Kw, u_interior_Kv, u

def compare_results_to_numerics(G, max_weight, modforms_rig, cuspforms_rig, eis_scaling_constants, u_QQbar, numerics_digit_prec, tol):
    """
    Compare the results that have been computed through the Belyi map with values obtained from a numerical method.
    More specifically, we test the rigorous expansions (and their u-v-factorization) at infinity of all modular objects
    by comparing them to the numerical values.
    We do the same for the Eisenstein scaling coefficients which tests the floating-point expansions at all cusps.
    """
    for weight in range(2,max_weight+1,2): #We only consider even weights
        if G.dimension_modular_forms(weight) != 0:
            modforms_num = get_modform_basis_approx(AutomorphicFormSpace(G,weight),numerics_digit_prec,prec_loss=10)
            for i in range(len(modforms_num)):
                if do_coefficients_match_the_numerics(modforms_rig[weight][i],modforms_num[i],tol,u_QQbar) == False:
                    raise ArithmeticError("We detected a modform coefficient that does not match the numerical values!")
            if G.dimension_cusp_forms(weight) != 0:
                cuspforms_num = get_cuspform_basis_approx(AutomorphicFormSpace(G,weight),numerics_digit_prec,prec_loss=10)
                for i in range(len(cuspforms_num)):
                    if do_coefficients_match_the_numerics(cuspforms_rig[weight][i],cuspforms_num[i],tol,u_QQbar) == False:
                        raise ArithmeticError("We detected a cuspform coefficient that does not match the numerical values!")
                eisforms_num, eis_scaling_constant_list_num = compute_eisenstein_series(cuspforms_num,modforms_num,return_scaling_constants=True)
                for i in range(len(eis_scaling_constant_list_num)):
                    for j in range(len(eis_scaling_constant_list_num[i])):
                        if does_result_match_numerics(eis_scaling_constants[weight][i][j],eis_scaling_constant_list_num[i][j],tol) == False:
                            print("i, j: ", i, j)
                            print("eis_scaling_constants[weight][i][j]: ", eis_scaling_constants[weight][i][j])
                            print("eis_scaling_constant_list_num[i][j]: ", eis_scaling_constant_list_num[i][j])
                            print("diff: ", abs(eis_scaling_constants[weight][i][j]-eis_scaling_constant_list_num[i][j]))
                            raise ArithmeticError("We detected a eis_scaling_constants that does not match the numerical values!")

def do_coefficients_match_the_numerics(f, f_numerics, tol, u_QQbar):
    """
    Compare coefficients of rigorously computed q-expansion (at infinity) to values obtained from a numerical method.
    """
    f_expansion = f.get_cusp_expansion(Cusp(1,0))
    f_expansion_factored = f.get_cusp_expansion(Cusp(1,0),factor_into_u_v=True)
    f_numerics_expansion = f_numerics.get_cusp_expansion(Cusp(1,0))
    CC = f_numerics_expansion.base_ring()
    u_CC = CC(u_QQbar)
    for i in range(f_expansion.degree()+1): #We could also get a 1/q part for j_G which is however always correct
        if does_result_match_numerics(f_expansion[i],f_numerics_expansion[i],tol) == False:
            print("i: ", i)
            print("CC(f_expansion[i]): ", CC(f_expansion[i]))
            print("f_numerics_expansion[i]: ", f_numerics_expansion[i])
            print("diff: ", abs(f_expansion[i]-f_numerics_expansion[i]))
            return False
        #Also test that the u-v-factoriazation works correctly
        if does_result_match_numerics(f_expansion_factored[i].change_ring(CC).subs(u=u_CC),f_numerics_expansion[i],tol) == False:
            print("i: ", i)
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
