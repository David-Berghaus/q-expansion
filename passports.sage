from psage.modform.maass.automorphic_forms import AutomorphicFormSpace

from belyi.number_fields import is_effectively_zero
from eisenstein.eisenstein_computation import compute_eisenstein_series
from point_matching.point_matching_arb_wrap import _get_echelon_normalization_from_label
from classes.fourier_expansion import get_hauptmodul_q_expansion_approx, get_cuspform_basis_approx, get_modform_basis_approx

def compute_passport_data_genus_zero(passport, rigorous_trunc_order, eisenstein_digit_prec, max_weight, compare_result_to_numerics=True, numerics_digit_prec=30, tol=1e-10):
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
    Output:
    ------
    belyi_map: BelyiMap class instance
    j_G: Rigorous q-expansion of the hauptmodul at infinity
    cuspforms: Rigorous q-expansion of the cuspforms at infinity
    modforms: Rigorous q-expansion of the modforms at infinity
    eis_scaling_constants: Constants such that the basis of modforms multiplied with these constants constructs the space of Eisenstein series
    """
    max_extension_field_degree = get_max_extension_field_degree(passport)
    B = BelyiMap(passport[0],max_extension_field_degree=max_extension_field_degree)
    if B._Kv.degree() != max_extension_field_degree:
        raise ArithmeticError("We have not considered the case of decaying numberfields yet!")
    G = B.G

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
                        if eis_scaling_constant_list[i][j] != 0 and eis_scaling_constant_list[i][j] != 1 and is_effectively_zero(eis_scaling_constant_list[i][j],eisenstein_digit_prec) == True:
                            eis_scaling_constant_list[i][j] = 0 #We have a numerical zero which we now set to a true zero
            else:
                eisforms_fl = modforms_fl #In this case the eisforms are equivalent to modforms
                eis_scaling_constant_list = [_get_echelon_normalization_from_label(i,len(eisforms_fl)) for i in range(len(eisforms_fl))]
            eis_scaling_constants[weight] = eis_scaling_constant_list

    if compare_result_to_numerics == True:
        compare_results_to_numerics(G,max_weight,modforms_rig,cuspforms_rig,eis_scaling_constants,B.u_QQbar,numerics_digit_prec,tol)#Verify results by comparing them to numerical values

    #Should we also specify the embeddings of v into CC for different passport elements?

    #We currently don't return the arb-results because sage does not support saving these...
    res = dict()
    res["G"] = B.G.as_permutation_group()
    res["Kv"] = B._Kv
    res["v"] = B._Kv.gen()
    res["u"] = B.u_QQbar
    res["u_str"] = B.get_u_str()
    res["curve"] = B._return_res_as_dict()
    res["q_expansions"] = dict()
    for weight in range(0,max_weight+1,2): #We only consider even weights
        res["q_expansions"][weight] = dict()
        if weight == 0:
            res["q_expansions"][weight]["hauptmodul_raw"] = j_G_rig.get_cusp_expansion(Cusp(1,0))
            res["q_expansions"][weight]["hauptmodul_pretty"] = j_G_rig.get_cusp_expansion(Cusp(1,0),factor_into_u_v=True)
            #res["q_expansions"][weight]["hauptmodul_arb"] = j_G_fl.get_cusp_expansion(Cusp(1,0))
        else:
            if G.dimension_modular_forms(weight) != 0:
                res["q_expansions"][weight]["modforms_raw"] = [modform.get_cusp_expansion(Cusp(1,0)) for modform in modforms_rig[weight]]
                res["q_expansions"][weight]["modforms_pretty"] = [modform.get_cusp_expansion(Cusp(1,0),factor_into_u_v=True) for modform in modforms_rig[weight]]
                #res["q_expansions"][weight]["modforms_arb"] = [modform.get_cusp_expansion(Cusp(1,0)) for modform in modforms_fl[weight]]
                res["q_expansions"][weight]["eisenstein_basis_factors"] = eis_scaling_constants[weight]
                if G.dimension_cusp_forms(weight) != 0:
                    res["q_expansions"][weight]["cuspforms_raw"] = [cuspform.get_cusp_expansion(Cusp(1,0)) for cuspform in cuspforms_rig[weight]]
                    res["q_expansions"][weight]["cuspforms_pretty"] = [cuspform.get_cusp_expansion(Cusp(1,0),factor_into_u_v=True) for cuspform in cuspforms_rig[weight]]
                    #res["q_expansions"][weight]["cuspforms_arb"] = [cuspform.get_cusp_expansion(Cusp(1,0)) for cuspform in cuspforms_fl[weight]]
    return res

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
                        if abs(eis_scaling_constants[weight][i][j]-eis_scaling_constant_list_num[i][j]) > tol:
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
        diff = abs(f_expansion[i]-f_numerics_expansion[i])
        if diff > tol:
            return False
        #Also test that the u-v-factoriazation works correctly
        diff = abs(f_expansion_factored[i].change_ring(CC).subs(u=u_CC)-f_numerics_expansion[i])
        if diff > tol:
            return False

    return True

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
