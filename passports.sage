load("run.sage")

from belyi.number_fields import is_effectively_zero
from eisenstein.eisenstein_computation import compute_eisenstein_series
from point_matching.point_matching_arb_wrap import _get_echelon_normalization_from_label

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

def compute_passport_data_genus_zero(passport, rigorous_trunc_order, eisenstein_digit_prec, max_weight):
    """
    Compute relevant data for a specified passport.
    Input:
    ------
    passport: List of (Galois conjugate) subgroups that are part of the passport.
    rigorous_trunc_order: Amount of terms in the rigorous computation of the q-expansion at infinity
    eisenstein_digit_prec: Digit precision of the numerical approximation of the Eisenstein series
    max_weight: The maximum weight of the computed modular forms
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
    eis_scaling_constants = dict()
    for weight in range(2,max_weight+1,2): #We only consider even weights
        if G.dimension_modular_forms(weight) != 0:
            modforms_fl = B.get_modforms(weight,eisenstein_trunc_orders,digit_prec=eisenstein_digit_prec,j_G=j_G_fl)
            if G.dimension_cusp_forms(weight) != 0:
                cuspforms_fl = B.get_cuspforms(weight,eisenstein_trunc_orders,digit_prec=eisenstein_digit_prec,j_G=j_G_fl)
                eisforms_fl, eis_scaling_constant_list = compute_eisenstein_series(cuspforms_fl,modforms_fl,return_scaling_constants=True)
                for i in range(len(eis_scaling_constant_list)):
                    for j in range(len(eis_scaling_constant_list[i])):
                        if eis_scaling_constant_list[i][j] != 0 and eis_scaling_constant_list[i][j] != 1 and is_effectively_zero(eis_scaling_constant_list[i][j],eisenstein_digit_prec) == True:
                            eis_scaling_constant_list[i][j] = 0 #We have a numerical zero which we now set to a true zero
            else:
                eisforms_fl = modforms_fl #In this case the eisforms are equivalent to modforms
                eis_scaling_constant_list = [_get_echelon_normalization_from_label(i,len(eisforms_fl)) for i in range(len(eisforms_fl))]
            eis_scaling_constants[weight] = eis_scaling_constant_list

    #To Do: Embeddings

    #To Do: Verify results by comparing them to numerical values

    #To Do: Write unit test for this function

    res = dict()
    res["belyi_map"] = B
    res["j_G"] = j_G_rig
    res["cuspforms"] = cuspforms_rig
    res["modforms"] = modforms_rig
    res["eis_scaling_constants"] = eis_scaling_constants
    return res