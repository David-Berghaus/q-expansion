import numpy as np

from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from psage.groups.permutation_alg import MyPermutation
from psage.modform.arithgroup.mysubgroup import MySubgroup

from point_matching.point_matching_arb_wrap import get_coefficients_cuspform_ir_arb_wrap, get_coefficients_haupt_ir_arb_wrap, get_coefficients_modform_ir_arb_wrap, digits_to_bits, get_horo_height_arb_wrap
from classes.gamma_2_subgroup import Gamma_2_Subgroup
from classes.modform_class import ModForm
from classes.belyi_map import BelyiMap
from eisenstein.haberland import compute_petersson_product_haberland
from eisenstein.eisenstein_computation import compute_eisenstein_series

S = AutomorphicFormSpace(Gamma0(1),weight=12)

from classes.fourier_expansion import get_cuspform_q_expansion_approx, get_modform_q_expansion_approx, get_hauptmodul_q_expansion_approx, get_cuspform_basis_approx, get_modform_basis_approx, to_reduced_row_echelon_form

load("subgroups.sage")
load("passports.sage")

def eisenstein_algebraics_test(passport, floating_expansions, weight, d):
    """
    Search for algebraic relations among eisenstein series of specified weight up to order d.
    """
    from eisenstein.eisenstein_computation import echelon_basis_to_eisenstein_basis, get_algebraic_relations
    modforms = floating_expansions[weight]["modforms_float"]
    eis_facts = passport["q_expansions"][weight]["eisenstein_basis_factors"]
    eis_forms_echelon = []
    for i in range(len(eis_facts)):
        eis_form_echelon = 0
        for j in range(len(modforms)):
            eis_form_echelon = modforms[j]*eis_facts[i][j] + eis_form_echelon
        eis_forms_echelon.append(eis_form_echelon)
    eis_forms = echelon_basis_to_eisenstein_basis(eis_forms_echelon)
    print("algebraic relations: ", get_algebraic_relations(eis_forms,d))
    return eis_forms

def haberland_precision_test():
    """
    Test precision of Petersson product evaluation for several examples.
    """
    passports = get_list_of_all_passports(11,0)
    for i in range(len(passports)):
        G = passports[i][0]
        B = BelyiMap(G)
        digit_prec1 = 150
        digit_prec2 = 170
        trunc_orders1 = B._get_trunc_orders_convergence(4,digit_prec1)
        trunc_orders2 = B._get_trunc_orders_convergence(4,digit_prec2)
        jG1 = B.get_hauptmodul_q_expansion_approx(trunc_orders1,digit_prec1)
        jG2 = B.get_hauptmodul_q_expansion_approx(trunc_orders2,digit_prec2)
        cs1 = B.get_cuspforms(4,trunc_orders1,j_G=jG1)
        cs2 = B.get_cuspforms(4,trunc_orders2,j_G=jG2)
        ms1 = B.get_modforms(4,trunc_orders1,j_G=jG1)
        ms2 = B.get_modforms(4,trunc_orders2,j_G=jG2)
        for ic in range(len(cs1)):
            cs1[ic]._set_constant_coefficients_to_zero_inplace()
            cs2[ic]._set_constant_coefficients_to_zero_inplace()
            c1 = cs1[ic]
            c2 = cs2[ic]
            for im in range(len(ms1)):
                m1 = ms1[im]
                m2 = ms2[im]
                p1 = compute_petersson_product_haberland(c1,m1)
                p2 = compute_petersson_product_haberland(c2,m2)
                print((p2-p1).abs())

def test_construction_of_higher_weight_forms():
    B = BelyiMap(U1)
    cuspforms = dict()
    modforms = dict()
    for weight in range(2,7,2):
        if U1.dimension_cusp_forms(weight) != 0:
            cuspforms[weight] = B.get_cuspforms(weight,10)
        if U1.dimension_modular_forms(weight) != 0:
            modforms[weight] = B.get_modforms(weight,10)
    # test = [modforms[2][0]**2,modforms[4][1],modforms[4][2]]
    # test = to_reduced_row_echelon_form(test)
    # print("modforms[4]: ", modforms[4])
    test = [cuspforms[4][0]*modforms[2][0],cuspforms[6][1]]
    test = to_reduced_row_echelon_form(test)
    print("cuspforms[6]: ", cuspforms[6])
    print("test: ", test)

def test(trunc_order, digit_prec):
    bit_prec = round(3.33*digit_prec)
    CBF = ComplexBallField(bit_prec)
    P.<x> = PowerSeriesRing(CBF)
    s = P(x).O(trunc_order)/P((x-16)**3).O(trunc_order)
    print("s[s.degree()]: ", s[s.degree()])
    r = s.reverse()
    print("r[r.degree()]: ", r[r.degree()])
    res = r.subs({x:1/j_invariant_qexp(trunc_order)})
    return res

def test2(G):
    B = BelyiMap(G)
    trunc_orders = B._get_trunc_orders_convergence(6,100)
    cs = B.get_cuspforms(6,trunc_orders,digit_prec=100)
    ms = B.get_modforms(6,trunc_orders,digit_prec=100)
    from eisenstein.eisenstein_computation import compute_eisenstein_series
    es, sc = compute_eisenstein_series(cs,ms,return_scaling_constants=True)
    return es, sc

#These notations are part of the Fiori, Franc paper. Note that we do not get the same permT...
G1 = MySubgroup(o2='(1 2)(3 4)(5 6)(7)',o3='(1)(2 3 5)(4 6 7)')
u_G1 = (-7)**(1/4)/7**2
H1 = MySubgroup(o2='(1 2)(3 4)(5 6)(7)',o3='(1)(2 3 5)(7 6 4)')
u_H1 = (-7**3)**(1/5)/7**2
# z3 = exp(2*pi*I/3)
# u_U1 = ((1763*z3 + 1255)*2**2*3/7**7)**(1/6)
U1 = MySubgroup(o2='(1 2)(3 4)(5)(6 7)',o3='(1)(2 3 5)(4 6 7)')

def oldform_example():
    #We consider an example of a group with a non-trivial supergroup
    #Note that this group has non-trivial conjugators between both G_prime and PSL2(Z) (see Strömberg's paper). 
    G = MySubgroup(o2='(1)(2 3)(4 5)(6 7)(8 9)',o3='(1 2 4)(3 6 8)(5 9 7)')
    print("Non-trivial supergroups: ", list(G.surgroups()))
    G_prime = list(G.surgroups())[0] #This might not always initialize a correct group...
    #Todo: Check which modular/cusp forms correspond to oldforms and how to identify them.

def has_subgroup_duplicate_cusps(G):
    """
    Check if subgroup has cusps with equal cusp-width.
    """
    cusp_widths = G.cusp_widths()
    return len(cusp_widths) != len(set(cusp_widths))

def numberfield_reduction_example():
    #Example how to reduce a numberfield and write an algebraic number in terms of powers of this reduced numberfield
    P.<x> = ZZ[]

    numberfield = P(17*x**5 + 124*x**4 + 69*x**3 + 420*x**2 + 3*x + 1) #Some random polynomial to define a numberfield
    r = numberfield.roots(ring=QQbar,multiplicities=False)
    numberfield_red = P(pari.polredabs(numberfield))
    r_red = numberfield_red.roots(ring=QQbar,multiplicities=False)

    gp("default(realprecision, 100)")
    alg_number_approximation = N(r[1],digits=100) #Obviously we usually start here and then determine 'numberfield'
    for i in range(len(r_red)):
        gp("x = " + N(r_red[i],digits=100).str())
        gp("alg_number_approximation = " + alg_number_approximation.str())
        gp_command = "lindep([alg_number_approximation,1"
        for j in range(1,5):
            gp_command += ",x^" + str(j)
        gp_command += ",Pi])"
        lindep_res = gp(gp_command).sage()
        print(r_red[i], lindep_res)
    print("We hence see that our algebraic number can be written as a powers of ", r_red[4])
