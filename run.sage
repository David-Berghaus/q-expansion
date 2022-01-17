import numpy as np

from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from psage.groups.permutation_alg import MyPermutation
from psage.modform.arithgroup.mysubgroup import MySubgroup

from point_matching.point_matching_arb_wrap import get_coefficients_cuspform_ir_arb_wrap, get_coefficients_haupt_ir_arb_wrap, get_coefficients_modform_ir_arb_wrap, digits_to_bits, get_horo_height_arb_wrap
from classes.gamma_2_subgroup import Gamma_2_Subgroup
from classes.modform_class import ModForm
from classes.belyi_map import BelyiMap
from point_matching.fft_test import test_fft

S = AutomorphicFormSpace(Gamma0(1),weight=12)

from classes.fourier_expansion import get_cuspform_q_expansion_approx, get_modform_q_expansion_approx, get_hauptmodul_q_expansion_approx, get_cuspform_basis_approx, get_modform_basis_approx, to_reduced_row_echelon_form

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
