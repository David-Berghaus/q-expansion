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

def test():
    B = BelyiMap(Gamma0(6))
    j_G = B.get_hauptmodul_q_expansion_approx(25,50,only_principal_cusp_expansion=False)
    B_fact = B._get_B(8)
    return B._get_regularized_modular_form_q_expansion(8,j_G,B_fact)

def numberfield_example1():
    """
    Let QQ(v) be an algebraic extension field with generator v.
    Let u be an expression that is the n-th root of an element in QQ(v) (which we know explicitly).
    We are looking for a way to perform exact arithmetic with expressions in u and v and afterwards recover the results.
    This function makes use of "number_field_elements_from_algebraics" (see https://ask.sagemath.org/question/40389/extension-field-adjoining-two-roots/?answer=40411#post-id-40411).
    """
    P.<x> = QQ[]
    p = P(x**3+x**2-x+1)
    #p = P(x^7 - 2*x^6 + 468*x^5 - 4568*x^4 + 102256*x^3 - 1569312*x^2 + 19226560*x + 45652288)
    Kv.<v> = NumberField(p,embedding=CC(-0.5,0.86)) #This is QQ(v) (the embedding is arbitrary)
    u = QQbar(123*v**5+64*v**4/13+15*v**3+1769*v**2+(1763*v+1255)/7).nth_root(6) #Arbitrary expression

    K, (u_K,v_K), phi = number_field_elements_from_algebraics([u,v],minimal=True) #phi is the embedding K -> QQbar
    print("u_K: ", u_K)
    print("v_K: ", v_K)
    K.<a> = NumberField(K.polynomial(),embedding=phi(K.gen()))
    print("K: ", K) #K contains both u and v so we can use it to perform arithmetic

    #Perform some arithmetic in K:
    some_expression_from_a_computation = (15*v_K/2+7)*u_K**3 #Arbitrary expression

    #Get the expression in v back:
    expression_to_recognize_in_v = some_expression_from_a_computation/u_K**3
    #res = Kv(expression_to_recognize_in_v) #fails
    roots_Kv = expression_to_recognize_in_v.minpoly().roots(ring=Kv,multiplicities=False) #This seems like it can be expensive. Is there a better way for this?
    diffs = [(CC(root_Kv)*CC(phi(u_K))**3-CC(phi(some_expression_from_a_computation))).abs() for root_Kv in roots_Kv]
    res = roots_Kv[diffs.index(min(diffs))] #Localize the right embedding
    assert(res == 15*v/2+7)

    return res

def numberfield_example2():
    """
    Let QQ(v) be an algebraic extension field with generator v.
    Let u be an expression that is the n-th root of an element in QQ(v) (which we know explicitly).
    We are looking for a way to perform exact arithmetic with expressions in u and v and afterwards recover the results.
    This function works in the numberfield of u and does everything exact (but can be slow).
    """
    import time

    P.<x> = QQ[]
    p = P(x**3+x**2-x+1)
    #p = P(x^12 - 2*x^6 + 468*x^5 - 4568*x^4 + 102256*x^3 - 1569312*x^2 + 19226560*x + 45652288)
    Kv.<v> = NumberField(p,embedding=CC(-0.5,0.86)) #This is QQ(v)
    u = QQbar(123*v**5+64*v**4/13+15*v**3+1769*v**2+(1763*v+1255)/7).nth_root(6) #Arbitrary expression
    u_minpoly = u.minpoly() #Using lindep([1,u1**6,u**12]) should be faster but we are prefering cleaner code here
    K.<u_K> = NumberField(u_minpoly,embedding=CC(u)) #This numberfield contains both u and v, so we can use it for arithmetic

    t = time.time()
    v_K = K(v) #Maybe use lindep([CC(v),1,CC(u_K**6),CC(u_K**(6*2))])
    print("Computed K(v) in: ", time.time()-t)
    tmp = (15*v_K/2+7)*u_K**3 #Arbitrary expression
    print("tmp: ", tmp)
    tmp /= u_K**3
    
    #We now want to recover tmp in terms of v
    t = time.time()
    res = Kv(tmp) #We could compute this efficiently through substitution of u but sage already seems to be quite fast
    print("Computed Kv(tmp) in: ", time.time()-t)
    assert(res == 15*v/2+7)

    return res

def numberfield_example3():
    """
    Let QQ(v) be an algebraic extension field with generator v.
    Let u be an expression that is the n-th root of an element in QQ(v) (which we know explicitly).
    We are looking for a way to perform exact arithmetic with expressions in u and v and afterwards recover the results.
    This function works in the numberfield of u and uses LLL for operations that can be expensive (without loosing rigorousity on the final result).
    """
    from belyi.number_fields import lindep
    CC = ComplexField(30000)
    cusp_width = 6

    P.<x> = QQ[]
    p = P(x**4+x**2-x+1)
    #p = P(x^12 - 2*x^6 + 468*x^5 - 4568*x^4 + 102256*x^3 - 1569312*x^2 + 19226560*x + 45652288)
    Kv.<v> = NumberField(p,embedding=CC(-0.5,0.86)) #This is QQ(v)
    u_interior = 123*v**5+64*v**4/13+15*v**3+1769*v**2+(1763*v+1255)/7
    u_alg = QQbar(u_interior).nth_root(cusp_width) #Arbitrary expression
    u_minpoly = u_alg.minpoly() #Using lindep([1,u1**6,u**12]) should be faster but we are prefering cleaner code here
    K.<u_K> = NumberField(u_minpoly,embedding=CC(u_alg)) #This numberfield contains both u and v, so we can use it for arithmetic

    #v_K = K(v) #Can be slow
    u_K_approx_pow = CC(u_K)**cusp_width
    lll_input = [CC(v)]
    for i in range(p.degree()):
        lll_input.append(u_K_approx_pow**i)
    L = lindep(lll_input)
    if L == None:
        raise ArithmeticError("Please increase precision")
    v_K = 0
    for i in range(1,len(L)):
        pow = (i-1)*cusp_width
        v_K += L[i]*u_K**pow
    v_K /= -L[0]

    tmp = (54*v_K**3+15*v_K/2+7)*u_K**3 #Arbitrary expression
    tmp /= u_K**3
    
    #We now want to recover tmp in terms of v
    #res = Kv(tmp) #Can be slow
    x = var("x")
    coeffs = list(P(tmp.polynomial().subs(x=x**(1/cusp_width))))
    res = 0 #We could simply use .subs but the result would then be in QQbar
    for i in range(len(coeffs)):
        res += coeffs[i]*u_interior**i #Better to use Horner here

    CC = ComplexField(64)
    print("res: ", res)
    print("CC(res): ", CC(res))
    print("CC(54*v**3+15*v/2+7)", CC(54*v**3+15*v/2+7))

    assert(res == 54*v**3+15*v/2+7)

    return res

def oldform_example():
    #We consider an example of a group with a non-trivial supergroup
    #Note that this group has non-trivial conjugators between both G_prime and PSL2(Z) (see Strömberg's paper). 
    G = MySubgroup(o2='(1)(2 3)(4 5)(6 7)(8 9)',o3='(1 2 4)(3 6 8)(5 9 7)')
    print("Non-trivial supergroups: ", list(G.surgroups()))
    G_prime = list(G.surgroups())[0] #This might not always initialize a correct group...
    #Todo: Check which modular/cusp forms correspond to oldforms and how to identify them.


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
