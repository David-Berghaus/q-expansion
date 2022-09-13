import numpy as np

from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from psage.groups.permutation_alg import MyPermutation
from psage.modform.arithgroup.mysubgroup import MySubgroup

from point_matching.point_matching_arb_wrap import get_coefficients_cuspform_ir_arb_wrap, get_coefficients_haupt_ir_arb_wrap, get_coefficients_modform_ir_arb_wrap, digits_to_bits, get_horo_height_arb_wrap, get_coefficients_gmres_modform_arb_wrap, get_modform_basis_gmres_arb_wrap
from classes.gamma_2_subgroup import Gamma_2_Subgroup
from classes.modform_class import ModForm
from classes.belyi_map import BelyiMap
from eisenstein.haberland import compute_petersson_product_haberland
from eisenstein.eisenstein_computation import compute_eisenstein_series

S = AutomorphicFormSpace(Gamma0(1),weight=12)

H_5 = MySubgroup(o2='(1 2)(3 4)(5 8)(6 9)(7 12)(10 14)(11 15)(13 17)(16 19)(18 21)(20 24)(22 28)(23 30)(25 32)(26 33)(27 36)(29 37)(31 38)(34 40)(35 41)(39 43)(42 45)(44 49)(46 51)(47 53)(48 55)(50 57)(52 58)(54 60)(56 59)',
    o3='(1 3 2)(4 5 9)(6 13 14)(7 15 8)(10 18 17)(11 19 12)(16 23 24)(20 31 32)(21 22 33)(25 39 38)(26 34 40)(27 41 28)(29 37 30)(35 45 36)(42 50 51)(43 44 53)(46 52 58)(47 59 60)(48 55 49)(54 57 56)')

from classes.fourier_expansion import get_cuspform_q_expansion_approx, get_modform_q_expansion_approx, get_hauptmodul_q_expansion_approx, get_cuspform_basis_approx, get_modform_basis_approx, to_reduced_row_echelon_form

load("subgroups.sage")
load("passports.sage")

def get_M_2_H_5(digit_prec, maxiter):
    from classes.fourier_expansion import FourierExpansion, c_vec_to_cusp_expansions
    from point_matching.point_matching_arb_wrap import _get_normalization_modforms
    
    S = AutomorphicFormSpace(H_5,2)

    if digit_prec < 60:
        raise ArithmeticError("Digit precision must be at least 60.")
    c_vec_m_9 = get_coefficients_gmres_modform_arb_wrap(S,digit_prec,multiplicity=11,label=10,maxiter=maxiter)
    c_vecs, M_0, labels = get_modform_basis_gmres_arb_wrap(S,digit_prec,multiplicity=9,labels=[i for i in range(9)],normalization_zeros=[1],return_M_and_labels=True,maxiter=maxiter)

    #Transform c_vecs for m_0-m_8 into FourierExpansion objects
    starting_order = 0
    bit_prec = digits_to_bits(digit_prec)
    basis = []
    for i in range(len(c_vecs)):
        normalization = _get_normalization_modforms(S,9,label=labels[i])
        c_vec_mcbd = c_vecs[i]._get_mcbd(bit_prec)
        cusp_expansions = c_vec_to_cusp_expansions(c_vec_mcbd,S,starting_order,normalization,M_0)
        base_ring = cusp_expansions[Cusp(1,0)].base_ring()
        basis.append(FourierExpansion(S.group(),S.weight(),cusp_expansions,"ModForm",base_ring))

    #Now for the m_9 form
    normalization = _get_normalization_modforms(S,11,label=10)
    c_vec_mcbd = c_vec_m_9._get_mcbd(bit_prec)
    cusp_expansions = c_vec_to_cusp_expansions(c_vec_mcbd,S,starting_order,normalization,M_0)
    base_ring = cusp_expansions[Cusp(1,0)].base_ring()
    basis.append(FourierExpansion(S.group(),S.weight(),cusp_expansions,"ModForm",base_ring))

    return basis

def test():
    # for d in range(100,120):
    #     try:
    #         print("Testing d = ", d)
    #         tmp = get_coefficients_modform_ir_arb_wrap(AutomorphicFormSpace(H_5,2),40,label=d-1,multiplicity=d,maxiter=10)
    #         return tmp
    #     except ArithmeticError: #This happens if we exceed maxiter
    #         continue
    for d in range(49,50):
        try:
            print("Testing d = ", d)
            tmp = get_coefficients_modform_ir_arb_wrap(AutomorphicFormSpace(Gamma0(28),12),50,label=d-1,multiplicity=d,maxiter=10)
            return tmp
        except ArithmeticError: #This happens if we exceed maxiter
            continue

def test_eig():
    """
    Try to compute the space S_4(H_5) by computing the eigenspace of eigenvalue 1 associated to A*b=b.
    """
    from numpy.linalg import eig
    import math
    from point_matching.point_matching_dp import get_V_matrix_cuspform_dp

    G = H_5
    M_0 = 123
    Y = 0.06
    nc = G.ncusps()
    weight = 4

    # G = Gamma0(1)
    # M_0 = 30
    # Y = 0.5
    # nc = G.ncusps()
    # weight = 12

    V = get_V_matrix_cuspform_dp(AutomorphicFormSpace(G,weight),M_0,Y)

    #Scale linear system of equations so that the Fourier coefficients correspond to the lambda=1 eigenspace.
    for cii in range(nc):
        for cjj in range(nc):
            V_view = V[cii*M_0:(cii+1)*M_0,cjj*M_0:(cjj+1)*M_0] #using a memory view is certainly not efficient here
            for n in range(1,M_0+1):
                V_view[n-1,:] /= Y**(weight//2)*math.exp(-2*math.pi*Y*n)
    
    #Compute eigenvalues and eigenvectors
    eigvals, eigvecs = eig(V)

    #Sort eigenvectors by eigenvalue size
    eigvals = np.abs(eigvals)
    idx = eigvals.argsort()[::-1]
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:,idx]

    return eigvals, eigvecs

    
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