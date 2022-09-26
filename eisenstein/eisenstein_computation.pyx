# ****************************************************************************
#       Copyright (C) 2022 David Berghaus <berghaus@th.physik.uni-bonn.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from itertools import combinations_with_replacement

from sage.matrix.matrix_space import MatrixSpace
from sage.misc.misc_c import prod
from sage.functions.other import conjugate
from sage.symbolic.constants import pi
from sage.modular.cusps import Cusp

from point_matching.point_matching_arb_wrap import bits_to_digits
from eisenstein.haberland import compute_petersson_product_haberland, get_exp_two_pi_i_a_m_dict, clear_memoized_caches
from belyi.number_fields import lindep, is_effectively_zero

def compute_eisenstein_series(cuspforms, modforms, return_scaling_constants=False, truncate_cusp_expansions_to_convergence=True):
    """
    Compute the orthogonal complement of the cuspforms in the space of modular forms (which corresponds to a basis of Eisenstein series).
    We assume that cuspforms/modforms are lists of basis functions in reduced row echelon form.
    This function then returns a basis of Eisenstein series in reduced row echelon form.
    """
    cuspforms_CC = [cuspform._convert_to_CC() for cuspform in cuspforms] #We need to work over ComplexFields because CBFs have hashing issues in haberland
    modforms_CC = [modform._convert_to_CC() for modform in modforms]
    if truncate_cusp_expansions_to_convergence == True: #Ignore all higher coeffs that are not required to reach convergence
        cuspforms_CC = [cuspform._get_convergence_truncated_instance() for cuspform in cuspforms_CC]
        modforms_CC = [modform._get_convergence_truncated_instance() for modform in modforms_CC]
    dim_S = len(cuspforms_CC)
    dim_M = len(modforms)
    dim_E = dim_M-dim_S

    CC = cuspforms_CC[0].base_ring
    rho = (CC(0,2*pi)/3).exp()
    a_values = [rho+1,CC(0,1),CC(1,1)]
    exp_two_pi_i_a_m_dict = get_exp_two_pi_i_a_m_dict(a_values,cuspforms_CC[0],modforms_CC[0],CC) #We precompute this dict to re-use it
    petersson_products = dict()
    for i in range(dim_S):
        petersson_products[i] = [conjugate(compute_petersson_product_haberland(cuspforms_CC[i],modforms_CC[j],clear_memoized_caches_bool=False,exp_two_pi_i_a_m_dict=exp_two_pi_i_a_m_dict)) for j in range(dim_M)]
    clear_memoized_caches() #Although we could in principle think about using these for different weights

    M_A, M_b = MatrixSpace(CC,dim_S,dim_S), MatrixSpace(CC,dim_S,1)
    A = M_A([petersson_products[i][j] for i in range(dim_S) for j in range(dim_E,dim_M)])
    b_vecs = [-M_b([petersson_products[i][j] for i in range(dim_S)]) for j in range(dim_E)]
    c_vecs = [A.solve_right(b_vecs[j]) for j in range(dim_E)]

    #We are now ready to construct the eisforms from the modforms
    eisforms = []
    scaling_constants = dict() #Constants by which the modform basis gets scaled to produce eisenstein series
    for j in range(dim_E):
        eisform = modforms_CC[j]
        normalization = [1 if i == j else 0 for i in range(dim_E)]
        scaling_constants[j] = normalization
        for i in range(dim_S):
            scaling_constant = c_vecs[j][i,0]
            eisform += modforms_CC[dim_E+i]*scaling_constant
            scaling_constants[j].append(scaling_constant)
        eisforms.append(eisform)

    #The code below is in principle more generic, but might choose different pivots for different precisions which makes it harder to compare results
    # digit_prec = bits_to_digits(CC.precision())
    # M_A = MatrixSpace(CC,dim_S,dim_M)
    # A = M_A([petersson_products[i] for i in range(dim_S)])
    # for i in range(A.nrows()):
    #     for j in range(A.ncols()):
    #         if is_effectively_zero(A[i,j].abs(),digit_prec):
    #             A[i,j] = CC(0) #Make sure we set numerical zeros to be effictively zero to perform select the right pivots
    # K = A.right_kernel_matrix()
    # c_vecs = [v for v in K]

    #We are now ready to construct the eisforms from the modforms
    # eisforms = []
    # scaling_constants = dict() #Constants by which the modform basis gets scaled to produce eisenstein series
    # for j in range(dim_E):
    #     scaling_constants[j] = [c_vecs[j][i] for i in range(dim_M)]
    #     eisform = modforms_CC[0]*scaling_constants[j][0]
    #     for i in range(1,dim_M):
    #         eisform += modforms_CC[i]*scaling_constants[j][i]
    #     eisforms.append(eisform)

    if return_scaling_constants == False:
        return eisforms
    else:
        return eisforms, scaling_constants

def echelon_basis_to_eisenstein_basis(eisforms, return_scaling_constants=False):
    """
    Given a basis of Eisenstein series in reduced row echelon form, transform basis to a form
    in which the zeroth coefficient is one for one cusp and zero for the other ones.
    """
    CC = eisforms[0].base_ring
    G = eisforms[0].G
    dim_E = len(eisforms)
    cusps = G.cusps()
    M_A = MatrixSpace(CC,len(cusps),dim_E)
    M_b = MatrixSpace(CC,len(cusps),1)

    A = M_A([[eisforms[i].get_cusp_expansion(c)[0] for i in range(dim_E)] for c in cusps])

    scaling_constants = {}
    for j in range(dim_E):
        b = M_b([1 if i == j else 0 for i in range(len(cusps))])
        c = A.solve_right(b)
        scaling_constants[j] = c.list()

    new_eisforms = []
    for j in range(dim_E):
        new_eisform = 0
        for i in range(dim_E):
            new_eisform = eisforms[i]*scaling_constants[j][i] + new_eisform
        new_eisforms.append(new_eisform)
    
    if return_scaling_constants == False:
        return new_eisforms
    else:
        return new_eisforms, scaling_constants

def get_monomial_products(list_of_vars, d):
    """
    Form all possible monomials of specified vars up to degree d.
    Example:
    form_monomial_products([x,y,z],3)
    ->
    [x^3, x^2*y, x^2*z, x*y^2, x*y*z, x*z^2, y^3, y^2*z, y*z^2, z^3]
    """
    return list(map(lambda x: prod(x), combinations_with_replacement(list_of_vars,d)))

def get_algebraic_relations(eisforms, d):
    """
    Search for algebraic relations between the first non-trivial coefficients of Eisenstein series.
    Note that we assume that the Eisenstein series are normalized according to the function 'echelon_basis_to_eisenstein_basis'.
    """
    first_non_trivial_coeffs = [eisform.get_cusp_expansion(Cusp(1,0))[1] for eisform in eisforms]
    monomial_products = get_monomial_products(first_non_trivial_coeffs,d)
    return lindep(monomial_products,check_if_result_is_invalid=True)
