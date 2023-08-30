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

    M_A, M_b, M_c = MatrixSpace(CC,dim_S,dim_M), MatrixSpace(CC,dim_S,1), MatrixSpace(CC,dim_M,1)
    digit_prec = bits_to_digits(CC.precision())
    A = M_A([petersson_products[i] for i in range(dim_S)])
    trivial_c_vecs = [] #In case we have zero columns in A, we can already use these vectors as Eisenstein series
    zero_cols = []
    for col_id in range(A.ncols()):
        if is_column_effectively_zero(A,col_id,digit_prec): #If a column is zero then this means that the corresponding modform is already an Eisenstein series
            trivial_c_vecs.append(M_c([1 if i == col_id else 0 for i in range(dim_M)]))
            zero_cols.append(col_id)
    #Drop all zero columns
    A = A.delete_columns(zero_cols)
    dim_E_reduced = dim_E-len(trivial_c_vecs) #Number of Eisenstein series that are left to be computed
    #Solve the remaining linear system
    A = set_effectively_zero_entries_to_zero(A,digit_prec)
    K = A.right_kernel_matrix()
    c_vecs = [v for v in K]

    scaling_constants = {} #Constants by which the modform basis gets scaled to produce eisenstein series
    label = 0
    for c_vec in trivial_c_vecs:
        scaling_constants[label] = [c_vec[i,0] for i in range(dim_M)]
        label += 1
    for c_vec in c_vecs:
        k = 0
        scaling_constants[label] = []
        for i in range(dim_M):
            if i in zero_cols:
                scaling_constants[label].append(0)
            else:
                scaling_constants[label].append(c_vec[k])
                k += 1
        label += 1
        
    #Sort scaling constants by valuation (that is, the position of the first non-zero entry)
    non_zero_pos = []
    for factors in scaling_constants.values():
        for (i,c) in enumerate(factors):
            if c != 0:
                non_zero_pos.append(i)
                break
    scaling_constants = dict(sorted(scaling_constants.items(), key=lambda item: non_zero_pos[item[0]]))
    
    # We are now ready to construct the eisforms from the modforms
    eisforms = []
    for factors in scaling_constants.values():
        eisform = modforms_CC[0]*factors[0]
        for i in range(1,dim_M):
            eisform += modforms_CC[i]*factors[i]
        eisforms.append(eisform)

    if return_scaling_constants == False:
        return eisforms
    else:
        return eisforms, scaling_constants

def is_column_effectively_zero(A, column_id, digit_prec, tol=0.8):
    """
    Check if a column of a matrix is effectively zero.
    """
    for i in range(A.nrows()):
        if not is_effectively_zero(A[i,column_id].abs(),int(tol*digit_prec)):
            return False
    return True

def set_effectively_zero_entries_to_zero(A, digit_prec, tol=0.8):
    for i in range(A.nrows()):
        for j in range(A.ncols()):
            if is_effectively_zero(A[i,j].abs(),digit_prec):
                A[i,j] = 0
    return A

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
