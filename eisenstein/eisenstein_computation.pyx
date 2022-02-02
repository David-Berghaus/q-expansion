from sage.matrix.matrix_space import MatrixSpace

from eisenstein.haberland import compute_petersson_product_haberland, get_exp_two_pi_i_a_m_dict, clear_memoized_caches
from sage.functions.other import conjugate
from sage.symbolic.constants import pi

def compute_eisenstein_series(cuspforms, modforms, return_scaling_constants=False):
    """
    Compute the orthogonal complement of the cuspforms in the space of modular forms (which corresponds to a basis of Eisenstein series).
    We assume that cuspforms/modforms are lists of basis functions in reduced row echelon form.
    This function then returns a basis of Eisenstein series in reduced row echelon form.
    """
    cuspforms_CC = [cuspform._convert_to_CC() for cuspform in cuspforms] #We need to work over ComplexFields because CBFs have hashing issues in haberland
    modforms_CC = [modform._convert_to_CC() for modform in modforms]
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

    if return_scaling_constants == False:
        return eisforms
    else:
        return eisforms, scaling_constants