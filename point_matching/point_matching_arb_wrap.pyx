# ****************************************************************************
#       Copyright (C) 2022 David Berghaus <berghaus@th.physik.uni-bonn.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import math
import numpy as np
from cysignals.signals cimport sig_on, sig_off

from sage.libs.arb.arb cimport *
from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_mat cimport *
from sage.rings.complex_arb cimport *
from sage.matrix.matrix_complex_ball_dense cimport *
from sage.rings.real_arb import RealBallField
from sage.rings.complex_arb import ComplexBallField
from sage.matrix.matrix_space import MatrixSpace
from sage.arith.misc import prime_factors
from sage.rings.real_double import RealDoubleElement

from psage.modform.maass.automorphic_forms_alg import get_M_for_holom

from arblib_helpers.acb_approx cimport *
from pullback.my_pullback cimport my_pullback_pts_arb_wrap, apply_moebius_transformation_arb_wrap
from classes.acb_mat_class cimport Acb_Mat, Acb_Mat_Win
from classes.block_factored_mat_class cimport Block_Factored_Mat, Block_Factored_Element, J_class, W_class
from classes.acb_dft_class cimport Acb_DFT
from classes.plu_class cimport PLU_Mat
from iterative_solvers.gmres_arb_wrap import gmres_mgs_arb_wrap
from iterative_solvers.iterative_refinement_arb_wrap import iterative_refinement_arb_wrap

cdef _get_J_block_matrix_arb_wrap(acb_mat_t J,int Ms,int Mf,int weight,int Q,coordinates,int bit_prec):
    cdef int coord_len = len(coordinates)
    cdef RR = RealBallField(bit_prec)
    cdef CC = ComplexBallField(bit_prec)
    cdef ComplexBall two_pi_i = CC(0,2*get_pi_ball(bit_prec))
    cdef ComplexBall z_horo, czd, fact, tmp, exp_one
    cdef ComplexBall weight_fact = CC(1,0)
    cdef RealBall c,d
    cdef int j, n
    cdef RealBall one_over_2Q = RR(1)/(2*Q)

    for j in range(coord_len):
        (z_horo,_,_,c,d) = coordinates[j]
        x_horo = z_horo.real()
        if weight != 0:
            czd = c*z_horo+d
            weight_fact = (czd.abs()/czd)**weight
        fact = weight_fact*one_over_2Q
        exp_one = (-two_pi_i*x_horo).exp()
        tmp = fact*((-two_pi_i*Ms*x_horo).exp()) #We could use exp_one here for performance
        acb_set(acb_mat_entry(J, 0, j), tmp.value)
        for n in range(Ms+1,Mf+1):
            acb_approx_mul(acb_mat_entry(J, n-Ms, j), acb_mat_entry(J, n-Ms-1, j), exp_one.value, bit_prec)

cpdef get_pi_ball(int bit_prec): #Since RR(pi) does not compile...
    RR = RealBallField(bit_prec)
    cdef RealBall pi_ball = RR(0)
    cdef arb_t arb_pi
    arb_init(arb_pi)
    arb_const_pi(arb_pi, bit_prec)
    arb_set(pi_ball.value,arb_pi)
    arb_clear(arb_pi)
    return pi_ball

cpdef digits_to_bits(int digits): #Approximates how many bits are required to achieve this amount of digits
    return math.ceil(digits*math.log(10,2)) #log(x,2) is log_2

cpdef bits_to_digits(int bit_prec):
    return math.ceil(bit_prec*math.log10(2))

cdef _get_W_block_matrix_arb_wrap(acb_mat_t W,int Ms,int Mf,int weight,coordinates,int bit_prec,trunc_W=True,mix_precisions=False):
    """
    Compute W matrix and store it (inplace) into W.
    We also provide additional parameters:
    trunc_W: Tries to (experimentally) truncate columns of W up to a sufficient order
    mix_precisions: Tries to (experimentally) compute entries of W to the least required precision.
                    We have currently disabled this by default because it does not seem to improve the performance significantly
    """
    cdef int weight_half = weight//2
    cdef int coord_len = len(coordinates)
    cdef RR = RealBallField(bit_prec)
    cdef CC = ComplexBallField(bit_prec)
    cdef ComplexBall two_pi_i = CC(0,2*get_pi_ball(bit_prec))
    cdef int j, l
    cdef ComplexBall z_horo, z_fund, tmp, exp_one
    cdef RealBall a, b, c, d
    cdef RealBall y_fund_fact = RR(1)

    #Hopefully we can remove this section once arb adds fast Horner schemes...
    cdef CF_low_prec = ComplexBallField(53)
    cdef double log10_pow_minus_D = -1.1*bits_to_digits(bit_prec)*math.log(10)

    for j in range(coord_len):
        (z_horo,a,b,c,d) = coordinates[j]
        z_fund = apply_moebius_transformation_arb_wrap(z_horo,a,b,c,d)
        if weight != 0:
            y_fund_fact = (z_fund.imag())**weight_half
        exp_one = (two_pi_i*z_fund).exp() #We could in principle re-use this for the normalized column to get some more performance...
        tmp = y_fund_fact*(exp_one**Ms)
        acb_set(acb_mat_entry(W, j, 0), tmp.value)
        if trunc_W == True: #Hopefully we can remove this section once arb adds fast Horner schemes...
            #This feature is naive and experimental. This should only be a temporary solution until Horner uses auto-truncation.
            suggested_trunc_order = int(log10_pow_minus_D/(RealDoubleElement(CF_low_prec(exp_one).abs().log())))
            trunc_order = min(suggested_trunc_order,Mf)
        else:
            trunc_order = Mf
        if mix_precisions == True:
            exp_one_log_2 = RealDoubleElement(CF_low_prec(exp_one).abs().log(2)) #Maybe there is a better way to extract the binary exponent...
            for l in range(Ms+1,trunc_order+1):
                computation_bit_prec = max(int(bit_prec+l*exp_one_log_2+weight_half*math.log2(l)+10),64)
                acb_approx_mul(acb_mat_entry(W, j, l-Ms), acb_mat_entry(W, j, l-Ms-1), exp_one.value, computation_bit_prec)
        else:
            for l in range(Ms+1,trunc_order+1):
                acb_approx_mul(acb_mat_entry(W, j, l-Ms), acb_mat_entry(W, j, l-Ms-1), exp_one.value, bit_prec)

cdef _compute_V_block_matrix_arb_wrap(acb_mat_t V_view,acb_mat_t J,int Ms,int Mf,int weight,coordinates,int bit_prec): #computes a V-block-matrix and stores it in V
    coord_len = len(coordinates)
    cdef acb_mat_t W
    acb_mat_init(W, coord_len, Mf-Ms+1)
    _get_W_block_matrix_arb_wrap(W,Ms,Mf,weight,coordinates,bit_prec)
    sig_on()
    acb_mat_approx_mul(V_view,J,W,bit_prec)
    sig_off()
    acb_mat_clear(W) #It would be more efficient to re-use these allocations...

cdef _compute_V_block_matrix_normalized_column_arb_wrap(acb_mat_t b_view,J_class J,int l_normalized,int weight,coordinates,int bit_prec): #Computes column of V corresponding to l_normalized
    coord_len = len(coordinates)
    cdef acb_mat_t W
    acb_mat_init(W, coord_len, 1)
    _get_W_block_matrix_arb_wrap(W,l_normalized,l_normalized,weight,coordinates,bit_prec)
    J.act_on_vec(b_view,W,bit_prec)
    acb_mat_clear(W) #It would be more efficient to re-use these allocations...

cdef Acb_Mat get_diagonal_terms(int Ms,int Mf,int weight,Y,int bit_prec):
    cdef int M = Mf-Ms+1
    cdef int weight_half, i
    cdef RR = RealBallField(bit_prec)
    cdef CC = ComplexBallField(bit_prec)
    cdef RealBall Y_pow_weight_half, two_pi, tmp
    cdef ComplexBall exp_one
    weight_half = weight//2
    Y_pow_weight_half = Y**weight_half
    two_pi = 2*get_pi_ball(bit_prec)
    cdef Acb_Mat diag = Acb_Mat(M, 1) #We could use real entries here...
    tmp = Y_pow_weight_half*((-two_pi*Ms*Y).exp())
    acb_set_arb(acb_mat_entry(diag.value,0,0),tmp.value)
    exp_one = CC((-two_pi*Y).exp(),0)
    for i in range(Ms+1,Mf+1):
        acb_approx_mul(acb_mat_entry(diag.value,i-Ms,0), acb_mat_entry(diag.value,i-Ms-1,0), exp_one.value, bit_prec)

    return diag

cdef Acb_Mat get_diagonal_inv_terms(int Ms,int Mf,int weight,Y,int bit_prec):
    cdef int M = Mf-Ms+1
    cdef int weight_half, i
    cdef RR = RealBallField(bit_prec)
    cdef CC = ComplexBallField(bit_prec)
    cdef RealBall Y_pow_minus_weight_half, two_pi, tmp
    cdef ComplexBall exp_one
    minus_weight_half = -weight//2
    Y_pow_minus_weight_half = Y**minus_weight_half
    two_pi = 2*get_pi_ball(bit_prec)
    cdef Acb_Mat diag_inv = Acb_Mat(M, 1) #We could use real entries here...
    tmp = Y_pow_minus_weight_half*((two_pi*Ms*Y).exp())
    acb_set_arb(acb_mat_entry(diag_inv.value,0,0),tmp.value)
    exp_one = CC((two_pi*Y).exp(),0)
    for i in range(Ms+1,Mf+1):
        acb_approx_mul(acb_mat_entry(diag_inv.value,i-Ms,0), acb_mat_entry(diag_inv.value,i-Ms-1,0), exp_one.value, bit_prec)

    return diag_inv

def get_inv_roots_of_unity(int N, int bit_prec):
    """
    Compute N inverse roots of unity. The computations are performed recursively.
    """
    cdef int i
    cdef Acb_Mat inv_roots_of_unity = Acb_Mat(1,N) #Maybe replace this with vector later
    cdef RR = RealBallField(bit_prec)
    cdef CC = ComplexBallField(bit_prec)
    cdef ComplexBall exp_one = CC(0,-2*get_pi_ball(bit_prec)/N).exp()

    acb_one(acb_mat_entry(inv_roots_of_unity.value,0,0))
    for i in range(1,N):
        acb_approx_mul(acb_mat_entry(inv_roots_of_unity.value,0,i), acb_mat_entry(inv_roots_of_unity.value,0,i-1), exp_one.value, bit_prec)
    
    return inv_roots_of_unity

cdef _subtract_diagonal_terms(acb_mat_t V_view,int Ms,int Mf,int weight,Y,int bit_prec):
    """
    Transforms V to V_tilde by subtracting the diagonal elements
    """
    cdef int M = Mf-Ms+1
    cdef int weight_half, i
    cdef RR = RealBallField(bit_prec)
    cdef CC = ComplexBallField(bit_prec)
    cdef RealBall Y_pow_weight_half, two_pi, tmp
    weight_half = weight//2
    Y_pow_weight_half = Y**weight_half
    two_pi = 2*get_pi_ball(bit_prec)
    for i in range(Ms,Mf+1):
        tmp = Y_pow_weight_half*((-two_pi*i*Y).exp())
        acb_sub_arb(acb_mat_entry(V_view,i-Ms,i-Ms),acb_mat_entry(V_view,i-Ms,i-Ms),tmp.value,bit_prec)

cpdef get_V_tilde_matrix_cuspform_arb_wrap(S,int M,int Q,Y,int bit_prec):
    cdef int weight = S.weight()
    cdef int Ms = 1
    cdef int Mf = M
    pb = my_pullback_pts_arb_wrap(S,1-Q,Q,Y,bit_prec)
    G = S.group()
    cdef int nc = G.ncusps()
    cdef Acb_Mat V = Acb_Mat(nc*M,nc*M)
    cdef int cii,cjj
    cdef Acb_Mat_Win V_view
    cdef Acb_Mat J
    for cii in range(nc):
        for cjj in range(nc):
            coordinates = pb[cii][cjj]['coordinates']
            coord_len = len(coordinates)
            V_view = V.get_window(cii*M,cjj*M,(cii+1)*M,(cjj+1)*M)
            if coord_len != 0:
                J = Acb_Mat(M, coord_len)
                _get_J_block_matrix_arb_wrap(J.value,Ms,Mf,weight,Q,coordinates,bit_prec) #we compute J here to re-use it later for normalized column
                _compute_V_block_matrix_arb_wrap(V_view.value,J.value,Ms,Mf,weight,coordinates,bit_prec)
            if cii == cjj:
                _subtract_diagonal_terms(V_view.value,Ms,Mf,weight,Y,bit_prec)
    return V

cpdef _get_l_normalized(cjj,normalization,is_cuspform):
    """
    Return indices i of c_i where c_i=1 (i.e. return index of normalized coefficients).
    'starting_index' refers to the first index of the expansion which is
    1 for cuspforms and 0 for modforms series.
    For cuspforms we currently only support a normalization basis in reduced echelon-form (as this always seems to work)
    while for modforms we also support normalizations like [1,1,0].
    """
    cdef int i
    if is_cuspform:
        starting_index = 1
    else:
        starting_index = 0
    l_normalized = []
    for i in range(len(normalization[cjj])):
        tmp = normalization[cjj][i]
        if tmp == 0 or tmp == 1:
            if tmp == 1:
                l_normalized.append(i+starting_index)
        else:
            raise ArithmeticError("Normalization list currently only supports entries zero or one.")
    if l_normalized == []:
        raise ArithmeticError("Could not determine l_normalized...")
    return l_normalized

def _get_echelon_normalization_from_label(label, multiplicity):
    """
    Get echelon normalization for the principal cusp, based on label and multiplicity.
    Examples:
    label==0, multiplicity==3 -> [1,0,0]
    label==1, multiplicity==3 -> [0,1,0]
    """
    if label < 0 or label+1 > multiplicity:
        raise ArithmeticError("Invalid label!")
    res = [0]*multiplicity
    res[label] = 1
    return res

def get_higher_genera_normalizations(S, is_cuspform):
    """
    Get normalizations for all forms in S, where S is an automorphic space of a subgroup with g>0.
    We determine the normalizations by verifying that the forms exist at a low precision.
    """
    G = S.group()
    if is_cuspform:
        dim = int(G.dimension_cusp_forms(S.weight()))
    else:
        dim = int(G.dimension_modular_forms(S.weight()))
    if dim == 0:
        raise NameError("The space of forms is of dimension zero for this weight!")
    normalizations = []

    #First determine the form with the highest valuation
    max_valuation = dim
    while True:
        imposed_zeros = []
        normalization = {}
        normalization[0] = _get_echelon_normalization_from_label(max_valuation-1,max_valuation)
        for i in range(1,G.ncusps()):
            normalization[i] = []
        if is_normalization_valid(S,normalization,imposed_zeros,is_cuspform):
            normalizations.append((normalization,imposed_zeros))
            break
        max_valuation += 1
    
    #Now determine the other normalizations
    if max_valuation == dim: #We can just use a Victor Miller basis
        normalizations = []
        for label in range(dim):
            normalizations.append(_get_victor_miller_normalization(S,is_cuspform,label=label))
        return normalizations

    zero_coeffs = [max_valuation] #Indices of the coefficients which we normalize to be zero
    d = max_valuation
    j = d-2 #Because we already know that the form with valuation d-1 exists
    while len(normalizations) < dim:
        normalization = {}
        for i in range(1,G.ncusps()):
            normalization[i] = []

        #First try reduced echelon normalization
        #We do this because it might happen that the undetermined coefficients are zero,
        #in which case our normalization works, but the conditioning of the linear system is bad.
        normalization[0] = _get_echelon_normalization_from_label(j,max_valuation)
        if is_normalization_valid(S,normalization,[],is_cuspform):
            zero_coeffs.append(j+1)
            normalizations.append((normalization,[]))
            j -= 1
            continue

        #Now try normalization with imposed zeros
        normalization[0] = _get_echelon_normalization_from_label(j,d)
        imposed_zeros = []
        for zero_coeff in zero_coeffs:
            imposed_zero = zero_coeff-d-1
            if imposed_zero > 0:
                imposed_zeros.append(imposed_zero)
        imposed_zeros.sort()
        if is_normalization_valid(S,normalization,imposed_zeros,is_cuspform):
            zero_coeffs.append(j+1)
            normalizations.append((normalization,imposed_zeros))
            j -= 1
        else:
            d = j
            j -= 1
            if j < 0:
                raise ArithmeticError("Could not determine normalizations...")
    if len(normalizations) != dim:
        raise ArithmeticError("Could not determine normalizations!")
    normalizations.reverse() #We want the normalizations in increasing order of valuation
    return normalizations

def is_normalization_valid(S, normalization, imposed_zeros, is_cuspform):
    """
    To Do:
    We could make this function faster by searching if forms exist using the double precision functionalities.
    """
    digit_prec = 30
    max_iter = 30
    bit_prec = digits_to_bits(digit_prec)
    M_0 = get_M_0(S,digit_prec,is_cuspform=is_cuspform)
    RBF = RealBallField(bit_prec)
    Y_fact = 0.9
    eps = RBF(1e-10)
    Y = get_horo_height_arb_wrap(S,RBF,M_0,is_cuspform=is_cuspform)*RBF(Y_fact)
    
    try:
        if is_cuspform:
            c_1 = get_coefficients_cuspform_ir_arb_wrap(S,digit_prec,max_iter=max_iter,normalization=normalization,imposed_zeros=imposed_zeros)._get_mcbd(bit_prec)
            c_2 = get_coefficients_cuspform_ir_arb_wrap(S,digit_prec,max_iter=max_iter,Y=Y,M_0=M_0,normalization=normalization,imposed_zeros=imposed_zeros)._get_mcbd(bit_prec)
        else:
            c_1 = get_coefficients_modform_ir_arb_wrap(S,digit_prec,max_iter=max_iter,normalization=normalization,imposed_zeros=imposed_zeros)._get_mcbd(bit_prec)
            c_2 = get_coefficients_modform_ir_arb_wrap(S,digit_prec,max_iter=max_iter,Y=Y,M_0=M_0,normalization=normalization,imposed_zeros=imposed_zeros)._get_mcbd(bit_prec)
    except ArithmeticError: #IR didn't converge -> invalid normalization
        return False
    coeff_threshold = c_1.dimensions()[0]//4 #We only check the first 1/4 of the coefficients
    avg_diff = sum([abs(x) for x in (c_1[:coeff_threshold]-c_2[:coeff_threshold]).list()])/coeff_threshold
    if avg_diff > eps:
        return False
    return True

def _get_victor_miller_normalization(S, is_cuspform, label=0):
    """
    Returns normalization for each cusp. For modforms the first expansion coefficient is c_0, for cuspforms the first expansion coefficient is c_1.
    """
    G = S.group()
    if is_cuspform:
        multiplicity = int(G.dimension_cusp_forms(S.weight()))
    else:
        multiplicity = int(G.dimension_modular_forms(S.weight()))
    normalization = dict()
    if multiplicity == 0:
        raise NameError("The space of forms is of dimension zero for this weight!")
    else:
        normalization[0] = _get_echelon_normalization_from_label(label, multiplicity)
    for i in range(1,G.ncusps()):
        normalization[i] = []
    imposed_zeros = []
    return normalization, imposed_zeros

def _get_normalization_hauptmodul(S):
    """
    Returns normalization for each cusp. For modforms the first expansion coefficient is c_0.
    We treat the 1/q term separately when computing the 'b-vector'.
    """
    G = S.group()
    normalization = dict()
    normalization[0] = [0] #This is our principal cusp
    for i in range(1,G.ncusps()):
        normalization[i] = []
    return normalization

cpdef get_V_tilde_matrix_b_arb_wrap(S,int M,int Q,Y,int bit_prec,bint is_cuspform):
    """
    Returns V_tilde,b of V_tilde*x=b where b corresponds to (minus) the column at c_l_normalized
    """
    cdef int weight = S.weight()
    G = S.group()
    pb = my_pullback_pts_arb_wrap(S,1-Q,Q,Y,bit_prec)
    cdef int nc = G.ncusps()
    normalization, _ = _get_victor_miller_normalization(S,is_cuspform)

    cdef Acb_Mat V = Acb_Mat(nc*M,nc*M)
    cdef Acb_Mat b = Acb_Mat(nc*M,1)
    cdef int cii,cjj
    cdef Acb_Mat_Win V_view, b_view
    cdef J_class J
    use_FFT = False #We currently don't support the matrix multiplication J*W through FFTs

    for cii in range(nc):
        for cjj in range(nc):
            coordinates = pb[cii][cjj]['coordinates']
            coord_len = len(coordinates)
            V_view = V.get_window(cii*M,cjj*M,(cii+1)*M,(cjj+1)*M)
            Msjj = len(normalization[cjj])
            if is_cuspform:
                Msjj += 1
            Mfjj = Msjj+M-1
            Msii = len(normalization[cii])
            if is_cuspform:
                Msii += 1
            Mfii = Msii+M-1
            if coord_len != 0:
                J = J_class(use_FFT)
                J._construct(M,Msii,Mfii,weight,Q,coordinates,bit_prec,use_FFT)
                _compute_V_block_matrix_arb_wrap(V_view.value,J.J.value,Msjj,Mfjj,weight,coordinates,bit_prec)
                if cjj == 0:
                    b_view = b.get_window(cii*M,0,(cii+1)*M,1)
                    l_normalized = _get_l_normalized(cjj,normalization,is_cuspform)
                    if len(l_normalized) != 1:
                        raise ArithmeticError("We have not implemented this scenario for cuspforms yet.")
                    l_normalized = l_normalized[0]
                    _compute_V_block_matrix_normalized_column_arb_wrap(b_view.value,J,l_normalized,weight,coordinates,bit_prec)
            if cii == cjj:
                _subtract_diagonal_terms(V_view.value,Msjj,Mfjj,weight,Y,bit_prec)
    sig_on()
    acb_mat_neg(b.value, b.value)
    sig_off()
    return V, b

cpdef get_V_tilde_matrix_factored_b_arb_wrap(S,int M,int Q,Y,int bit_prec,bint use_FFT,bint use_splitting,bint is_cuspform,normalizations,labels=None):
    """
    Returns V_tilde,b of V_tilde*x=b where b corresponds to (minus) the column at c_l_normalized.
    V_tilde is not explicitly computed but instead consists of block-matrices of the form J*W
    """
    cdef int weight = S.weight()
    G = S.group()
    pb = my_pullback_pts_arb_wrap(S,1-Q,Q,Y,bit_prec)
    cdef int nc = G.ncusps()
    if labels == None:
        if is_cuspform:
            multiplicity = int(G.dimension_cusp_forms(S.weight()))
        else:
            multiplicity = int(G.dimension_modular_forms(S.weight()))
        labels = range(multiplicity)

    cdef Block_Factored_Mat block_factored_mat = Block_Factored_Mat(nc)
    if use_FFT == True:
        DFT_precomp = Acb_DFT(2*Q, bit_prec)
        inv_roots_of_unity = get_inv_roots_of_unity(2*Q, bit_prec)
    else:
        DFT_precomp = None
        inv_roots_of_unity = None
    V_factored = block_factored_mat.A
    diag_factored = block_factored_mat.diag
    diag_inv_factored = block_factored_mat.diag_inv
    b_vecs = [Acb_Mat(nc*M,1) for _ in range(len(labels))]
    cdef int cii, cjj
    cdef Acb_Mat_Win b_view
    cdef J_class J
    cdef W_class W
    cdef Block_Factored_Element block_factored_element

    for cii in range(nc):
        for cjj in range(nc):
            coordinates = pb[cii][cjj]['coordinates']
            coord_len = len(coordinates)
            Msjj = len(normalizations[0][cjj])
            if is_cuspform:
                Msjj += 1
            Mfjj = Msjj+M-1
            Msii = len(normalizations[0][cii])
            if is_cuspform:
                Msii += 1
            Mfii = Msii+M-1
            if coord_len != 0:
                j_values = pb[cii][cjj]['j_values']
                V_factored[cii][cjj] = Block_Factored_Element(J_class(use_FFT), W_class(use_splitting))
                block_factored_element = V_factored[cii][cjj]
                J = block_factored_element.J
                J._construct(M,Msii,Mfii,weight,Q,coordinates,bit_prec,use_FFT,j_values=j_values,DFT_precomp=DFT_precomp,inv_roots_of_unity=inv_roots_of_unity)
                W = block_factored_element.W
                W._construct(M,Msjj,Mfjj,weight,coordinates,bit_prec,use_splitting)
                if cjj == 0:
                    for ni in range(len(normalizations)):
                        normalization = normalizations[ni]
                        b_view = b_vecs[ni].get_window(cii*M,0,(cii+1)*M,1)
                        l_normalized = _get_l_normalized(cjj,normalization,is_cuspform)
                        if len(l_normalized) != 1: #This would be a normalization like [1,1,0]
                            raise ArithmeticError("We have not implemented this scenario for cuspforms yet.")
                        _compute_V_block_matrix_normalized_column_arb_wrap(b_view.value,J,l_normalized[0],weight,coordinates,bit_prec)
                        sig_on()
                        acb_mat_neg(b_view.value, b_view.value)
                        sig_off()
            if cii == cjj:
                diag_factored[cii] = get_diagonal_terms(Msjj,Mfjj,weight,Y,bit_prec)
                diag_inv_factored[cii] = get_diagonal_inv_terms(Msjj,Mfjj,weight,Y,bit_prec)

    block_factored_mat._init_parameters_for_dp_construction(S,M,Q,Y,normalizations[0],pb,is_cuspform)

    return block_factored_mat, b_vecs

cpdef get_V_tilde_matrix_factored_b_haupt_arb_wrap(S,int M,int Q,Y,int bit_prec,bint use_FFT,bint use_splitting):
    """
    Returns V_tilde,b of V_tilde*x=b where b corresponds to (minus) the column at c_l_normalized.
    V_tilde is not explicitly computed but instead consists of block-matrices of the form J*W.
    This function is devoted to computing the hauptmodul for genus zero surfaces with normalization c_{-1} = 1, c_0 = 0 at the principal cusp
    """
    cdef int weight = S.weight()
    if weight != 0:
        raise NameError("This function only works for weight zero!")
    G = S.group()
    if G.genus() != 0:
        raise NameError("This function only works for genus zero surfaces!")
    pb = my_pullback_pts_arb_wrap(S,1-Q,Q,Y,bit_prec)
    normalization = _get_normalization_hauptmodul(S)
    cdef int nc = G.ncusps()

    cdef Block_Factored_Mat block_factored_mat = Block_Factored_Mat(nc)
    if use_FFT == True:
        DFT_precomp = Acb_DFT(2*Q, bit_prec)
        inv_roots_of_unity = get_inv_roots_of_unity(2*Q, bit_prec)
    else:
        DFT_precomp = None
        inv_roots_of_unity = None
    V_factored = block_factored_mat.A
    diag_factored = block_factored_mat.diag
    diag_inv_factored = block_factored_mat.diag_inv
    cdef Acb_Mat b = Acb_Mat(nc*M,1)
    cdef int cii, cjj
    cdef Acb_Mat_Win b_view
    cdef J_class J
    cdef W_class W
    cdef Block_Factored_Element block_factored_element

    for cii in range(nc):
        for cjj in range(nc):
            coordinates = pb[cii][cjj]['coordinates']
            coord_len = len(coordinates)
            Msjj = len(normalization[cjj])
            Msii = len(normalization[cii])
            Mfjj = Msjj+M-1
            Mfii = Msii+M-1
            if coord_len != 0:
                j_values = pb[cii][cjj]['j_values']
                V_factored[cii][cjj] = Block_Factored_Element(J_class(use_FFT), W_class(use_splitting))
                block_factored_element = V_factored[cii][cjj]
                J = block_factored_element.J
                J._construct(M,Msii,Mfii,weight,Q,coordinates,bit_prec,use_FFT,j_values=j_values,DFT_precomp=DFT_precomp,inv_roots_of_unity=inv_roots_of_unity)
                W = block_factored_element.W
                W._construct(M,Msjj,Mfjj,weight,coordinates,bit_prec,use_splitting)
                if cjj == 0:
                    b_view = b.get_window(cii*M,0,(cii+1)*M,1)
                    _compute_V_block_matrix_normalized_column_arb_wrap(b_view.value,J,-1,weight,coordinates,bit_prec)
            if cii == cjj:
                diag_factored[cii] = get_diagonal_terms(Msjj,Mfjj,weight,Y,bit_prec)
                diag_inv_factored[cii] = get_diagonal_inv_terms(Msjj,Mfjj,weight,Y,bit_prec)

    sig_on()
    acb_mat_neg(b.value, b.value)
    sig_off()

    block_factored_mat._init_parameters_for_dp_construction(S,M,Q,Y,normalization,pb,False)

    return block_factored_mat, b

def get_M_0(S, digit_prec, is_cuspform=True):
    """
    Get truncation order of modular form. M_0 = M(Y_0) where Y_0 is the the point in the fundamental region with the smallest height.
    """
    Y_0 = S.group().minimal_height()
    if is_cuspform == True: #For cuspforms, the coefficients grow like O(n^(weight/2))
        weight = S.weight()
    else:
        if S.weight() != 0:
            weight = 2*S.weight() #For modforms, the coefficients grow like O(n^(weight-1))
        else:
            weight = 12 #TO DO: ADD PROPER ASYMPTOTIC FORMULA FOR HAUPTMODUL
    M_0 = math.ceil(get_M_for_holom(Y_0,weight,digit_prec))
    return M_0

cpdef get_horo_height_arb_wrap(S, RBF, M_0, prec_loss=None, is_cuspform=True):
    """
    Get a choice of the horocycle height. 'prec_loss' denotes the expected loss of precision of the last (M_0th) coefficient.
    """
    if prec_loss == None: #We only want the first few coefficients to have full precision, so we just choose some Y < Y_0
        return RBF(S.group().minimal_height()*0.8)
    else:
        Y_1 = S.group().minimal_height()*0.8
        Y_2 = prec_loss*math.log(10)/(2*math.pi*M_0)
        Y = min(Y_1,Y_2)
        return RBF(Y)

cpdef get_Q(Y, weight, digit_prec, is_cuspform=True):
    """
    Get amount of sampling points 'Q', based on choice of horocycle height 'Y'.
    """
    if is_cuspform == False:
        if weight != 0:
            weight = 2*weight #For modforms, the coefficients grow like O(n^(weight-1))
        else:
            weight = 12 #TO DO: ADD PROPER ASYMPTOTIC FORMULA FOR HAUPTMODUL
    Q_min = math.ceil(get_M_for_holom(Y,weight,digit_prec))+1

    #Now choose Q in a way such that it has small prime factors to make life easier for the FFT
    max_prime_factors = []
    for i in range(10): #Look at the next ten Q-values to check which has the smallest prime factors
        max_prime_factors.append(max(prime_factors(Q_min+i)))
    if min(max_prime_factors) <= 13: #Found a decent choice for Q
        Q = Q_min+max_prime_factors.index(min(max_prime_factors))
    else: #Keep searching until we are satisfied
        Q = Q_min+10
        while max(prime_factors(Q)) > 13:
            Q += 1

    return Q

cpdef get_coefficients_cuspform_arb_wrap(S,int digit_prec,Y=0,int M_0=0,int Q=0,prec_loss=None):
    """
    Computes expansion coefficients with direct methods (i.e. explicitly constructs V_tilde and performs a LU-decomposition)
    """
    bit_prec = digits_to_bits(digit_prec)
    RBF = RealBallField(bit_prec)
    if M_0 == 0:
        M_0 = get_M_0(S,digit_prec)
    if float(Y) == 0: #This comparison does not seem to be defined for arb-types...
        Y = get_horo_height_arb_wrap(S,RBF,M_0,prec_loss=prec_loss)
    if Q == 0:
        Q = get_Q(Y,S.weight(),digit_prec,is_cuspform=True)
    print("Y = ", Y)
    print("M_0 = ", M_0)
    print("Q = ", Q)
    print("ncusps = ", S.group().ncusps())
    cdef Acb_Mat V,b
    V,b = get_V_tilde_matrix_b_arb_wrap(S,M_0,Q,Y,bit_prec,True)
    sig_on()
    acb_mat_approx_solve(b.value,V.value,b.value,bit_prec)
    sig_off()
    return b.get_window(0,0,M_0,1)

cpdef get_coefficients_gmres_cuspform_arb_wrap(S,int digit_prec,Y=0,int M_0=0,int Q=0,label=0,prec_loss=None,use_FFT=True,use_splitting=False):
    """ 
    Computes expansion coefficients using GMRES, preconditioned with low_prec LU-decomposition
    """
    use_scipy_lu = False #For GMRES we cannot use the scipy LU because we need to cast the LU matrix to working precision
    bit_prec = digits_to_bits(digit_prec)
    RBF = RealBallField(bit_prec)
    CBF = ComplexBallField(bit_prec)
    if M_0 == 0:
        M_0 = get_M_0(S,digit_prec)
    if float(Y) == 0: #This comparison does not seem to be defined for arb-types...
        Y = get_horo_height_arb_wrap(S,RBF,M_0,prec_loss=prec_loss)
    if Q == 0:
        Q = get_Q(Y,S.weight(),digit_prec,is_cuspform=True)
    print("Y = ", Y)
    print("M_0 = ", M_0)
    print("Q = ", Q)
    print("ncusps = ", S.group().ncusps())
    cdef Block_Factored_Mat V
    cdef Acb_Mat b, res
    cdef PLU_Mat plu

    V, b_vecs = get_V_tilde_matrix_factored_b_arb_wrap(S,M_0,Q,Y,bit_prec,use_FFT,use_splitting,True,labels=[label])
    b = b_vecs[0]
    tol = RBF(10.0)**(-digit_prec+1)

    V_dp = V.construct_sc_np()
    plu = PLU_Mat(V_dp,53,use_scipy_lu)

    x_gmres_arb_wrap = gmres_mgs_arb_wrap(V, b, bit_prec, tol, PLU=plu)

    res = x_gmres_arb_wrap[0]
    V.diag_inv_scale_vec(res, res, bit_prec)

    return res.get_window(0,0,M_0,1)

cpdef get_coefficients_gmres_non_scaled_cuspform_arb_wrap(S,int digit_prec,Y=0,int M_0=0,int Q=0,prec_loss=None):
    """
    Constructs V_tilde explicitly (uses the non-scaled version) and afterwards uses GMRES to solve the resolving linear system of equations.
    WARNING:
    ONLY USE THE FUNCTION FOR TESTING BECAUSE IT IS VERY INEFFICIENT!
    """
    bit_prec = digits_to_bits(digit_prec)
    RBF = RealBallField(bit_prec)
    if M_0 == 0:
        M_0 = get_M_0(S,digit_prec)
    if float(Y) == 0: #This comparison does not seem to be defined for arb-types...
        Y = get_horo_height_arb_wrap(S,RBF,M_0,prec_loss=prec_loss)
    if Q == 0:
        Q = get_Q(Y,S.weight(),digit_prec,is_cuspform=True)
    print("Y = ", Y)
    print("M_0 = ", M_0)
    print("Q = ", Q)
    print("ncusps = ", S.group().ncusps())
    cdef Acb_Mat V,b
    V,b = get_V_tilde_matrix_b_arb_wrap(S,M_0,Q,Y,bit_prec,True)
    tol = RBF(10.0)**(-digit_prec+1)
    x_gmres_arb_wrap = gmres_mgs_arb_wrap(V, b, bit_prec, tol, maxiter=10**4)
    res = x_gmres_arb_wrap[0]
    return res.get_window(0,0,M_0,1)

cpdef get_coefficients_cuspform_ir_arb_wrap(S,int digit_prec,Y=0,int M_0=0,int Q=0,return_M=False,label=0,prec_loss=None,use_FFT=True,use_splitting=True,use_scipy_lu=True,max_iter=None,normalization=None,imposed_zeros=None):
    """ 
    Computes expansion coefficients of cuspform using classical iterative refinement
    """
    bit_prec = digits_to_bits(digit_prec)
    RBF = RealBallField(bit_prec)
    CBF = ComplexBallField(bit_prec)
    if M_0 == 0:
        M_0 = get_M_0(S,digit_prec)
    if float(Y) == 0: #This comparison does not seem to be defined for arb-types...
        Y = get_horo_height_arb_wrap(S,RBF,M_0,prec_loss=prec_loss)
    if Q == 0:
        Q = get_Q(Y,S.weight(),digit_prec,is_cuspform=True)
    print("Y = ", Y)
    print("M_0 = ", M_0)
    print("Q = ", Q)
    print("ncusps = ", S.group().ncusps())
    cdef Block_Factored_Mat V
    cdef Acb_Mat b, res
    cdef PLU_Mat plu

    if normalization is None and imposed_zeros is None:
        normalization, imposed_zeros = _get_victor_miller_normalization(S,True,label=label)

    V, b_vecs = get_V_tilde_matrix_factored_b_arb_wrap(S,M_0,Q,Y,bit_prec,use_FFT,use_splitting,True,[normalization],labels=[label])
    b = b_vecs[0]
    tol = RBF(10.0)**(-digit_prec+1)

    V_dp = V.construct_sc_np()
    if imposed_zeros != None and len(imposed_zeros) > 0: #Delete the rows and columns corresponding to the imposed zeros
        V_dp = np.delete(V_dp, imposed_zeros, 0)
        V_dp = np.delete(V_dp, imposed_zeros, 1)
    plu = PLU_Mat(V_dp,53,use_scipy_lu)

    res = iterative_refinement_arb_wrap(V, b, bit_prec, tol, plu, maxiter=max_iter, imposed_zeros=imposed_zeros)
    for imposed_zero in imposed_zeros:
        acb_zero(acb_mat_entry(res.value,imposed_zero,0))

    V.diag_inv_scale_vec(res, res, bit_prec)

    if return_M == False:
        return res
    else:
        return res, M_0

cpdef get_coefficients_modform_ir_arb_wrap(S,int digit_prec,Y=0,int M_0=0,int Q=0,return_M=False,label=0,prec_loss=None,use_FFT=True,use_splitting=True,use_scipy_lu=True,max_iter=None,normalization=None,imposed_zeros=None):
    """ 
    Computes Fourier-expansion coefficients of modforms using classical iterative refinement
    """
    bit_prec = digits_to_bits(digit_prec)
    RBF = RealBallField(bit_prec)
    CBF = ComplexBallField(bit_prec)
    if M_0 == 0:
        M_0 = get_M_0(S,digit_prec,is_cuspform=False)
    if float(Y) == 0: #This comparison does not seem to be defined for arb-types...
        Y = get_horo_height_arb_wrap(S,RBF,M_0,is_cuspform=False,prec_loss=prec_loss)
    if Q == 0:
        Q = get_Q(Y,S.weight(),digit_prec,is_cuspform=False)
    print("Y = ", Y)
    print("M_0 = ", M_0)
    print("Q = ", Q)
    print("ncusps = ", S.group().ncusps())
    cdef Block_Factored_Mat V
    cdef Acb_Mat b, res
    cdef PLU_Mat plu

    if normalization is None and imposed_zeros is None:
        normalization, imposed_zeros = _get_victor_miller_normalization(S,False,label=label)

    V, b_vecs = get_V_tilde_matrix_factored_b_arb_wrap(S,M_0,Q,Y,bit_prec,use_FFT,use_splitting,False,[normalization],labels=[label])
    b = b_vecs[0]
    tol = RBF(10.0)**(-digit_prec+1)

    V_dp = V.construct_sc_np()
    if imposed_zeros != None and len(imposed_zeros) > 0: #Delete the rows and columns corresponding to the imposed zeros
        V_dp = np.delete(V_dp, imposed_zeros, 0)
        V_dp = np.delete(V_dp, imposed_zeros, 1)
    plu = PLU_Mat(V_dp,53,use_scipy_lu)

    res = iterative_refinement_arb_wrap(V, b, bit_prec, tol, plu, maxiter=max_iter, imposed_zeros=imposed_zeros)
    for imposed_zero in imposed_zeros:
        acb_zero(acb_mat_entry(res.value,imposed_zero,0))

    V.diag_inv_scale_vec(res, res, bit_prec)

    if return_M == False:
        return res
    else:
        return res, M_0

cpdef get_coefficients_haupt_ir_arb_wrap(S,int digit_prec,Y=0,int M_0=0,int Q=0,only_principal_expansion=True,return_M=False,prec_loss=None,use_FFT=True,use_splitting=True,use_scipy_lu=True):
    """ 
    Computes expansion coefficients of hauptmodul using classical iterative refinement
    """
    bit_prec = digits_to_bits(digit_prec)
    RBF = RealBallField(bit_prec)
    CBF = ComplexBallField(bit_prec)
    if M_0 == 0:
        M_0 = get_M_0(S,digit_prec,is_cuspform=False)
    if float(Y) == 0: #This comparison does not seem to be defined for arb-types...
        Y = get_horo_height_arb_wrap(S,RBF,M_0,is_cuspform=False,prec_loss=prec_loss)
    if Q == 0:
        Q = get_Q(Y,S.weight(),digit_prec,is_cuspform=False)
    print("Y = ", Y)
    print("M_0 = ", M_0)
    print("Q = ", Q)
    print("ncusps = ", S.group().ncusps())
    cdef Block_Factored_Mat V
    cdef Acb_Mat b, res
    cdef PLU_Mat plu

    V, b = get_V_tilde_matrix_factored_b_haupt_arb_wrap(S,M_0,Q,Y,bit_prec,use_FFT,use_splitting)
    tol = RBF(10.0)**(-digit_prec+1)

    V_dp = V.construct_sc_np()
    plu = PLU_Mat(V_dp,53,use_scipy_lu)

    res = iterative_refinement_arb_wrap(V, b, bit_prec, tol, plu)

    V.diag_inv_scale_vec(res, res, bit_prec)

    if only_principal_expansion == True:
        if return_M == False:
            return res.get_window(0,0,M_0,1)
        else:
            return res.get_window(0,0,M_0,1), M_0
    else:
        if return_M == False:
            return res
        else:
            return res, M_0

cpdef get_coefficients_cuspform_ir_restarting_arb_wrap(S,digit_precs,return_M=False,label=0,use_FFT=True,use_splitting=True,use_scipy_lu=True):
    """ 
    Computes expansion coefficients of cuspform using classical iterative refinement.
    This function uses the precisions specified in 'digit_precs' to gradually increase the size of the system of linear equations
    and restart the iterative refinement process.

    Remarks:
    This function currently does not seem to be significantly faster because the construction of the W matrices takes so long.
    For this reason, we put this function on to hold until arb releases a faster Horner scheme.
    """
    G = S.group()
    cdef Block_Factored_Mat V
    cdef Acb_Mat b, res, acb_mat_wrap
    cdef PLU_Mat plu
    cdef ComplexBall acb_wrap
    M_0_values = [get_M_0(S,digit_prec) for digit_prec in digit_precs]
    c_vec_for_next_iter = None

    for (i,digit_prec) in enumerate(digit_precs):
        bit_prec = digits_to_bits(digit_prec)
        RBF = RealBallField(bit_prec)
        CBF = ComplexBallField(bit_prec)
        M_0 = M_0_values[i]
        Y = get_horo_height_arb_wrap(S,RBF,M_0)
        Q = get_Q(Y,S.weight(),digit_prec,is_cuspform=True)
        print("Y = ", Y)
        print("M_0 = ", M_0)
        print("Q = ", Q)
        print("ncusps = ", S.group().ncusps())

        V, b_vecs = get_V_tilde_matrix_factored_b_arb_wrap(S,M_0,Q,Y,bit_prec,use_FFT,use_splitting,True,labels=[label])
        b = b_vecs[0]
        tol = RBF(10.0)**(-digit_prec+1)

        V_dp = V.construct_sc_np()
        plu = PLU_Mat(V_dp,53,use_scipy_lu)

        if i == 0:
            starting_prec = 0
        else:
            starting_prec = digit_precs[i-1]
        res = iterative_refinement_arb_wrap(V, b, bit_prec, tol, plu, x0=c_vec_for_next_iter, starting_prec=starting_prec)
        
        if i < len(digit_precs)-1: #Inflate c_vec_for_next_iter to a higher precision and to a higher M by filling unknown coeffs with zeros
            res_mcbd = res._get_mcbd(bit_prec)
            next_bit_prec = digits_to_bits(digit_precs[i+1])
            next_M_0 = M_0_values[i+1]
            c_vec_for_next_iter = Acb_Mat(G.ncusps()*next_M_0,1)
            acb_mat_wrap = c_vec_for_next_iter
            for cii in range(G.ncusps()):
                for i in range(M_0):
                    acb_wrap = res_mcbd[i+cii*M_0][0]
                    acb_approx_set(acb_mat_entry(acb_mat_wrap.value,i+cii*next_M_0,0),acb_wrap.value)
        else:
            V.diag_inv_scale_vec(res, res, bit_prec) #Only do this in the last iteration

    if return_M == False:
        return res
    else:
        return res, M_0

cpdef get_cuspform_basis_ir_arb_wrap(S,int digit_prec,Y=0,int M_0=0,int Q=0,return_M_and_labels=False,labels=None,prec_loss=None,use_FFT=True,use_splitting=True,use_scipy_lu=True):
    """
    Compute a basis of cuspforms of AutomorphicFormSpace 'S' to 'digit_prec' digits precision.
    """
    bit_prec = digits_to_bits(digit_prec)
    RBF = RealBallField(bit_prec)
    CBF = ComplexBallField(bit_prec)
    if M_0 == 0:
        M_0 = get_M_0(S,digit_prec)
    if float(Y) == 0: #This comparison does not seem to be defined for arb-types...
        Y = get_horo_height_arb_wrap(S,RBF,M_0,prec_loss=prec_loss)
    if Q == 0:
        Q = get_Q(Y,S.weight(),digit_prec,is_cuspform=True)
    print("Y = ", Y)
    print("M_0 = ", M_0)
    print("Q = ", Q)
    print("ncusps = ", S.group().ncusps())
    cdef Block_Factored_Mat V
    cdef PLU_Mat plu

    V, b_vecs = get_V_tilde_matrix_factored_b_arb_wrap(S,M_0,Q,Y,bit_prec,use_FFT,use_splitting,True,labels=labels)
    tol = RBF(10.0)**(-digit_prec+1)

    V_dp = V.construct_sc_np()
    plu = PLU_Mat(V_dp,53,use_scipy_lu)

    res_vec = []
    for i in range(len(b_vecs)):
        res = iterative_refinement_arb_wrap(V, b_vecs[i], bit_prec, tol, plu)
        V.diag_inv_scale_vec(res, res, bit_prec)
        res_vec.append(res)
    
    if return_M_and_labels == False:
        return res_vec
    else:
        if labels == None:
            multiplicity = S.group().dimension_cusp_forms(S.weight())
            labels = range(multiplicity)
        return res_vec, M_0, labels

cpdef get_modform_basis_ir_arb_wrap(S,int digit_prec,Y=0,int M_0=0,int Q=0,return_M_and_labels=False,labels=None,prec_loss=None,use_FFT=True,use_splitting=True,use_scipy_lu=True):
    """
    Compute a basis of modular forms of AutomorphicFormSpace 'S' to 'digit_prec' digits precision.
    """
    bit_prec = digits_to_bits(digit_prec)
    RBF = RealBallField(bit_prec)
    CBF = ComplexBallField(bit_prec)
    if M_0 == 0:
        M_0 = get_M_0(S,digit_prec,is_cuspform=False)
    if float(Y) == 0: #This comparison does not seem to be defined for arb-types...
        Y = get_horo_height_arb_wrap(S,RBF,M_0,is_cuspform=False,prec_loss=prec_loss)
    if Q == 0:
        Q = get_Q(Y,S.weight(),digit_prec,is_cuspform=False)
    print("Y = ", Y)
    print("M_0 = ", M_0)
    print("Q = ", Q)
    print("ncusps = ", S.group().ncusps())
    cdef Block_Factored_Mat V
    cdef PLU_Mat plu

    V, b_vecs = get_V_tilde_matrix_factored_b_arb_wrap(S,M_0,Q,Y,bit_prec,use_FFT,use_splitting,False,labels=labels)
    tol = RBF(10.0)**(-digit_prec+1)

    V_dp = V.construct_sc_np()
    plu = PLU_Mat(V_dp,53,use_scipy_lu)

    res_vec = []
    for i in range(len(b_vecs)):
        res = iterative_refinement_arb_wrap(V, b_vecs[i], bit_prec, tol, plu)
        V.diag_inv_scale_vec(res, res, bit_prec)
        res_vec.append(res)

    if return_M_and_labels == False:
        return res_vec
    else:
        if labels == None:
            multiplicity = S.group().dimension_modular_forms(S.weight())
            labels = range(multiplicity)
        return res_vec, M_0, labels