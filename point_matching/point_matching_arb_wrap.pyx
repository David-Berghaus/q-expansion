import math
from cysignals.signals cimport sig_on, sig_off

from sage.libs.arb.arb cimport *
from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_mat cimport *
from sage.rings.complex_arb cimport *
from sage.matrix.matrix_complex_ball_dense cimport *
from sage.rings.real_arb import RealBallField
from sage.rings.complex_arb import ComplexBallField
from sage.matrix.matrix_space import MatrixSpace

from psage.modform.maass.automorphic_forms_alg import get_M_for_holom

from arblib_helpers.acb_approx cimport *
from pullback.my_pullback cimport my_pullback_pts_arb_wrap, apply_moebius_transformation_arb_wrap
from classes.acb_mat_class cimport Acb_Mat, Acb_Mat_Win
from classes.block_factored_mat_class cimport Block_Factored_Mat, Block_Factored_Element
from classes.plu_class cimport PLU_Mat
from iterative_solvers.gmres_arb_wrap import gmres_mgs_arb_wrap

cdef _get_J_block_matrix_arb_wrap(acb_mat_t J,int Ms,int Mf,int weight,int Q,coordinates,int bit_prec):
    cdef int coord_len = len(coordinates)
    cdef RR = RealBallField(bit_prec)
    cdef CC = ComplexBallField(bit_prec)
    cdef ComplexBall two_pi_i = CC(0,2*get_pi_ball(bit_prec))
    cdef ComplexBall z_horo, czd, weight_fact, fact, tmp, exp_one
    cdef RealBall c,d
    cdef int j, n
    cdef RealBall one_over_2Q = RR(1)/(2*Q)

    for j in range(coord_len):
        (z_horo,_,_,c,d,_) = coordinates[j]
        x_horo = z_horo.real()
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
    return math.ceil(digits*math.log(10)/math.log(2)-4)

cdef _get_W_block_matrix_arb_wrap(acb_mat_t W,int Ms,int Mf,int weight,coordinates,int bit_prec):
    cdef int weight_half = weight//2
    cdef int coord_len = len(coordinates)
    cdef RR = RealBallField(bit_prec)
    cdef CC = ComplexBallField(bit_prec)
    cdef ComplexBall two_pi_i = CC(0,2*get_pi_ball(bit_prec))
    cdef int j, l
    cdef ComplexBall z_horo, z_fund, tmp, exp_one
    cdef RealBall a, b, c, d, y_fund_fact
    for j in range(coord_len):
        (z_horo,a,b,c,d,_) = coordinates[j]
        z_fund = apply_moebius_transformation_arb_wrap(z_horo,a,b,c,d)
        y_fund_fact = (z_fund.imag())**weight_half
        exp_one = (two_pi_i*z_fund).exp()
        tmp = y_fund_fact*((two_pi_i*Ms*z_fund).exp()) #We could use exp_one here for performance
        acb_set(acb_mat_entry(W, j, 0), tmp.value)
        for l in range(Ms+1,Mf+1):
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

cdef _compute_V_block_matrix_normalized_column_arb_wrap(acb_mat_t b_view,acb_mat_t J,int cii,int cjj,int l_normalized,int weight,Y,coordinates,int bit_prec): #Computes column of V corresponding to l_normalized
    coord_len = len(coordinates)
    cdef acb_mat_t W
    acb_mat_init(W, coord_len, 1)
    _get_W_block_matrix_arb_wrap(W,l_normalized,l_normalized,weight,coordinates,bit_prec)
    sig_on()
    acb_mat_approx_mul(b_view,J,W,bit_prec)
    sig_off()
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

cpdef get_V_tilde_matrix_arb_wrap(S,int M,Y,int bit_prec):
    cdef int weight = S.weight()
    cdef int Ms = 1
    cdef int Mf = M
    cdef int Q = M+8
    pb = my_pullback_pts_arb_wrap(S,1-Q,Q,Y,bit_prec)
    G = S.group()
    cdef int nc = G.ncusps()
    cdef Acb_Mat V = Acb_Mat(nc*M,nc*M)
    cdef int cii,cjj
    cdef Acb_Mat_Win V_view
    cdef Acb_Mat J
    for cii in range(nc):
        for cjj in range(nc):
            coordinates = pb[cii][cjj]
            coord_len = len(coordinates)
            V_view = V.get_window(cii*M,cjj*M,(cii+1)*M,(cjj+1)*M)
            if coord_len != 0:
                J = Acb_Mat(M, coord_len)
                _get_J_block_matrix_arb_wrap(J.value,Ms,Mf,weight,Q,coordinates,bit_prec) #we compute J here to re-use it later for normalized column
                _compute_V_block_matrix_arb_wrap(V_view.value,J.value,Ms,Mf,weight,coordinates,bit_prec)
            if cii == cjj:
                _subtract_diagonal_terms(V_view.value,Ms,Mf,weight,Y,bit_prec)
    return V

cdef _get_l_normalized(cjj,normalization):
    cdef int i
    cdef int l_normalized = 0
    for i in range(len(normalization[cjj])):
        if normalization[cjj][i] != 0:
            l_normalized = i+1 #we set c_l_normalized = 1
    if l_normalized == 0:
        raise NameError("Could not determine l_normalized...")
    return l_normalized

cdef _get_normalization(S):
    """
    Returns normalization for each cusp.
    A feature to be implemented in the future is to test if the normalization
    works correctly (with a double-precision computation).
    """
    G = S.group()
    cdef int multiplicity = G.dimension_cusp_forms(S.weight()) #!!!Might not always work (consider using dimension_cuspforms from MySubgroup)
    normalization = dict()
    if multiplicity == 0:
        raise NameError("The space of cuspforms is of dimension zero for this weight!")
    if multiplicity == 1:
        normalization[0] = [1]
    if multiplicity == 2:
        print("Careful, this normalization might not work for all groups!")
        print("")
        normalization[0] = [1,0] #this corresponds to c1=1, c2=0 for the first cusp
    if multiplicity > 2:
        raise NameError("This case has not been implemented yet")
    for i in range(1,G.ncusps()):
        normalization[i] = []
    return normalization

cpdef get_V_tilde_matrix_b_arb_wrap(S,int M,Y,int bit_prec):
    """
    Returns V_tilde,b of V_tilde*x=b where b corresponds to (minus) the column at c_l_normalized
    """
    cdef int weight = S.weight()
    G = S.group()
    cdef int Q = M+8
    pb = my_pullback_pts_arb_wrap(S,1-Q,Q,Y,bit_prec)
    cdef int nc = G.ncusps()
    normalization = _get_normalization(S)

    cdef Acb_Mat V = Acb_Mat(nc*M,nc*M)
    cdef Acb_Mat b = Acb_Mat(nc*M,1)
    cdef int cii,cjj
    cdef Acb_Mat_Win V_view, b_view
    cdef Acb_Mat J

    for cii in range(nc):
        for cjj in range(nc):
            coordinates = pb[cii][cjj]
            coord_len = len(coordinates)
            V_view = V.get_window(cii*M,cjj*M,(cii+1)*M,(cjj+1)*M)
            Msjj = len(normalization[cjj])+1
            Mfjj = Msjj+M-1
            Msii = len(normalization[cii])+1
            Mfii = Msii+M-1
            if coord_len != 0:
                J = Acb_Mat(M,coord_len)
                _get_J_block_matrix_arb_wrap(J.value,Msii,Mfii,weight,Q,coordinates,bit_prec) #we compute J here to re-use it later for normalized column
                _compute_V_block_matrix_arb_wrap(V_view.value,J.value,Msjj,Mfjj,weight,coordinates,bit_prec)
                if cjj == 0:
                    b_view = b.get_window(cii*M,0,(cii+1)*M,1)
                    l_normalized = _get_l_normalized(cjj,normalization)
                    _compute_V_block_matrix_normalized_column_arb_wrap(b_view.value,J.value,cii,cjj,l_normalized,weight,Y,coordinates,bit_prec)
            if cii == cjj:
                _subtract_diagonal_terms(V_view.value,Msjj,Mfjj,weight,Y,bit_prec)
    sig_on()
    acb_mat_neg(b.value, b.value)
    sig_off()
    return V, b

cpdef get_V_tilde_matrix_factored_b_arb_wrap(S,int M,Y,int bit_prec):
    """
    Returns V_tilde,b of V_tilde*x=b where b corresponds to (minus) the column at c_l_normalized.
    V_tilde is not explicitly computed but instead consists of block-matrices of the form J*W
    """
    cdef int weight = S.weight()
    G = S.group()
    cdef int Q = M+8
    pb = my_pullback_pts_arb_wrap(S,1-Q,Q,Y,bit_prec)
    cdef int nc = G.ncusps()
    normalization = _get_normalization(S)

    cdef Block_Factored_Mat block_factored_mat = Block_Factored_Mat(nc)
    V_factored = block_factored_mat.A
    diag_factored = block_factored_mat.diag
    diag_inv_factored = block_factored_mat.diag_inv
    cdef Acb_Mat b = Acb_Mat(nc*M,1)
    cdef int cii, cjj
    cdef Acb_Mat_Win b_view
    cdef Acb_Mat J, W
    cdef Block_Factored_Element block_factored_element

    for cii in range(nc):
        for cjj in range(nc):
            coordinates = pb[cii][cjj]
            coord_len = len(coordinates)
            Msjj = len(normalization[cjj])+1
            Mfjj = Msjj+M-1
            Msii = len(normalization[cii])+1
            Mfii = Msii+M-1
            if coord_len != 0:
                V_factored[cii][cjj] = Block_Factored_Element(Acb_Mat(M,coord_len), Acb_Mat(coord_len, Mfjj-Msjj+1))
                block_factored_element = V_factored[cii][cjj]
                J = block_factored_element.J
                W = block_factored_element.W
                _get_J_block_matrix_arb_wrap(J.value,Msii,Mfii,weight,Q,coordinates,bit_prec)
                _get_W_block_matrix_arb_wrap(W.value,Msjj,Mfjj,weight,coordinates,bit_prec)
                if cjj == 0:
                    b_view = b.get_window(cii*M,0,(cii+1)*M,1)
                    l_normalized = _get_l_normalized(cjj,normalization)
                    _compute_V_block_matrix_normalized_column_arb_wrap(b_view.value,J.value,cii,cjj,l_normalized,weight,Y,coordinates,bit_prec)
            if cii == cjj:
                diag_factored[cii] = get_diagonal_terms(Msjj,Mfjj,weight,Y,bit_prec)
                diag_inv_factored[cii] = get_diagonal_inv_terms(Msjj,Mfjj,weight,Y,bit_prec)

    sig_on()
    acb_mat_neg(b.value, b.value)
    sig_off()
    return block_factored_mat, b

cpdef get_coefficients_arb_wrap(S,int digit_prec,Y=0,int M=0):
    """
    Computes expansion coefficients with direct methods (i.e. explicitly constructs V_tilde and performs a LU-decomposition)
    """
    bit_prec = digits_to_bits(digit_prec)
    RR = RealBallField(bit_prec)
    if float(Y) == 0: #This comparison does not seem to be defined for arb-types...
        Y = RR(S.group().minimal_height()*0.8)
    if M == 0:
        weight = S.weight()
        M = math.ceil(get_M_for_holom(Y,weight,digit_prec))
    print("Y = ", Y)
    print("M = ", M)
    print("dimen = ", S.group().ncusps()*M)
    cdef Acb_Mat V,b
    V,b = get_V_tilde_matrix_b_arb_wrap(S,M,Y,bit_prec)
    sig_on()
    acb_mat_approx_solve(b.value,V.value,b.value,bit_prec)
    sig_off()
    return b.get_window(0,0,M,1)

cpdef get_coefficients_gmres_arb_wrap(S,int digit_prec,Y=0,int M=0):
    """ 
    Computes expansion coefficients with using GMRES, preconditioned with low_prec LU-decomposition
    """
    bit_prec = digits_to_bits(digit_prec)
    RBF = RealBallField(bit_prec)
    CBF = ComplexBallField(bit_prec)
    if float(Y) == 0: #This comparison does not seem to be defined for arb-types...
        Y = RBF(S.group().minimal_height()*0.8)
    if M == 0:
        weight = S.weight()
        M = math.ceil(get_M_for_holom(Y,weight,digit_prec))
    print("Y = ", Y)
    print("M = ", M)
    print("dimen = ", S.group().ncusps()*M)
    cdef Block_Factored_Mat V
    cdef Acb_Mat b, res
    cdef PLU_Mat plu

    V, b = get_V_tilde_matrix_factored_b_arb_wrap(S,M,Y,bit_prec)
    tol = RBF(10.0)**(-digit_prec+1)

    V_dp = V.construct_sc_np()
    plu = PLU_Mat(V_dp,prec=53)

    x_gmres_arb_wrap = gmres_mgs_arb_wrap(V, b, bit_prec, tol, PLU=plu)

    res = x_gmres_arb_wrap[0]
    V.diag_inv_scale_vec(res, res, bit_prec)

    return res.get_window(0,0,M,1)