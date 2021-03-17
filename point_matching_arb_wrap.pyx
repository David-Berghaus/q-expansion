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

from acb_mat_approx cimport *
from my_pullback cimport my_pullback_pts_arb_wrap, apply_moebius_transformation_arb_wrap

cdef _get_J_block_matrix_arb_wrap(acb_mat_t J,int Ms,int Mf,int weight,int Q,coordinates,int bit_prec):
    cdef int coord_len = len(coordinates)
    cdef RR = RealBallField(bit_prec)
    cdef CC = ComplexBallField(bit_prec)
    cdef ComplexBall two_pi_i = CC(0,2*get_pi_ball(bit_prec))
    cdef ComplexBall z_horo, czd, weight_fact, fact, tmp
    cdef RealBall c,d
    cdef int j, n
    cdef RealBall one_over_2Q = RR(1)/(2*Q)
    for j in range(coord_len):
        (z_horo,_,_,c,d,_) = coordinates[j]
        x_horo = z_horo.real()
        czd = c*z_horo+d
        weight_fact = (czd.abs()/czd)**weight
        fact = weight_fact*one_over_2Q
        for n in range(Ms,Mf+1):
            tmp = fact*((-two_pi_i*n*x_horo).exp())
            acb_set(acb_mat_entry(J, n-Ms, j), tmp.value)

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
    cdef ComplexBall z_horo, z_fund, tmp
    cdef RealBall a, b, c, d, y_fund_fact
    for j in range(coord_len):
        (z_horo,a,b,c,d,_) = coordinates[j]
        z_fund = apply_moebius_transformation_arb_wrap(z_horo,a,b,c,d)
        y_fund_fact = (z_fund.imag())**weight_half
        for l in range(Ms,Mf+1):
            tmp = y_fund_fact*((two_pi_i*l*z_fund).exp())
            acb_set(acb_mat_entry(W, j, l-Ms), tmp.value)

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

cdef _subtract_diagonal_terms(acb_mat_t V_view,int Ms,int Mf,int weight,Y,int bit_prec): #transforms V to V_tilde by subtracting the diagonal elements
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

cpdef get_V_matrix_arb_wrap(S,int M,Y,int bit_prec):
    cdef int weight = S.weight()
    cdef int Ms = 1
    cdef int Mf = M
    cdef int Q = M+8
    pb = my_pullback_pts_arb_wrap(S,1-Q,Q,Y,bit_prec)
    cdef int nc = S.group().ncusps()
    G = S.group()
    cdef RR = RealBallField(bit_prec)
    cdef CC = ComplexBallField(bit_prec)
    cdef Matrix_complex_ball_dense V = MatrixSpace(CC,nc*M,nc*M).zero()
    cdef int cii,cjj
    cdef acb_mat_t V_view, J
    for cii in range(nc):
        for cjj in range(nc):
            coordinates = pb[cii][cjj]
            if len(coordinates) != 0:
                coord_len = len(coordinates)
                acb_mat_init(J, M, coord_len)
                _get_J_block_matrix_arb_wrap(J,Ms,Mf,weight,Q,coordinates,bit_prec) #we compute J here to re-use it later for normalized column
                acb_mat_window_init(V_view,V.value,cii*M,cjj*M,(cii+1)*M,(cjj+1)*M)
                _compute_V_block_matrix_arb_wrap(V_view,J,Ms,Mf,weight,coordinates,bit_prec)
                acb_mat_clear(J) #It would be more efficient to re-use these allocations...
                acb_mat_window_clear(V_view)
    return V

cpdef get_V_tilde_matrix_arb_wrap(S,int M,Y,int bit_prec):
    cdef int weight = S.weight()
    cdef int Ms = 1
    cdef int Mf = M
    cdef int Q = M+8
    pb = my_pullback_pts_arb_wrap(S,1-Q,Q,Y,bit_prec)
    cdef int nc = S.group().ncusps()
    G = S.group()
    cdef RR = RealBallField(bit_prec)
    cdef CC = ComplexBallField(bit_prec)
    cdef Matrix_complex_ball_dense V = MatrixSpace(CC,nc*M,nc*M).zero()
    cdef int cii,cjj
    cdef acb_mat_t V_view, J
    for cii in range(nc):
        for cjj in range(nc):
            coordinates = pb[cii][cjj]
            coord_len = len(coordinates)
            acb_mat_window_init(V_view,V.value,cii*M,cjj*M,(cii+1)*M,(cjj+1)*M)
            if coord_len != 0:
                acb_mat_init(J, M, coord_len)
                _get_J_block_matrix_arb_wrap(J,Ms,Mf,weight,Q,coordinates,bit_prec) #we compute J here to re-use it later for normalized column
                _compute_V_block_matrix_arb_wrap(V_view,J,Ms,Mf,weight,coordinates,bit_prec)
                acb_mat_clear(J) #It would be more efficient to re-use these allocations...
            if cii == cjj:
                _subtract_diagonal_terms(V_view,Ms,Mf,weight,Y,bit_prec)
            acb_mat_window_clear(V_view)
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

cpdef get_V_tilde_matrix_b_arb_wrap(S,int M,Y,int bit_prec): #Returns V_tilde,b of V_tilde*x=b where b corresponds to (minus) the column at c_l_normalized
    cdef int weight = S.weight()
    G = S.group()
    cdef int multiplicity = G.dimension_cusp_forms(weight) #!!!Might not always work (consider using dimension_cuspforms from MySubgroup)
    cdef int Q = M+8
    pb = my_pullback_pts_arb_wrap(S,1-Q,Q,Y,bit_prec)
    cdef RR = RealBallField(bit_prec)
    cdef CC = ComplexBallField(bit_prec)
    cdef int nc = G.ncusps()
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
    for i in range(1,nc):
        normalization[i] = []
    cdef Matrix_complex_ball_dense V = MatrixSpace(CC,nc*M,nc*M).zero()
    cdef Matrix_complex_ball_dense b = MatrixSpace(CC,nc*M,1).zero()
    cdef int cii,cjj
    cdef acb_mat_t V_view, J, b_view
    for cii in range(nc):
        for cjj in range(nc):
            coordinates = pb[cii][cjj]
            coord_len = len(coordinates)
            acb_mat_window_init(V_view,V.value,cii*M,cjj*M,(cii+1)*M,(cjj+1)*M)
            Msjj = len(normalization[cjj])+1
            Mfjj = Msjj+M-1
            Msii = len(normalization[cii])+1
            Mfii = Msii+M-1
            if coord_len != 0:
                acb_mat_init(J,M,coord_len)
                _get_J_block_matrix_arb_wrap(J,Msii,Mfii,weight,Q,coordinates,bit_prec) #we compute J here to re-use it later for normalized column
                _compute_V_block_matrix_arb_wrap(V_view,J,Msjj,Mfjj,weight,coordinates,bit_prec)
                if cjj == 0:
                    acb_mat_window_init(b_view,b.value,cii*M,0,(cii+1)*M,1)
                    l_normalized = _get_l_normalized(cjj,normalization)
                    _compute_V_block_matrix_normalized_column_arb_wrap(b_view,J,cii,cjj,l_normalized,weight,Y,coordinates,bit_prec)
                    acb_mat_window_clear(b_view)
                acb_mat_clear(J)
            if cii == cjj:
                _subtract_diagonal_terms(V_view,Msjj,Mfjj,weight,Y,bit_prec)
            acb_mat_window_clear(V_view)
    sig_on()
    acb_mat_neg(b.value, b.value)
    sig_off()
    return V,b

def get_coefficients_arb_wrap(S,int digit_prec,Y=0,int M=0):
    bit_prec = digits_to_bits(digit_prec)
    RR = RealBallField(bit_prec)
    if float(Y) == 0: #This comparison does not seem to be defined for arb-types...
        Y = RR(S.group().minimal_height()*0.8)
    if M == 0:
        weight = S.weight()
        M = math.ceil(get_M_for_holom(Y,weight,digit_prec))
    print("Y = ", Y)
    print("M = ", M)
    cdef Matrix_complex_ball_dense V,b
    V,b = get_V_tilde_matrix_b_arb_wrap(S,M,Y,bit_prec)
    sig_on()
    acb_mat_approx_solve(b.value,V.value,b.value,bit_prec)
    sig_off()
    return b[:M]

# def test_low_precision_inverse(S): #Test if matrix can be inverted with much lower precision so that we can later use it as a preconditioner
#     bit_prec = digits_to_bits(500)
#     RR = RealBallField(bit_prec)
#     Y = RR(S.group().minimal_height()*0.8)
#     M = 376
#     cdef Matrix_complex_ball_dense V,b
#     V,b = get_V_tilde_matrix_b_arb_wrap(S,M,Y,12,1,bit_prec)
#     print "V computed"
#     cdef acb_mat_t inv_low, inv_high
#     acb_mat_init(inv_low, M, M)
#     acb_mat_init(inv_high, M, M)
#     acb_mat_approx_inv(inv_low,V.value,64)
#     print "low inverse computed"
#     acb_mat_approx_inv(inv_high,V.value,bit_prec)
#     print "high inverse computed"
#     acb_printd(acb_mat_entry(inv_low,0,0), 16)
#     acb_printd(acb_mat_entry(inv_high,0,0), 16)
#     print('')
#     acb_printd(acb_mat_entry(inv_low,0,M-1), 16)
#     acb_printd(acb_mat_entry(inv_high,0,M-1), 16)