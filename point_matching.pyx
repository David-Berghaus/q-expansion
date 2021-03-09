import math
import numpy as np
cimport numpy as np

cdef extern from "complex.h":
    cdef double cimag(double complex)
    cdef double creal(double complex)
    cdef double cabs(double complex)
    cdef double complex cexp(double complex)

from psage.modform.maass.automorphic_forms_alg import get_M_for_holom

from my_pullback cimport my_pullback_pts_dp

cdef _get_J_block_matrix_dp(int M,int weight,int Q,coordinates):
    cdef int coord_len = len(coordinates)
    J = np.empty(shape=(M,coord_len),dtype=np.complex_,order="F") #Fortran order allows more efficient access here
    cdef np.complex128_t [:, :] J_view = J
    cdef double complex two_pi_i = 2*math.pi*1j
    cdef double complex z_horo, czd, weight_fact, fact
    cdef int j, n, c, d
    cdef double one_over_2Q = 1.0/(2*Q)
    for j in range(coord_len):
        (z_horo,_,_,c,d,_) = coordinates[j]
        x_horo = creal(z_horo)
        czd = c*z_horo+d
        weight_fact = (cabs(czd)/czd)**weight
        fact = weight_fact*one_over_2Q
        for n in range(1,M+1):
            J_view[n-1,j] = fact*cexp(-two_pi_i*n*x_horo)
    return J

cdef _get_W_block_matrix_dp(int Ms,int Mf,int weight,coordinates):
    cdef int weight_half = weight//2
    cdef int coord_len = len(coordinates)
    W = np.empty(shape=(coord_len,Mf-Ms+1),dtype=np.complex_)
    cdef np.complex128_t [:, :] W_view = W
    cdef double complex two_pi_i = 2*math.pi*1j
    cdef int j, l
    cdef double complex z_fund
    cdef double y_fund_fact
    for j in range(coord_len):
        (_,_,_,_,_,z_fund) = coordinates[j]
        y_fund_fact = cimag(z_fund)**weight_half
        for l in range(Ms,Mf+1):
            W_view[j,l-Ms] = y_fund_fact*cexp(two_pi_i*l*z_fund)
    return W

cdef _get_W_block_matrix_normalized_column_dp(int l_normalized,int weight,coordinates): #column corresponding to c_l_normalized = 1
    cdef int weight_half = weight//2
    cdef int coord_len = len(coordinates)
    W_col = np.empty(coord_len,dtype=np.complex_)
    cdef np.complex128_t [:] W_col_view = W_col
    cdef double complex two_pi_i = 2*math.pi*1j
    cdef int j
    cdef double complex z_fund
    cdef double y_fund_fact
    for j in range(coord_len):
        (_,_,_,_,_,z_fund) = coordinates[j]
        y_fund_fact = cimag(z_fund)**weight_half
        W_col_view[j] = y_fund_fact*cexp(two_pi_i*l_normalized*z_fund)
    return W_col

cdef _compute_V_block_matrix_dp(V_view,J,int cii,int cjj,int Ms,int Mf,int weight,double Y,coordinates): #computes a V-block-matrix and stores it in V
    W = _get_W_block_matrix_dp(Ms,Mf,weight,coordinates)
    np.matmul(J,W,out=V_view)

cdef _compute_V_tilde_block_matrix_dp(V_view,J,int cii,int cjj,int Ms,int Mf,int weight,double Y,coordinates): #computes a V_tilde-block-matrix and stores it in V
    _compute_V_block_matrix_dp(V_view,J,cii,cjj,Ms,Mf,weight,Y,coordinates)
    cdef int M = Mf-Ms+1
    cdef int weight_half, i
    cdef double Y_pow_weight_half, two_pi
    if cii == cjj:
        weight_half = weight//2
        Y_pow_weight_half = Y**weight_half
        two_pi = 2*math.pi
        for i in range(Ms,M+1):
            V_view[i-1,i-Ms] -= Y_pow_weight_half*cexp(-two_pi*i*Y)

cdef _compute_V_tilde_block_matrix_normalized_column_dp(V_column_view,J,int cii,int cjj,int l_normalized,int weight,double Y,coordinates): #computes a normalized column of V_tilde-block-matrix and stores it in V
    W_col = _get_W_block_matrix_normalized_column_dp(l_normalized,weight,coordinates)
    np.matmul(J,W_col,out=V_column_view)
    if cii == cjj:
        weight_half = weight//2
        Y_pow_weight_half = Y**weight_half
        two_pi = 2*math.pi
        V_column_view[l_normalized-1] -= Y_pow_weight_half*cexp(-two_pi*l_normalized*Y)

cpdef get_V_tilde_matrix_dp(S,int M,double Y,int weight):
    cdef int Ms = 1
    cdef int Mf = M
    cdef int Q = M+8
    pb = my_pullback_pts_dp(S,1-Q,Q,Y)
    cdef int nc = S.group().ncusps()
    G = S.group()
    V = np.zeros(shape=(nc*M,nc*M),dtype=np.complex_)
    cdef int cii,cjj
    for cii in range(nc):
        for cjj in range(nc):
            coordinates = pb[cii][cjj]
            if len(coordinates) != 0:
                J = _get_J_block_matrix_dp(M,weight,Q,coordinates) #we compute J here to re-use it later for normalized column
                V_view = V[cii*M:(cii+1)*M,cjj*M:(cjj+1)*M] #using a memory view is certainly not efficient here
                _compute_V_tilde_block_matrix_dp(V_view,J,cii,cjj,Ms,Mf,weight,Y,coordinates)
    return V

cpdef get_V_tilde_matrix_b_dp(S,int M,double Y,int weight,int multiplicity): #Returns V_tilde,b of V_tilde*x=b where b corresponds to (minus) the column at c_l_normalized
    cdef int Q = M+8
    pb = my_pullback_pts_dp(S,1-Q,Q,Y)
    G = S.group()
    cdef int nc = G.ncusps()
    normalization = dict()
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
    V = np.zeros(shape=(nc*M,nc*M),dtype=np.complex_)
    b = np.zeros(nc*M,dtype=np.complex_)
    cdef int cii,cjj,Ms,Mf,l_normalized
    for cii in range(nc):
        for cjj in range(nc):
            coordinates = pb[cii][cjj]
            if len(coordinates) != 0:
                J = _get_J_block_matrix_dp(M,weight,Q,coordinates) #we compute J here to re-use it later for normalized column
                V_view = V[cii*M:(cii+1)*M,cjj*M:(cjj+1)*M] #using a memory view is certainly not efficient here
                Ms = len(normalization[cjj])+1
                Mf = Ms+M-1
                _compute_V_tilde_block_matrix_dp(V_view,J,cii,cjj,Ms,Mf,weight,Y,coordinates)
                if cjj == 0:
                    b_view = b[cii*M:(cii+1)*M]
                    for i in range(len(normalization[cjj])):
                        if normalization[cjj][i] != 0:
                            l_normalized = i+1
                    _compute_V_tilde_block_matrix_normalized_column_dp(b_view,J,cii,cjj,l_normalized,weight,Y,coordinates)
    np.negative(b,out=b)
    return V,b

cpdef get_V_matrix_dp(S,int M,double Y,int weight):
    cdef int Ms = 1
    cdef int Mf = M
    cdef int Q = M+8
    pb = my_pullback_pts_dp(S,1-Q,Q,Y)
    cdef int nc = S.group().ncusps()
    G = S.group()
    V = np.zeros(shape=(nc*M,nc*M),dtype=np.complex_)
    cdef int cii,cjj
    for cii in range(nc):
        for cjj in range(nc):
            coordinates = pb[cii][cjj]
            if len(coordinates) != 0:
                J = _get_J_block_matrix_dp(M,weight,Q,coordinates) #we compute J here to re-use it later for normalized column
                V_view = V[cii*M:(cii+1)*M,cjj*M:(cjj+1)*M] #using a memory view is certainly not efficient here
                _compute_V_block_matrix_dp(V_view,J,cii,cjj,Ms,Mf,weight,Y,coordinates)
    return V

def get_coefficients(S,int weight,int multiplicity,double Y=0,int M=0,prec=14):
    if Y == 0:
        Y = S.group().minimal_height()*0.8
    if M == 0:
        M = math.ceil(1.2*get_M_for_holom(Y,weight,prec))
    print("Y = ", Y)
    print("M = ", M)
    # V = get_V_tilde_matrix_dp(S,M,Y,weight)
    # Nrows, Ncols = V.shape 
    # A = V[:(Nrows-1),1:]
    # b = -V[:(Nrows-1),0]
    V,b = get_V_tilde_matrix_b_dp(S,M,Y,weight,multiplicity)
    return np.linalg.solve(V,b)

# cdef get_V_tilde_element(int n,int l,int cii,int cjj,int Q,int weight,double Y,coordinates):
#     cdef double complex two_pi_i_n = 2*math.pi*1j*n
#     cdef double complex two_pi_i_l = 2*math.pi*1j*l
#     cdef int weight_half = weight//2
#     cdef double complex V_el = 0.
#     cdef double complex weight_fact, z_horo, z_pb
#     cdef int a,b,c,d
#     for (z_horo,a,b,c,d,z_pb) in coordinates:
#         weight_fact = (cabs(c*z_horo+d)/(c*z_horo+d))**weight
#         V_el += weight_fact*cimag(z_pb)**weight_half*cexp(two_pi_i_l*z_pb)*cexp(-two_pi_i_n*creal(z_horo))
#     cdef double one_over_2Q = 1.0/(2*Q)
#     V_el *= one_over_2Q
#     if cii == cjj:
#         if n == l:
#             V_el -= Y**weight_half*cexp(-2*math.pi*n*Y)
#     return V_el

# cpdef get_V_tilde_matrix_by_element_dp(S,int M,double Y,int weight):
#     cdef int Ms = 1
#     cdef int Mf = M
#     cdef int Q = M+8
#     pb = my_pullback_pts_dp(S,1-Q,Q,Y)
#     cdef int nc = S.group().ncusps()
#     G = S.group()
#     V = np.zeros(shape=(nc*M,nc*M),dtype=np.complex_)
#     cdef int cii,cjj,n,l
#     for cii in range(nc):
#         for cjj in range(nc):
#             V_view = V[cii*M:(cii+1)*M,cjj*M:(cjj+1)*M]
#             coordinates = pb[cii][cjj]
#             if len(coordinates) != 0:
#                 for n in range(1,M+1):
#                     for l in range(1,M+1):
#                         V_view[n-1,l-1] = get_V_tilde_element(n,l,cii,cjj,Q,weight,Y,coordinates)
#     return V

# cpdef get_coefficients_test(S,int M,double Y,int weight):
#     V = get_V_tilde_matrix_by_element_dp(S,M,Y,weight)
#     Nrows, Ncols = V.shape 
#     A = V[:(Nrows-1),1:]
#     b = -V[:(Nrows-1),0]
#     return np.linalg.solve(A,b)