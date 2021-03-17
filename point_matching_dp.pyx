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

cdef _get_J_block_matrix_dp(int Ms,int Mf,int weight,int Q,coordinates):
    cdef int coord_len = len(coordinates)
    J = np.empty(shape=(Mf-Ms+1,coord_len),dtype=np.complex_,order="F") #Fortran order allows more efficient access here
    cdef np.complex128_t [:, :] J_view = J
    cdef double complex two_pi_i = 2*math.pi*1j
    cdef double complex z_horo, czd, weight_fact, fact
    cdef double c, d
    cdef int j, n
    cdef double one_over_2Q = 1.0/(2*Q)
    for j in range(coord_len):
        (z_horo,_,_,c,d,_) = coordinates[j]
        x_horo = creal(z_horo)
        czd = c*z_horo+d
        weight_fact = (cabs(czd)/czd)**weight
        fact = weight_fact*one_over_2Q
        for n in range(Ms,Mf+1):
            J_view[n-Ms,j] = fact*cexp(-two_pi_i*n*x_horo)
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

cdef _compute_V_block_matrix_dp(V_view,J,int Ms,int Mf,int weight,coordinates): #computes a V-block-matrix and stores it in V
    W = _get_W_block_matrix_dp(Ms,Mf,weight,coordinates)
    np.matmul(J,W,out=V_view)

cdef _compute_V_block_matrix_normalized_column_dp(b_view,J,int cii,int cjj,int l_normalized,int weight,double Y,coordinates): #Computes column of V corresponding to l_normalized
    W = _get_W_block_matrix_dp(l_normalized,l_normalized,weight,coordinates)
    np.matmul(J,W,out=b_view)

cdef _subtract_diagonal_terms(V_view,int Ms,int Mf,int weight,double Y): #transforms V to V_tilde by subtracting the diagonal elements
    cdef int weight_half, i
    cdef double Y_pow_weight_half, two_pi
    weight_half = weight//2
    Y_pow_weight_half = Y**weight_half
    two_pi = 2*math.pi
    for i in range(Ms,Mf+1):
        V_view[i-Ms,i-Ms] -= Y_pow_weight_half*cexp(-two_pi*i*Y)

cpdef get_V_tilde_matrix_dp(S,int M,double Y):
    cdef int weight = S.weight()
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
            V_view = V[cii*M:(cii+1)*M,cjj*M:(cjj+1)*M] #using a memory view is certainly not efficient here
            if len(coordinates) != 0:
                J = _get_J_block_matrix_dp(Ms,Mf,weight,Q,coordinates)
                _compute_V_block_matrix_dp(V_view,J,Ms,Mf,weight,coordinates)
            if cii == cjj:
                _subtract_diagonal_terms(V_view,Ms,Mf,weight,Y)
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

cpdef get_V_tilde_matrix_b_dp(S,int M,double Y): #Returns V_tilde,b of V_tilde*x=b where b corresponds to (minus) the column at c_l_normalized
    cdef int weight = S.weight()
    G = S.group()
    cdef int multiplicity = G.dimension_cusp_forms(weight) #!!!Might not always work (consider using dimension_cuspforms from MySubgroup)
    cdef int Q = M+8
    pb = my_pullback_pts_dp(S,1-Q,Q,Y)
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
    V = np.zeros(shape=(nc*M,nc*M),dtype=np.complex_)
    b = np.zeros(shape=(nc*M,1),dtype=np.complex_)
    cdef int cii,cjj
    for cii in range(nc):
        for cjj in range(nc):
            coordinates = pb[cii][cjj]
            V_view = V[cii*M:(cii+1)*M,cjj*M:(cjj+1)*M] #using a memory view is certainly not efficient here
            Msjj = len(normalization[cjj])+1
            Mfjj = Msjj+M-1
            Msii = len(normalization[cii])+1
            Mfii = Msii+M-1
            if len(coordinates) != 0:
                J = _get_J_block_matrix_dp(Msii,Mfii,weight,Q,coordinates) #we compute J here to re-use it later for normalized column
                _compute_V_block_matrix_dp(V_view,J,Msjj,Mfjj,weight,coordinates)
                if cjj == 0:
                    b_view = b[cii*M:(cii+1)*M] #Weird python would create a (nc*M,) shape out of b_view[cii*M:(cii+1)*M,0]...
                    l_normalized = _get_l_normalized(cjj,normalization)
                    _compute_V_block_matrix_normalized_column_dp(b_view,J,cii,cjj,l_normalized,weight,Y,coordinates)
            if cii == cjj:
                _subtract_diagonal_terms(V_view,Msjj,Mfjj,weight,Y)
    np.negative(b,out=b)
    return V,b

cpdef get_V_matrix_dp(S,int M,double Y):
    cdef int weight = S.weight()
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
                J = _get_J_block_matrix_dp(Ms,Mf,weight,Q,coordinates) #we compute J here to re-use it later for normalized column
                V_view = V[cii*M:(cii+1)*M,cjj*M:(cjj+1)*M] #using a memory view is certainly not efficient here
                _compute_V_block_matrix_dp(V_view,J,Ms,Mf,weight,coordinates)
    return V

def get_coefficients_dp(S,double Y=0,int M=0,prec=14):
    if Y == 0:
        Y = S.group().minimal_height()*0.9
    if M == 0:
        weight = S.weight()
        M = math.ceil(get_M_for_holom(Y,weight,prec))
    print("Y = ", Y)
    print("M = ", M)
    V,b = get_V_tilde_matrix_b_dp(S,M,Y)
    c = np.linalg.solve(V,b)
    return c[:M] #Expansion coefficients at first cusp