from cysignals.signals cimport sig_on, sig_str, sig_off
from cysignals.memory cimport sig_free,sig_malloc
import numpy as np
import math
cimport numpy as np
from scipy.sparse import csc_matrix, linalg as sla

cdef extern from "complex.h":
    cdef double cimag(double complex)
    cdef double creal(double complex)
    cdef double cabs(double complex)
    cdef double complex cexp(double complex)

from sage.libs.arb.arb cimport *
from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_mat cimport *
from sage.matrix.matrix_complex_ball_dense cimport *

from arblib_helpers.acb_approx cimport *
from classes.acb_mat_class cimport Acb_Mat

def get_V_tilde_matrix_sc_cuspform_np(S,int M,int Q,Y,normalization,pb):
    cdef double Y_dp = float(Y)
    cdef double pi = math.pi
    cdef double two_pi = 2*pi
    cdef double one_over_2Q = 1.0/(2*Q)
    cdef int weight = S.weight()
    cdef int weight_half = weight/2
    cdef double Y_pow_weight = Y_dp**weight_half
    G = S.group()
    cdef int nc = G.ncusps()
    V = np.zeros(shape=(nc*M,nc*M),dtype=np.complex_)
    fft_input = np.zeros(2*Q,dtype=np.complex_)
    fft_output = np.zeros(2*Q,dtype=np.complex_)
    cdef int cii,cjj,i,j
    cdef double c, d
    cdef double complex z_horo, z_fund, czd

    for cii in range(nc):
        for cjj in range(nc):
            coordinates = pb[cii][cjj]['coordinates_dp']
            coord_len = len(coordinates)
            j_values = pb[cii][cjj]['j_values']
            q_values = np.zeros(2*Q,dtype=np.complex_) #We use zero-padded values
            D_R = np.zeros(shape=(2*Q,1),dtype=np.complex_)
            y_facts = np.zeros(2*Q,dtype=np.double) #We use zero-padded values

            Msjj = len(normalization[cjj])+1
            Mfjj = Msjj+M-1
            Msii = len(normalization[cii])+1
            Mfii = Msii+M-1

            for i in range(coord_len):
                j_pos = j_values[i] #position in zero-padded vector
                (z_horo,c,d,z_fund) = coordinates[i]
                czd = c*z_horo+d
                D_R_fact = one_over_2Q*(cabs(czd)/czd)**weight
                if Msii != 0:
                    D_R_fact *= cexp(-one_over_2Q*pi*Msii*(2*j_pos+2*Q+1)*1j)
                D_R[j_pos,0] = D_R_fact
                q_values[j_pos] = cexp(two_pi*(z_fund*1j+Y_dp)) #We already divide by exp(-2*pi*Y) here
                y_facts[j_pos] = (cimag(z_fund)/Y_dp)**weight_half #We already divide by Y^(k/2) here
            max_abs_pos = np.abs(q_values).argmax() #Return position of entry with largest absolute value
            M_trunc = min(int(math.ceil(math.log(1e-16/y_facts[max_abs_pos])/math.log(cabs(q_values[max_abs_pos])))), M)

            W = np.zeros(shape=(2*Q,M_trunc),dtype=np.complex_,order='F')
            #Now compute first column of W
            if Msjj == 0: #q^0 = 1
                W[:,0] = y_facts[:]
            elif Msjj == 1:
                W[:,0] = y_facts[:]*q_values[:]
            else:
                W[:,0] = y_facts[:]*(q_values[:]**Msjj)
            #Compute remaining columns recursively
            for i in range(1,M_trunc):
                W[:,i] = W[:,i-1]*q_values[:]

            #Left multiply diagonal matrix
            W *= D_R

            #Perform FFT over all columns
            fft_res = np.fft.fft(W,axis=0)

            #Only consider first M rows
            res = fft_res[:M,:]

            #Compute D_L
            D_L = np.empty(shape=(M,1),dtype=np.complex_)
            for i in range(M):
                D_L[i,0] = cexp(-pi*1j*i*(0.5-Q)/Q)
            
            #Left multiply diagonal matrix
            res *= D_L

            #Filter out all elements < 1e-16
            res.real[abs(res.real) < 1e-16] = 0.0
            res.imag[abs(res.imag) < 1e-16] = 0.0

            #Store result in V
            V[cii*M:(cii+1)*M,cjj*M:cjj*M+M_trunc] = res

    #Subtract unit diagonal
    for i in range(nc*M):
        V[i,i] -= 1

    return V

cdef class PLU_Mat():
    """
    Computes and stores the LU decomposition of input matrix A at precision prec
    """
    def __cinit__(self, A, int prec):
        cdef int i, j, nrows, ncols
        cdef double tmp_double
        cdef Acb_Mat acb_mat_cast
        cdef double complex [::1, :] L, U #The result of scipy lu will be column major
        cdef int [:] piv

        if isinstance(A, Acb_Mat):
            acb_mat_cast = A
            nrows = acb_mat_nrows(acb_mat_cast.value)
            ncols = acb_mat_ncols(acb_mat_cast.value)
            sig_str("Arb exception")
            acb_mat_init(self.value, nrows, ncols)
            sig_off()
            self.P = <long*>sig_malloc(nrows*sizeof(long))
            sig_on()
            lu_successful = acb_mat_approx_lu(self.P, self.value, acb_mat_cast.value, prec)
            sig_off()
            if lu_successful == 0:
                raise ArithmeticError("LU-decomposition failed, please use higher precision!")
                
        elif isinstance(A, np.ndarray):
            nrows, ncols = A.shape
            A_sp = csc_matrix(A)
            lu_sp = sla.splu(A_sp, permc_spec="NATURAL", diag_pivot_thresh=0, options={"SymmetricMode":True}) #See: https://github.com/scipy/scipy/issues/7700#issuecomment-441131694

            sig_str("Arb exception")
            acb_mat_init(self.value, nrows, ncols)
            sig_off()
            self.P = <long*>sig_malloc(nrows*sizeof(long))
   
            L = lu_sp.L.A
            for i in range(1, nrows):
                for j in range(i):
                    tmp_double = creal(L[i,j])
                    if abs(tmp_double) > 1e-17: #If result is less than rounding error, we can set it to zero for faster multiplication later
                        arf_set_d(arb_midref(acb_realref(acb_mat_entry(self.value,i,j))), tmp_double)
                    tmp_double = cimag(L[i,j])
                    if abs(tmp_double) > 1e-17: #If result is less than rounding error, we can set it to zero for faster multiplication later
                        arf_set_d(arb_midref(acb_imagref(acb_mat_entry(self.value,i,j))), tmp_double)

            U = lu_sp.U.A
            for i in range(nrows):
                for j in range(i, ncols):
                    tmp_double = creal(U[i,j])
                    if abs(tmp_double) > 1e-17: #If result is less than rounding error, we can set it to zero for faster multiplication later
                        arf_set_d(arb_midref(acb_realref(acb_mat_entry(self.value,i,j))), tmp_double)
                    tmp_double = cimag(U[i,j])
                    if abs(tmp_double) > 1e-17: #If result is less than rounding error, we can set it to zero for faster multiplication later
                        arf_set_d(arb_midref(acb_imagref(acb_mat_entry(self.value,i,j))), tmp_double)

            P = self.P
            piv = lu_sp.perm_r
            for i in range(nrows):
                P[i] = piv[i]

    def __dealloc__(self):
        acb_mat_clear(self.value)
        sig_free(self.P)


