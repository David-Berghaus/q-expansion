from cysignals.signals cimport sig_on, sig_off
import numpy as np
import math
cimport numpy as np

cdef extern from "complex.h":
    cdef double cimag(double complex)
    cdef double creal(double complex)
    cdef double cabs(double complex)
    cdef double complex cexp(double complex)

from sage.libs.arb.acb_mat cimport *

from arblib_helpers.acb_approx cimport *
from classes.acb_mat_class cimport Acb_Mat, Acb_Mat_Win
from point_matching.point_matching_arb_wrap cimport _get_J_block_matrix_arb_wrap, _get_W_block_matrix_arb_wrap

cdef class J_class():
    def __init__(self, use_FFT):
        self.use_FFT = use_FFT
        self.is_initialized = False #This flag indicates if J has been fully initialized or not

    def _construct_non_fft(self,int M,int Ms,int Mf,int weight,int Q,coordinates,int bit_prec):
        """
        We compute the action of J by constructing it as a matrix.
        """
        self.J = Acb_Mat(M,len(coordinates))
        _get_J_block_matrix_arb_wrap(self.J.value,Ms,Mf,weight,Q,coordinates,bit_prec)
    
    def _construct_fft(self,int M,int Ms,int Mf,int weight,int Q,coordinates,int bit_prec):
        raise ArithmeticError("NOT IMPLEMENTED YET!")
    
    def _construct(self,int M,int Ms,int Mf,int weight,int Q,coordinates,int bit_prec,bint use_FFT):
        if use_FFT == False and self.use_FFT == False:
            self._construct_non_fft(M,Ms,Mf,weight,Q,coordinates,bit_prec)
        elif use_FFT == True and self.use_FFT == True:
            self._construct_fft(M,Ms,Mf,weight,Q,coordinates,bit_prec)
        else:
            raise ArithmeticError("Wrong initialization!")
        self.is_initialized = True

    cdef act_on_vec(self, acb_mat_t b, acb_mat_t x, int prec):
        cdef Acb_Mat two_Q_vec
        cdef Acb_Mat_Win post_fft_vec

        if self.is_initialized == True:
            if self.use_FFT == False:
                sig_on()
                acb_mat_approx_mul(b, self.J.value, x, prec)
                sig_off()
            else:
                raise ArithmeticError("Not implemented yet")
        else:
            raise ArithmeticError("J has not been properly initialized yet!")

cdef class W_class():
    def __init__(self, use_Horner):
        self.use_Horner = use_Horner
        self.is_initialized = False #This flag indicates if W has been fully initialized or not

    def _construct_non_horner(self,int M,int Ms,int Mf,int weight,coordinates,int bit_prec):
        self.W = Acb_Mat(len(coordinates),M)
        _get_W_block_matrix_arb_wrap(self.W.value,Ms,Mf,weight,coordinates,bit_prec)
    
    def _construct_horner(self,int M,int Ms,int Mf,int weight,coordinates,int bit_prec):
        raise ArithmeticError("NOT IMPLEMENTED YET!")
    
    def _construct(self,int M,int Ms,int Mf,int weight,coordinates,int bit_prec,bint use_Horner):
        self._nrows = len(coordinates)
        if use_Horner == False and self.use_Horner == False:
            self._construct_non_horner(M,Ms,Mf,weight,coordinates,bit_prec)
        elif use_Horner == True and self.use_Horner == True:
            self._construct_horner(M,Ms,Mf,weight,coordinates,bit_prec)
        else:
            raise ArithmeticError("Wrong initialization!")
        self.is_initialized = True

    cdef act_on_vec(self, acb_mat_t b, acb_mat_t x, int prec):
        if self.is_initialized == True:
            if self.use_Horner == False:
                sig_on()
                acb_mat_approx_mul(b, self.W.value, x, prec)
                sig_off()
            else:
                raise ArithmeticError("Not implemented yet")
        else:
            raise ArithmeticError("W has not been properly initialized yet!")
    
    def nrows(self):
        if self.is_initialized == True:
            return self._nrows
        else:
            raise ArithmeticError("W has not been properly initialized yet!")

cdef class Block_Factored_Element():
    def __cinit__(self, J_class J, W_class W):
        self.J = J
        self.W = W

cdef class Block_Factored_Mat():
    """
    This class stores a block-factored version of matrix V_tilde that can be used to compute
    the action of V_tilde on a vector without computing V_tilde explicitly.
    The class contains two variables:
    A : A two-dimensional list of size ncusps which later gets filled with matrices (J,W)
    diag : A list of size ncusps which later gets filled with Mx1 matrices that are diagonal elements
    diag_inv : A list of size ncusps which contains the inverses of diag
    """
    def __init__(self, int nc):
        A = [[None for cjj in range(nc)] for cii in range(nc)]
        diag = [None for cii in range(nc)]
        diag_inv = [None for cii in range(nc)]
        self.A = A
        self.diag = diag
        self.diag_inv = diag_inv
        self.nc = nc
        self.max_len = 0

        #Now to some parameters that we need to construct V_sc_dp
        self.S = None
        self.M = None
        self.Q = None
        self.Y = None
        self.normalization = None #We only need one representative of a 'normalization' but V can also be used for other types of normalizations
        self.pb = None
        self.is_cuspform = None
        self.parameters_for_dp_initialized = False
    
    def _init_parameters_for_dp_construction(self,S,int M,int Q,Y,normalization,pb,is_cuspform):
        """
        We construct V_sc_dp more or less from scratch so we need these parameters.
        """
        self.S = S
        self.M = M
        self.Q = Q
        self.Y = Y
        self.normalization = normalization #We only need one representative of a 'normalization' but V can also be used for other types of normalizations
        self.pb = pb
        self.is_cuspform = is_cuspform
        self.parameters_for_dp_initialized = True

    cpdef construct_sc_np(self):
        """
        Construct V_tilde_sc by using FFTs and truncate at 2e-16.
        We use this matrix to compute a preconditioner.
        We could in principle re-use some computations performed with arb by converting them to doubles but it seems faster
        to recompute everything from scratch using only doubles.
        """
        if self.parameters_for_dp_initialized == False:
            raise ArithmeticError("Class has not been properly initialized for this operation yet!")
        S = self.S
        cdef int M = self.M
        cdef int Q = self.Q
        cdef double Y_dp = float(self.Y)
        normalization = self.normalization
        pb = self.pb
        is_cuspform = self.is_cuspform
        cdef double pi = math.pi
        cdef double two_pi = 2*pi
        cdef double one_over_2Q = 1.0/(2*Q)
        cdef int weight = S.weight()
        cdef int weight_half = weight/2
        cdef double Y_pow_weight = Y_dp**weight_half
        G = S.group()
        cdef int nc = G.ncusps()
        V = np.zeros(shape=(nc*M,nc*M),dtype=np.complex_) #To reduce memory consumption, consider using https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.lil_matrix.html
        fft_input = np.zeros(2*Q,dtype=np.complex_)
        fft_output = np.zeros(2*Q,dtype=np.complex_)
        cdef int cii,cjj,i,j
        cdef double c, d
        cdef double complex z_horo, z_fund, czd, tmp
        if is_cuspform == True:
            coeff_shift = 1 #Cuspforms have c_0 = 1
        else:
            coeff_shift = 0

        for cii in range(nc):
            for cjj in range(nc):
                coordinates = pb[cii][cjj]['coordinates_dp']
                coord_len = len(coordinates)
                if coord_len != 0:
                    j_values = pb[cii][cjj]['j_values']
                    q_values = np.zeros(2*Q,dtype=np.complex_) #We use zero-padded values
                    D_R = np.zeros(shape=(2*Q,1),dtype=np.complex_)
                    y_facts = np.zeros(2*Q,dtype=np.double) #We use zero-padded values

                    Msjj = len(normalization[cjj])+coeff_shift
                    Mfjj = Msjj+M-1
                    Msii = len(normalization[cii])+coeff_shift
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
                    M_trunc = min(int(math.ceil(math.log(2e-16/y_facts[max_abs_pos])/math.log(cabs(q_values[max_abs_pos])))), M)

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
                    D_L[0,0] = 1
                    tmp = cexp(-pi*1j*(0.5-Q)/Q)
                    D_L[1,0] = tmp
                    for i in range(2,M):
                        D_L[i,0] = D_L[i-1,0]*tmp
                    
                    #Left multiply diagonal matrix
                    res *= D_L

                    #Filter out all elements < 2e-16 to make matrix more sparse
                    res.real[abs(res.real) < 2e-16] = 0.0
                    res.imag[abs(res.imag) < 2e-16] = 0.0

                    #Store result in V
                    V[cii*M:(cii+1)*M,cjj*M:cjj*M+M_trunc] = res

        #Subtract unit diagonal
        for i in range(nc*M):
            V[i,i] -= 1

        return V

    cpdef _get_max_len(self):
        """
        Returns the largest row-length of W (i.e. the highest amount of matching-points per cusp-pair)
        """
        cdef Block_Factored_Element block_factored_element
        if self.max_len == 0: #Not initialized yet so we need to compute it
            if self.diag[0] == None:
                raise NameError("Matrix is not properly initialized yet!")
            nc = self.nc
            A = self.A
            lengths = []

            for cii in range(nc):
                for cjj in range(nc):
                    if A[cii][cjj] != None:
                        block_factored_element = A[cii][cjj]
                        lengths.append((block_factored_element.W).nrows())

            self.max_len = max(lengths)
        
        return self.max_len

    cpdef act_on_vec_sc(self, Acb_Mat b, Acb_Mat x, int prec):
        """
        Computes Block_Factored_Mat*Diag_inv*x = b
        """
        cdef int nc, cii, cjj, M
        nc = self.nc
        A = self.A
        diag = self.diag
        diag_inv = self.diag_inv
        if x.value == b.value:
            raise NameError("Aliasing not allowed here!")
        if diag[0] == None:
            raise NameError("Matrix is not properly initialized yet!")
        M = diag[0].nrows() #diag[0] cannot be None if matrix is initialized
        cdef Acb_Mat tmp = Acb_Mat(self._get_max_len(), 1)
        cdef Acb_Mat_Win tmp_view
        cdef Acb_Mat tmp2 = Acb_Mat(M, 1)
        cdef J_class J
        cdef W_class W
        cdef Acb_Mat diag_cast
        cdef Acb_Mat_Win b_cast, x_cast
        cdef Block_Factored_Element block_factored_element

        #Slice vectors x & b into blocks
        x_blocks = []
        b_blocks = []
        for cii in range(nc):
            x_blocks.append(x.get_window(cii*M,0,(cii+1)*M,1))
            b_blocks.append(b.get_window(cii*M,0,(cii+1)*M,1))

        for cii in range(nc):
            b_cast = b_blocks[cii]
            x_cast = x_blocks[cii]
            #b[cii] = -x[cii]
            sig_on()
            acb_mat_set(b_cast.value, x_cast.value)
            sig_off()
            sig_on()
            acb_mat_neg(b_cast.value, b_cast.value)
            sig_off()
            for cjj in range(nc):
                if A[cii][cjj] != None:
                    block_factored_element = A[cii][cjj]
                    J = block_factored_element.J
                    W = block_factored_element.W
                    x_cast = x_blocks[cjj]
                    diag_cast = diag_inv[cjj]
                    #tmp2 = Diag_inv[cjj]*x[cjj]
                    acb_mat_approx_left_mul_diag(tmp2.value, diag_cast.value, x_cast.value, prec)
                    tmp_view = tmp.get_window(0, 0, W.nrows(), 1)
                    #tmp = W[cii][cjj]*tmp2
                    W.act_on_vec(tmp_view.value, tmp2.value, prec)
                    #tmp2 = J[cii][cjj]*tmp
                    J.act_on_vec(tmp2.value, tmp_view.value, prec)
                    #b[cii] += tmp2
                    sig_on()
                    acb_mat_approx_add(b_cast.value, b_cast.value, tmp2.value, prec)
                    sig_off()

    cpdef act_on_vec(self, Acb_Mat b, Acb_Mat x, int prec, is_scaled):
        """
        Computes action of self on vector x and stores result in b
        """
        if is_scaled == True:
            return self.act_on_vec_sc(b, x, prec)
        elif is_scaled == False:
            raise ArithmeticError("This functionality is not supported!")

    cpdef act_on_vec_win_sc(self, Acb_Mat_Win b, Acb_Mat_Win x, int prec):
        """
        Computes Block_Factored_Mat*Diag_inv*x = b
        """
        cdef int nc, cii, cjj, M
        nc = self.nc
        A = self.A
        diag = self.diag
        diag_inv = self.diag_inv
        if x.value == b.value:
            raise NameError("Aliasing not allowed here!")
        if diag[0] == None:
            raise NameError("Matrix is not properly initialized yet!")
        M = diag[0].nrows() #diag[0] cannot be None if matrix is initialized
        cdef Acb_Mat tmp = Acb_Mat(self._get_max_len(), 1)
        cdef Acb_Mat_Win tmp_view
        cdef Acb_Mat tmp2 = Acb_Mat(M, 1)
        cdef Acb_Mat diag_cast
        cdef J_class J
        cdef W_class W
        cdef Acb_Mat_Win b_cast, x_cast
        cdef Block_Factored_Element block_factored_element

        #Slice vectors x & b into blocks
        x_blocks = []
        b_blocks = []
        for cii in range(nc):
            x_blocks.append(x.get_window(cii*M,0,(cii+1)*M,1))
            b_blocks.append(b.get_window(cii*M,0,(cii+1)*M,1))

        for cii in range(nc):
            b_cast = b_blocks[cii]
            x_cast = x_blocks[cii]
            #b[cii] = -x[cii]
            sig_on()
            acb_mat_set(b_cast.value, x_cast.value)
            sig_off()
            sig_on()
            acb_mat_neg(b_cast.value, b_cast.value)
            sig_off()
            for cjj in range(nc):
                if A[cii][cjj] != None:
                    block_factored_element = A[cii][cjj]
                    J = block_factored_element.J
                    W = block_factored_element.W
                    x_cast = x_blocks[cjj]
                    diag_cast = diag_inv[cjj]
                    #tmp2 = Diag_inv[cjj]*x[cjj]
                    acb_mat_approx_left_mul_diag(tmp2.value, diag_cast.value, x_cast.value, prec)
                    tmp_view = tmp.get_window(0, 0, W.nrows(), 1)
                    #tmp = W[cii][cjj]*tmp2
                    W.act_on_vec(tmp_view.value, tmp2.value, prec)
                    #tmp2 = J[cii][cjj]*tmp
                    J.act_on_vec(tmp2.value, tmp_view.value, prec)
                    #b[cii] += tmp2
                    sig_on()
                    acb_mat_approx_add(b_cast.value, b_cast.value, tmp2.value, prec)
                    sig_off()

    cpdef act_on_vec_win(self, Acb_Mat_Win b, Acb_Mat_Win x, int prec, is_scaled):
        """
        Computes action of self on vector-window x and stores result in b
        """
        if is_scaled == True:
            return self.act_on_vec_win_sc(b, x, prec)
        elif is_scaled == False:
            raise ArithmeticError("This functionality is not supported!")

    cpdef nrows(self):
        """
        Returns full dimension of block-factored matrix
        """
        nc = self.nc
        diag = self.diag
        if diag[0] == None:
            raise NameError("Matrix is not properly initialized yet!")
        M = diag[0].nrows() #diag[0] cannot be None if matrix is initialized
        return M*nc

    cpdef diag_scale_vec(self, Acb_Mat res, Acb_Mat x, int prec):
        """
        Computes Diag*x where Diag consists of diagonal block-matrices
        """
        cdef int nc, cii, M
        nc = self.nc
        diag = self.diag
        if diag[0] == None:
            raise NameError("Matrix is not properly initialized yet!")
        M = diag[0].nrows() #diag[0] cannot be None if matrix is initialized
        cdef Acb_Mat diag_cast
        cdef Acb_Mat_Win res_cast, x_cast

        #Slice vectors res & x into blocks
        res_blocks = []
        x_blocks = []
        for cii in range(nc):
            res_blocks.append(res.get_window(cii*M,0,(cii+1)*M,1))
            x_blocks.append(x.get_window(cii*M,0,(cii+1)*M,1))

        for cii in range(nc):
            res_cast = res_blocks[cii]
            x_cast = x_blocks[cii] 
            diag_cast = diag[cii]
            acb_mat_approx_left_mul_diag(res_cast.value, diag_cast.value, x_cast.value, prec)

    cpdef diag_inv_scale_vec(self, Acb_Mat res, Acb_Mat x, int prec):
        """
        Computes Diag_inv*x where Diag_inv consists of diagonal block-matrices
        """
        cdef int nc, cii, M
        nc = self.nc
        diag_inv = self.diag_inv
        if diag_inv[0] == None:
            raise NameError("Matrix is not properly initialized yet!")
        M = diag_inv[0].nrows() #diag_inv[0] cannot be None if matrix is initialized
        cdef Acb_Mat diag_cast
        cdef Acb_Mat_Win res_cast, x_cast

        #Slice vectors res & x into blocks
        res_blocks = []
        x_blocks = []
        for cii in range(nc):
            res_blocks.append(res.get_window(cii*M,0,(cii+1)*M,1))
            x_blocks.append(x.get_window(cii*M,0,(cii+1)*M,1))

        for cii in range(nc):
            res_cast = res_blocks[cii]
            x_cast = x_blocks[cii] 
            diag_cast = diag_inv[cii]
            acb_mat_approx_left_mul_diag(res_cast.value, diag_cast.value, x_cast.value, prec)    





