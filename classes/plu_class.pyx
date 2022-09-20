# ****************************************************************************
#       Copyright (C) 2022 David Berghaus <berghaus@th.physik.uni-bonn.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

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
from classes.acb_mat_class cimport Acb_Mat, Acb_Mat_Win

cdef class PLU_Mat():
    """
    Computes and stores the LU decomposition of input matrix A at precision prec
    """
    def __cinit__(self, A, int prec, use_scipy_lu):
        cdef int i, j, nrows, ncols
        cdef double tmp_double
        cdef Acb_Mat acb_mat_cast
        cdef double complex [::1, :] L, U #The result of scipy lu will be column major
        cdef int [:] piv

        if isinstance(A, Acb_Mat):
            self.use_scipy_lu = False #No scipy involved here anyways
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
            self.use_scipy_lu = use_scipy_lu
            nrows, ncols = A.shape
            A_sp = csc_matrix(A)
            lu_sp = sla.splu(A_sp, permc_spec="NATURAL", diag_pivot_thresh=0, options={"SymmetricMode":True}) #See: https://github.com/scipy/scipy/issues/7700#issuecomment-441131694

            if use_scipy_lu == True:
                self.lu_sp = lu_sp
            else: #Convert result to an arb matrix
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

    def solve(self, Acb_Mat x, Acb_Mat b, bit_prec):
        """
        Solve A*x=b by using a precomputed LU-decomposition.
        """
        if self.use_scipy_lu == False:
            self.solve_arb(x,b,bit_prec)
        else:
            self.solve_scipy(x,b)

    def solve_arb(self, Acb_Mat x, Acb_Mat b, int bit_prec):
        """
        Solve A*x=b by using a LU-decomposition that is stored as an arb matrix.
        """
        sig_on()
        acb_mat_approx_solve_lu_precomp(x.value,self.P,self.value,b.value,bit_prec)
        sig_off()
    
    def solve_arb_win(self, Acb_Mat_Win x, Acb_Mat_Win b, int bit_prec):
        """
        Solve A*x=b by using a LU-decomposition that is stored as an arb matrix.
        """
        sig_on()
        acb_mat_approx_solve_lu_precomp(x.value,self.P,self.value,b.value,bit_prec)
        sig_off()
    
    def solve_scipy(self, Acb_Mat x, Acb_Mat b):
        """
        Solve A*x=b by using the LU-decomposition that is stored as a sparse scipy matrix.
        We first multiply b by 2^scaling_exponent to put all results in double exponent range.
        Afterwards we perform the solve using sparse double arithmetic.
        We then convert the result back to arb and scale them back.
        DISCLAIMER: THIS FUNCTION ASSUMES THAT ALL ENTRIES OF x AND b ARE UNIFORMLY DISTRIBUTED SO THAT THEY CAN BE FIT INTO DOUBLES
        AFTER SCALING! IF THIS IS NOT THE CASE, PLEASE MAKE USE OF THE ARB-ROUTINES INSTEAD!
        """
        cdef int i
        largest_exponent = -9223372036854775807
        cdef int nrows = b.nrows()
        for i in range(nrows): #Detect smallest exponent of real parts of b entries
            b_i_exponent = arb_get_exponent(acb_realref(acb_mat_entry(b.value,i,0))) #In principle we need to considers abs here but real is sufficient for our needs
            if b_i_exponent > largest_exponent:
                largest_exponent = b_i_exponent
        scaling_exponent = -largest_exponent
        for i in range(nrows): #Scale exponents of b
            acb_mul_2exp_si(acb_mat_entry(b.value,i,0),acb_mat_entry(b.value,i,0),scaling_exponent)
        b_scaled_np = np.zeros(nrows,dtype=np.complex_)
        cdef double complex [:] b_scaled_np_view = b_scaled_np
        for i in range(nrows): #Convert b to numpy
            b_scaled_np_view[i] = arf_get_d(arb_midref(acb_realref(acb_mat_entry(b.value,i,0))),ARF_RND_NEAR) + arf_get_d(arb_midref(acb_imagref(acb_mat_entry(b.value,i,0))),ARF_RND_NEAR)*1j
        x_scaled_np = self.lu_sp.solve(b_scaled_np)
        cdef double complex [:] x_scaled_np_view = x_scaled_np
        for i in range(nrows): #Convert x to arb
            acb_set_d_d(acb_mat_entry(x.value,i,0),creal(x_scaled_np_view[i]),cimag(x_scaled_np_view[i]))
        for i in range(nrows): #Scale exponents back
            acb_mul_2exp_si(acb_mat_entry(b.value,i,0),acb_mat_entry(b.value,i,0),-scaling_exponent)
            acb_mul_2exp_si(acb_mat_entry(x.value,i,0),acb_mat_entry(x.value,i,0),-scaling_exponent)

    def __dealloc__(self):
        if self.use_scipy_lu == False:
            acb_mat_clear(self.value)
            sig_free(self.P)


