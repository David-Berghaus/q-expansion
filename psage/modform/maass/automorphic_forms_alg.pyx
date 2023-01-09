# cython: profile=False
# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2010 Fredrik Strömberg <stroemberg@mathematik.tu-darmstadt.de>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
r"""
Cython algorithms for Harmoic weak Maass forms.
Used by routines in atomorphic_forms.py

NOTICE: This code is part of psage https://github.com/fredstro/psage/blob/master/psage and has been copied with permission of the license holder.
We copy it here to remove the dependency on the installation of psage which can be a bit tricky sometimes.
"""
from libc.stdint cimport uint64_t
from psage.rings.mpfr_nogil cimport *
from cysignals.memory cimport sig_free,sig_malloc,check_allocarray
from cysignals.signals cimport sig_on,sig_off,sig_check

import logging
log = logging.getLogger(__name__)
cdef extern from "stdio.h":
    cdef extern void printf(char *fmt,...) nogil
    
#from sage.libs.mpfr cimport *
cdef mpc_rnd_t rnd
cdef mp_rnd_t rnd_re
import cython
from cython.parallel cimport parallel, prange
rnd = MPC_RNDNN
rnd_re = MPFR_RNDN
from sage.rings.complex_mpc cimport MPComplexNumber,MPComplexField_class
from sage.rings.complex_mpc import MPComplexField
from sage.rings.real_mpfr cimport RealNumber,RealField_class
from sage.rings.real_mpfr import RealField
from sage.rings.complex_number cimport ComplexNumber
import sage.structure.element
#from sage.libs.pari.gen cimport GEN
#cimport sage.structure.element
from psage.modules.vector_real_mpfr_dense cimport Vector_real_mpfr_dense
from psage.functions.inc_gamma cimport incgamma_hint_c,incgamma_nint_c,incgamma_pint_c
from psage.rings.mpc_extras cimport _mpc_div,_mpc_mul,_mpc_set,_mpc_sub,_mpc_mul_fr,_mpc_add

from sage.structure.element cimport Element, ModuleElement, RingElement
from sage.all import MatrixSpace,vector,ComplexField
cdef extern from "math.h" nogil:
    cdef double fabs(double)
    cdef double fmax(double,double)
    cdef int ceil(double) 

from sage.functions.all import ceil as pceil
from sage.functions.all import floor
#load "
from maass_forms_alg import *
#from maass_forms_alg cimport *

from psage.modform.arithgroup.mysubgroups_alg import normalize_point_to_cusp_mpfr as normalize_point_to_cusp
from psage.modform.arithgroup.mysubgroups_alg cimport _apply_sl2z_map_mpfr
from psage.modform.arithgroup.mysubgroups_alg import apply_sl2z_map,pullback_to_psl2z_mat #,pullback_to_G
#from maass_forms_alg import smallest_inf_norm
from sage.modular.arithgroup.congroup_sl2z import SL2Z
import mpmath
from sage.all import SageObject,copy,deepcopy
from pullback_algorithms import pullback_pts_mpc,j_fak_int_mpfr,pullback_pts_mp
from psage.matrix.matrix_complex_dense cimport Matrix_complex_dense
from psage.modules.vector_complex_dense cimport Vector_complex_dense
from pullback_algorithms cimport pullback_pts_mpc_new_c
from pullback_algorithms cimport pullback_pts_mpc_new_c_sym

cimport openmp

openmp.omp_set_dynamic(1)#openmp.omp_set_num_threads(2)

def get_Y_and_M_for_hwmf(G,PP,weight,ndig):
    r"""
    Find a good Y and M for computing coefficients with precision 10^-ndig
    """
     # generalized_level
    Y0=G.minimal_height()*mpmath.mpf(0.98)
    # then get corresponding M
    #M0=get_M_for_hwmf(Y0,Kmax,Cmax,weight,ndig)
    M0=get_M_for_hwmf(Y0,weight,ndig,PP)
    return [Y0,M0]


def get_M_for_holom(Y,weight,prec=10):
    r""" Get the correct M needed to approximate a  holomorphic
    modular form. 
    """
    # to get the range of M we check we make some crude
    # estimates to begin with
    if weight <= 0:
        raise NotImplementedError,"We have not implemented holomorphic forms of weight <= 0!"
    elif(abs(weight-mpmath.mp.mpf(1))>mpmath.mp.eps()):
        t1=mpmath.mp.pi()*mpmath.mp.mpf(2)*mpmath.mpf(Y)/mpmath.mp.mpf(weight-1)
    else:
        t1=mpmath.mp.mpf(0.2)
        #print "t1=",t1
    x=mpmath.mp.pi()*mpmath.mpf(Y)
    q=mpmath.exp(-x)
    t0=(prec-mpmath.mp.ln(mpmath.mp.mpf(1)-q))/x
    mstart=max(ceil(t0)+1,10) # A crude estimate of good M
    #print "mstart=",mstart
    eps=mpmath.mp.mpf(0.5)*mpmath.power(10,-prec)
    for M in range(mstart,3*mstart):
        t0=mpmath.mp.ln(M)/mpmath.mp.mpf(M)
        #print "t0(",M,")=",t0
        if(t0 > t1):
            continue
        t2=err_est_holo(Y,M,weight)
        #print "t2(",M,")=",t2
        if(t2<eps):
            #print "t2:%s < eps:%s for M:%s" %(t2,eps,M)
            break
    if(M>=3*mstart-1):
        t2=err_est_holo(Y,M,weight)
        if(t2>eps):
            raise ArithmeticError,"Could not find good M for Y:%s!" % Y
    #print "Got M=",M
    return M