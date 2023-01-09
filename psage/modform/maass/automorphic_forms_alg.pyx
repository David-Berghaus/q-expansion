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
import mpmath
from sage.functions.all import ceil

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

def err_est_holo(Y,M,weight):
    r"""
    The theoretical error estimate for the truncation error.
    """
    arg=mpmath.mp.pi()*mpmath.mp.mpf(2)*mpmath.mp.mpf(Y)
    mm=mpmath.mp.mpf(M+1)
    wt=mpmath.mpf(weight+1)/mpmath.mp.mpf(2)
    res=mpmath.mp.power(arg,-wt)
    res=2*res*mpmath.mp.gammainc(wt,arg*mm)
    return res