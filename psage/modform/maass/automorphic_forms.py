# -*- coding: utf-8 -*-
#*****************************************************************************
#  Copyright (C) 2010 Fredrik Strömberg <fredrik314@gmail.com>
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
Implements Spaces of automorphic forms, for example Harmonic Weak Maass forms.
AUTHORS:
- Fredrik Strömberg

NOTICE: This code is part of psage https://github.com/fredstro/psage/blob/master/psage and has been copied with permission of the license holder.
We copy it here to remove the dependency on the installation of psage which can be a bit tricky sometimes.
"""

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division

from past.builtins import cmp
from builtins import str
from builtins import range
import mpmath
from sage.all import SageObject,Parent,ln,latex,random,divisors,ModularForms,prime_divisors,real,imag,PowerSeriesRing,\
    PolynomialRing,CyclotomicField,dimension_cusp_forms,dimension_modular_forms,CuspForms,ZZ,RealField,DirichletGroup,\
    Gamma0,trivial_character,Infinity,QQ,valuation,is_even,kronecker,SL2Z,CC,sign,copy,Integer,RR,identity_matrix, \
    matrix,MPComplexField,MatrixSpace,log_b,exp
from mpmath import mpf
from psage.modform.arithgroup.mysubgroup import *
from .automorphic_forms_alg import *
from sage.all import I,dumps,loads,ComplexField,LaurentPolynomialRing,next_prime,lcm
from sage.rings.fast_arith import prime_range
from sage.functions.other import real

class AutomorphicFormSpace(Parent):
    r"""
    General class of automorphic forms.
    Subclasses shoul
    d specialize to various types. 
    """
    def __init__(self,G,weight=0,multiplier="",character=0,holomorphic=False,weak=True,almost_holomorphic=False,cuspidal=False,unitary_action=0,dprec=15,prec=53,verbose=0,**kwds):
        r""" Initialize the space of automorphic forms.
        """
        self._from_group = None # try to keep the group used to construct the MyGroup instance
        if isinstance(G,(MySubgroup_class,HeckeTriangleGroup)):
            self._group=G
        elif is_int(G):
            self._from_group = Gamma0(G)
            self._group=MySubgroup(self._from_group)
        elif str(type(G)).find("gamma")>0 or str(type(G)).find("SL2Z")>0:
            self._from_group = G
            try:
                self._group=MySubgroup(G)
            except TypeError:
                raise TypeError("Incorrect input!! Need subgroup of PSL2Z! Got :{0}".format(G))
        else:
            raise TypeError("Could not convert G:{0} to a group!".format(G))
        self._unitary_action=unitary_action
        self._sym_type=None
        ## Define the character
        # And then the multiplier system, which should be an instance of MultiplierSystem or Subclass thereof, or None.
        if isinstance(character,sage.modular.dirichlet.DirichletCharacter):
            self._character = character
        elif is_int(character) and self._group.is_Gamma0():
            DG = DirichletGroup(self.level())
            if character >0 and character <len(DG): 
                self._character = DG[character]
            else:
                if verbose>0:
                    print("got character={0} as input!".format(character))
                self._character = trivial_character(1)
        elif character==0:
            self._character = trivial_character(1)
        else:
            raise TypeError("Could not find character {0} on group {1}".format(character,self._group))

        if not multiplier or multiplier=='':
            self._multiplier = TrivialMultiplier(self._group,character=self._character)
        elif isinstance(multiplier,MultiplierSystem):
            self._multiplier=multiplier
            self._character = multiplier._character
        else:
            raise TypeError("Incorrect multiplier! Got: {0}".format(multiplier))
        self._rdim=self._multiplier._dim
        # We assume weights are given as rational (integer of  half-integers)
        try:
            self._weight=QQ(weight)
        except:
            raise TypeError(" Need weights as rational numbers! Got:{0}".format(weight))
        # Check consistency of multiplier
        if not self._multiplier.is_consistent(self._weight):
            #print "mul=",self._multiplier
            #print "wt=",self._weight,type(self._weight)
            #print "even=",self._multiplier._character.is_even()
            #print "test=",self._multiplier.is_consistent(self._weight)
            #return self._multiplier
            raise ValueError(" The specified multiplier is not compatible with the given weight! \n multiplier:{0}, weight:{1}".format(self._multiplier,self._weight))

        self._dprec=dprec
        self._prec=prec
        if dprec>15:
            self._mp_ctx=mpmath.mp
            if prec==53:
                self._prec=int(3.4*dprec)+1
        else:
            self._mp_ctx=mpmath.fp
        #self._weight=mpmath.mpf(weight)
        self._holomorphic=holomorphic
        self._almost_holomorphic=almost_holomorphic
        self._weak=weak
        self._verbose=verbose
        self._cuspidal=cuspidal
        self._alphas={}
        self._dimension = -1
        # if we are interested in dimension of the corresponding cusp form space
        self._dimension_cusp_forms = -1
        self._dimension_modular_forms = -1
        self._basis_numerical = None
        # A properly working typing would make this unnecessary
        self._is_automorphic_form_space=True
        self._scaled=False
        if(self._multiplier.is_trivial()):
            self._rep=False
        else:
            self._rep=True
        # testing for types are tedious when in the interactive setting... 
        self._is_space_of_automorphic_functions=True
        self._do_mpmath = kwds.get("do_mpmath",0)
        
    def __repr__(self):
        r"""
        Return the string representation of self.
        EXAMPLES::
          sage: S=AutomorphicFormSpace(Gamma0(4));S
          Space of Automorphic Forms of weight = 0  with trivial multiplier  and character: Dirichlet character modulo 4 of conductor 1 mapping 3 |--> 1 on the group G:
          Arithmetic Subgroup of PSL2(Z) with index 6. Given by: 
          perm(S)=(1,2)(3,4)(5,6)
          perm(ST)=(1,3,2)(4,5,6)
          Constructed from G=Congruence Subgroup Gamma0(4)
          sage: S=AutomorphicFormSpace(Gamma0(4),multiplier=theta_multiplier,weight=1/2);S
          Space of Automorphic Forms of weight = 1/2  with theta multiplier  on the group G:
          Arithmetic Subgroup of PSL2(Z) with index 6. Given by: 
          perm(S)=(1,2)(3,4)(5,6)
          perm(ST)=(1,3,2)(4,5,6)
          Constructed from G=Congruence Subgroup Gamma0(4)
        
        """
        s="Space of "
        if self.is_holomorphic():
            if self.is_cuspidal():
                s+="Cusp Forms "
            else:
                s+="Modular Forms "
        elif str(type(self)).find("HarmonicWeakMaassFormSpace")>0:
            s+="Harmonic Weak Maass Forms "
        else:
            s+="Automorphic Forms "
        s+="of weight = "+str(self._weight)+" "
        if str(self._multiplier).find("theta_multiplier")>0:
            s+=" with theta multiplier "
        elif not self._multiplier.is_trivial():
            s+=" with multiplier:\n"+str(self._multiplier)
        else:
            s+=" with trivial multiplier "
        #if(self._character<>trivial and self._character<>trivial_character(self._group.level())):            
        #    s+=" and character: "+str(self._character)+" "
        s+=" on "
        if(self._from_group):
            s+=str(self._from_group)
        else:
            if self.group().index()==1:
                s+=" SL(2,Z)"
            else:
                s+="the group G\n"+str(self._group)
        return s


    def __reduce__(self):
        r"""
        """
        # we can't pickle functions so we store the names instead
        #if(isinstance(self._multiplier,type(trivial))):
        #    multiplier_s=self._multiplier.func_name            
        #else:
        #    multiplier_s=self._multiplier
        #if(isinstance(self._character,type(trivial))):
        #    character_s=self._character.func_name            
        #else:
        #    character_s=self._character
        if(self._from_group):
            G = self._from_group
        else:
            G = self._group
        return(AutomorphicFormSpace,(G,self._weight,self.multiplier(),self._holomorphic,self._weak,self._almost_holomorphic,self._cuspidal,self._dprec,self._verbose))

    def __eq__(self,other):
        r"""
        Compare self to other.
        """
        if self._verbose>0:
            print("in AutomorphicFormSpace.__eq__")
        if(not isinstance(other,type(self))):
            return False
        if(self._weight != other._weight):
            return False
        if(self._group != other._group):
            return False
        if(self._multiplier != other._multiplier):
            return False
        if(self._character != other._character):
            return False
        if(self._holomorphic != other._holomorphic):
            return False
        if(self._weak != other._weak):
            return False
        if(self._cuspidal != other._cuspidal):
            return False
        #eq = eq and (self._dprec == other._weak)
        #    return False
        #print "eq=",eq
        return True

    
    def __ne__(self,other):
        r"""
        Compare self to other.
        """
        return not self.__eq__(other)
    # return various properties of self

    def group(self):
        r""" Return the group of self.
        """
        return self._group

    def weight(self):
        r""" Return the weight of self.
        """
        return self._weight

    def character(self):
        r""" Return the character of self.
        """
        return self._multiplier._character

    def multiplier(self):
        r""" Return the multiplier of self.
        """
        return self._multiplier

    def sym_type(self):
        if self._sym_type==None:
            if hasattr(self._multiplier,"_sym_type"):
                self._sym_type=self._multiplier._sym_type
            else:
                self._sym_type=0
        return self._sym_type

    def prec(self,prec=None):
        if prec != None:
            self._prec=prec
        return self._prec

    def dprec(self,dprec=None):
        if dprec != None:
            self._dprec=dprec
        return self._dprec

    def is_holomorphic(self):
        r"""
        Return True if self is holomorphic, otherwise False. 
        """
        return self._holomorphic

    def is_almost_holomorphic(self):
        r"""
        Return True if self is almost holomorphic, otherwise False. 
        """
        return self._almost_holomorphic            

    def is_cuspidal(self):
        r"""
        Return True if self is cuspidal, otherwise False. 
        """
        return self._cuspidal



    def is_weak(self):
        r"""
        Return True if self is a weak form, i.e. has pole(s) at oo, otherwise False. 
        """
        return self._weak

    def is_harmonic(self):
        r"""
            Return True if self is a harmonic weak maassform, i.e. has pole(s) at oo, otherwise False. 
        """
        return self._harmonic


    def level(self):
        r""" Return the level of self (if self.group is a congruence subgroup).
        """
        if not self._group.is_congruence():
            raise ValueError("Level is only defined for congruence subgroups!")
        return self._group.generalised_level()

        
    def rep(self,A):
        r"""
        Calculate the representation = multiplier (including character) of self on the matrix A
        """
        v= self._multiplier(A)
        # duality and character should be builtin in the multiplier
        #    v = 1/v
        #v=v*self._character(A[1,1])
        return v

    def alpha(self,i):
        r"""
        Compute the translation at cusp nr. i, i.e. v(T_i)=e(alpha(i))
        -''i'' -- index of the cusp
        
        """
        RF=RealField(self._prec)
        CF=ComplexField(self._prec)
        if(i not in self._alphas):
            if(self._multiplier == None or self._multiplier.is_trivial()):
                self._alphas[i]=[RF(0),CF(1)]
            elif self.multiplier().ambient_rank()==1:
                #p=self._group._cusps[i]
                A=self._group._cusp_data[i]['stabilizer']
                tmp = self.rep(A)
                if hasattr(tmp,"complex_embedding"):
                    v=tmp.complex_embedding(self._prec)
                else:
                    v=CF(tmp)
                #ar=mpmath.arg(v)/mpmath.pi()/mpmath.mpf(2)
                ar=v.argument()/RF.pi()/RF(2)
                self._alphas[i]=[ar,v]
            else: ## This is nw a matrix-valued representation
                A=self._group._cusp_data[i]['stabilizer']
                mat = self.rep(A)[0]
                tmp = []
                for n in range(mat.nrows()):
                    a = mat[n,n]
                    if hasattr(a,"complex_embedding"):
                        v=a.complex_embedding(self._prec)
                    else:
                        v=CF(tmp)
                    ar=v.argument()/RF.pi()/RF(2)
                    ## Find the exact alpha
                    l = a.multiplicative_order()
                    if l < Infinity:
                        z = CyclotomicField(l).gens()[0]
                        for k in range(l):
                            if z**k == v:
                                break
                        tmp.append((ar,v,k,l))
                    else:
                        tmp.append((ar,v))
                self._alphas[i]=tmp
        return self._alphas[i]

    def alphas(self):
        return self._alphas

    def set_alphas(self,prec=None):
        r""" Compute the vector containing the shifts at the various cusps.
        """
        precold=self._prec
        if prec!=None and prec != self._prec:
            self._prec=prec
        for i in range(self._group._ncusps):
            self.alpha(i)
        self._prec=precold
    def dimension(self):
        r"""
        Return the dimension of self if we can compute it.
        """
        if self._dimension>=0:
            return self._dimension
        k=self._weight
        if self._holomorphic and self._rdim==1:
            if is_int(k) and self._multiplier.is_trivial():
                # Use the builtin sage routines
                xi = self.character()
                self._dimension_cusp_forms =  dimension_cusp_forms(xi,k)
                self._dimension_mod_forms =  dimension_modular_forms(xi,k)
            else:
                # In weight 1/2 use Serre-Stark
                if k==QQ(1)/QQ(2):
                    [d_mod,d_cusp]=self._dimension_of_weight_one_half_space()
                    if(self.is_cuspidal()):
                        return d_cusp
                    else:
                        return d_mod                
                elif k<QQ(1)/QQ(2):
                    self._dimension = 0
                    # Else use Cohen-Oesterle: (as described in The Web of Modularity by Ono)
                elif k>QQ(3)/QQ(2):
                    dimension_cusp_forms = self._difference_of_dimensions(k)
                    dimension_mod_forms = -self._difference_of_dimensions(QQ(2-k))
                elif k==QQ(3)/QQ(2):
                    [d_mod,d_cusp]=self._dimension_of_weight_one_half_space()                
                    dimension_cusp_forms = self._difference_of_dimensions(k)+d_mod
                    dimension_mod_forms = d_cusp - self._difference_of_dimensions(QQ(1)/QQ(2))
                # print "dim_cusp_forms=",dimension_cusp_forms
                # print "dim_mod_forms=",dimension_mod_forms            
                self._dimension_cusp_forms = dimension_cusp_forms
                self._dimension_mod_forms = dimension_mod_forms
                if(self.is_cuspidal()):
                    self._dimension = self._dimension_cusp_forms
                else:
                    self._dimension = self._dimension_mod_forms
        elif self._holomorphic:
            #if self._group.level()==1:
            #[d_mod,d_cusp]=self._dimension_of_vector_valued_forms()
            self._dimension = self._dimension_of_vector_valued_forms()
            
        else:
            raise NotImplementedError
                #self._dimension=-1
            
        return self._dimension

    def _difference_of_dimensions(self,k):
        r"""
        Use Cohen and Oesterle to compute the difference: dim(S_k)-dim(M_{2-k})
        """
        N = self._group.generalised_level()
        cond = self._character.conductor()
        #k = QQ(RR(self._weight)) 
        kk = ZZ(k - QQ(1)/QQ(2))
        r2 = valuation(N,2)
        s2 = valuation(cond,2)
        if r2>=4: #   zeta_k_l_chi = lambda_k_l_chi
            if 2*s2 <= r2:
                if is_even(r2):
                    rp = r2/QQ(2)
                    zeta_k_l_chi = 2**(rp) + 2**(rp-1)
                else:
                    rp = (r2-1)/QQ(2)
                    zeta_k_l_chi = 2*2**(rp)
            elif 2*s2 > r2:
                # TODO: Check formula!
                raise ArithmeticErro("This case is not implemented!")
                # zeta_k_l_chi = 2*2**(rp-sp)
        elif r2==3:
            zeta_k_l_chi = 3
        elif r2==2:
            zeta_k_l_chi = 0
            ## Condition (C)
            for p in prime_divisors(N):
                if (p % 4) == 3:
                    rp = valuation(N,p)
                    sp = valuation(cond,p)
                    if is_odd(rp) or (rp>0 and rp < 2*sp):
                        zeta_k_l_chi = 2
                        break
            if zeta_k_l_chi== 0: # not (C)
                if is_even(kk):
                    if s2==0:
                        zeta_k_l_chi = QQ(3)/QQ(2)
                    elif s2==2:
                        zeta_k_l_chi = QQ(5)/QQ(2)
                else:
                    if s2==0:
                        zeta_k_l_chi = QQ(5)/QQ(2)
                    elif s2==2:
                        zeta_k_l_chi = QQ(3)/QQ(2)
        if zeta_k_l_chi<=0:
            raise ArithmeticError("Could not compute zeta(k,l,chi)!")
        fak = QQ(1)
        for p in prime_divisors(N):
            fak = fak* QQ(1+QQ(1)/QQ(p))
        fak2 = QQ(1)
        for p in prime_divisors(N):
            if(p>2):
                rp = valuation(N,p)
                sp = valuation(cond,p)
                if(rp < 2*sp):
                    lam = QQ(2*p**(rp-sp))
                else:
                    if(is_even(rp)):
                        rprim=QQ(rp)/QQ(2)
                        lam = QQ(p**rprim)+QQ(p**(rprim-1))
                    else:
                        rprim=QQ(rp-1)/QQ(2)
                        lam = QQ(2*p**rprim)
                fak2 = fak2 * lam
            # S_k(Gamma0(N),chi) - M_{2-k}(Gamma0(N),chi)
        diff_dims = fak*QQ(k-1)*QQ(N)/QQ(12)-QQ(zeta_k_l_chi)/QQ(2)*fak2
        return diff_dims

    def _dimension_of_weight_one_half_space(self):
        r"""
        Computes the dimension of M and S_{1/2}(4N,chi) where 4N and chi are the level and character of self.
        
        """
        O = self._Omega()
        dim_mod_forms = len(O)
        # check number of totally even characters for the cusp forms
        nn = 0
        for (xi,t) in O:
            l = xi.decomposition()
            for xip in l:
                if xip(-1)==1:
                    nn = nn+1
        dim_cusp_forms = len(O) - nn
        return [dim_mod_forms,dim_cusp_forms]

    def _Omega(self):
        r"""
        Computes the set of pairs (psi,t) satisfying:
        (1) r**2 * t | self.level()/4
        (2) chi(n)= psi(n)*kronecker(t,n) for (n,self.level())=1
        """
        N = ZZ(QQ(self.level())/QQ(4))
        D = DirichletGroup(self.level())
        chi = self._character
        Omega=[]
        for t in divisors(N):
            for psi in D:
                r = psi.conductor()
                s = ZZ(r*r*t)
                if(not s.divides(N)):
                    continue
                ok = True
                for n in range(1,self.level()):
                    if(psi(n)*kronecker(t,n) != chi(n)):
                        ok = False
                        break
                if(ok):
                    Omega.append((psi,t))
        return Omega

    def _dimension_of_vector_valued_forms(self):
        r"""
        Calculates the dimension of self is self is a space of automorphic forms  on SL2(Z).
        """
        if self._dimension>=0:
            return self._dimension
        if self._weight < 2 or self.level() != 1:
            return -1
        try:
            if hasattr(self._multiplier,"dimension_cusp_forms"):
                if self._cuspidal:
                    self._dimension = self._multiplier.dimension_cusp_forms(self._weight)
                else:
                    self._dimension = self._multiplier.dimension_modular_forms(self._weight)
                return self._dimension
        except:
            pass
        term0 = QQ(self._rdim*(self._weight-1))/QQ(12)
        S,T=SL2Z.gens()
        R = S*T; R2=R*R
        if self._rdim>1:
            wS=self._multiplier(S)[0].trace()
            wR=self._multiplier(R)[0].trace()        
            wR2=self._multiplier(R2)[0].trace()
            evs = self._multiplier(T)[0].diagonal()
            wT=sum(evs) #self._multiplier(T).trace()
            alphas = [log(CC(x))/CC(2/pi) for x in evs]
        else:
            wS=self._multiplier(S)
            wR=self._multiplier(R)        
            wR2=self._multiplier(R2)            
            wT=self._multiplier(T)
        z2=CyclotomicField(4).gen()
        term1 = wS/QQ(4)*z2*z2**(self._weight-1)
        z3=CyclotomicField(6).gen()        
        term2 = wR/QQ(3)/sqrt(3)*z2*z3**(self._weight-1)
        term3 = wR/QQ(3)/sqrt(3)*z2*z3**(2*(self._weight-1))
        term4=0
        k0=0
        for a in alphas:
            if a != 0:
                term4+=a-0.5
            else:
                k0+=1
        term5 = k0/2*sign(self._weight-1)
        if self._cuspidal:
            term6 = 0
        else:
            term6 = k0
        if k0 != 0:
            if self._weight == 1:
                raise ArithmeticError("Need to compute the scattering determinant!")

        dim = term0 + term1 + term2 + term3 + term4 + term5 + term6
        return dim
        
    ## cuspidal subspace
    def cuspidal_subspace(self):
        r"""
        Construct the cuspidal subspace of self.        
        """
        S = copy(self) # copy self
        S._weak = False # not weak
        S._cuspidal = True # and cuspidal
        #S._holmorphic = True # and holmorphic in H
        #S=HalfIntegralWeightForms(G,self._weight,self._multiplier,character=self._character,holomorphic=self._holomorphic,weak=self._weak,cuspidal=True,dprec=self._dprec,verbose=self._verbose,construct=self._construct)
        return S



    def set_normalization(self,C=None):
        r"""
        -''C'' -- a dictionary of set coefficients in the form C[d][(i,n)]=c
        """
        #print "C0=",C
        N=dict()
        N['comp_dim']=1
        if isinstance(C,dict):
            N['comp_dim']=max(1,len(list(C.keys())))
        else:
            N['comp_dim']=max(1,len(C))
        N['SetCs']=dict()
        N['cuspidal']=self._cuspidal
        N['weak']=self._weak        
        nc=self._group.ncusps()
        for j in range(N['comp_dim']):
            N['SetCs'][j]=dict()
        #if(P.has_key((0,0)) and H._holo):
        #print "holo"
        #for j in range(N['comp_dim']):
        #    N['SetCs'][j][(0,0)]=0            
        if(N['cuspidal']):
            for icusp in range(nc):
                v=self.alpha(icusp)[1]
                if(v==1):
                    for j in range(N['comp_dim']):
                        N['SetCs'][j][(icusp,0)]=0
        if(not N['weak']):
            for icusp in range(nc):
                al=self.alpha(icusp)[0]
                if(al<-mpmath.eps()):
                    for j in range(N['comp_dim']):
                        N['SetCs'][j][(icusp,0)]=0
        if isinstance(C,dict) and C != {}:
            for i in list(C.keys()):
                for (r,n) in list(C[i].keys()):
                    N['SetCs'][i][(r,n)]=C[i][(r,n)]
        elif isinstance(C,list) and C != []:
            for i in range(len(C)):
                for (r,n) in list(C[i].keys()):
                    N['SetCs'][i][(r,n)]=C[i][(r,n)]
        return N

    def set_normalization_vv(self,P={},C={},c_t="pp"):
        r"""
        Set normalization for vector-valued case
        """
        N=dict()
        if isinstance(P,list):
            Pl=P
        else:
            Pl=[P]
        if isinstance(C,list):
            Cl=C
        else:
            Cl=[C]        
        if len(Pl)>0:
            N['comp_dim']=len(Pl)
        elif len(Cl)>0:
            N['comp_dim']=len(Cl)
        else:
            raise ValueError("Need either principal parts of set coefficients!")
        if len(Cl)>0:
            if len(Cl)!=len(Pl):
                raise ValueError("Need same number of principal parts and coefficients to set!")
            keys = list(Cl[0].keys())
            for j in range(1,N['comp_dim']):
                if list(Cl[j].keys())!=keys:
                    raise ValueError("Need to set the same coefficients! (or call the method more than once)")
        else:
            Cl=[]
            for j in range(N['comp_dim']):
                Cl.append(C)
        if self._verbose>0:
            print("Pl=",Pl)
            print("Cl=",Cl)
        N['Vals']=list()
        N['Vals']=list()
        N['SetCs']=list()
        for i in range(N['comp_dim']):
            N['Vals'].append({})
            N['Vals'].append({})
            N['SetCs'].append([])
            ## First look at all zeroth-coefficients and see if we are forced to set some to zero
            for j in range(self.multiplier().weil_module().rank()):
                a=self.multiplier().weil_module().basis()[j]
                x=self.multiplier().Qv[j]
                #N['Vals'][i][(0,j)]=dict()
                
                if x==0:
                    if c_t=="pp" and (0,j) in Pl[i]:
                        N['SetCs'][i].append((j,0))
                        N['Vals'][i][(j,0)]=Pl[i][(j,0)]
                    elif self._cuspidal:
                        N['SetCs'][i].append((j,0))
                        N['Vals'][i][(j,0)]=0 #P[(0,0)]
                    
                elif x<0 and self._holomorphic:
                    N['SetCs'][i].append((j,0))
                    N['Vals'][i][(j,0)]=0 #P[(0,0)]            

            if isinstance(Cl[i],dict):
                for (r,n) in list(Cl[i].keys()):
                    if(N['SetCs'][i].count((r,n))==0):
                        N['SetCs'][i].append((r,n))
                        N['Vals'][i][(r,n)]=Cl[i][(r,n)] 
        return N



    def get_Y_and_M(self,digs=10,principal_part=[{}]):
        r"""
        Get good values for Y and M to truncate with an error of prec digits
        """
        ## todo : more alternatives
        ## we use the largest term of the principal part
        ## to estimate the Y and M
        #print "pp0=",principal_part
        if not isinstance(principal_part,list):
            pp = [principal_part['+']]
        elif len(principal_part)==1:
            if '+' in principal_part[0]:
                pp = [principal_part[0]['+']]
            else:
                pp = [principal_part[0]]
        else:
            pp=list()
            for p in principal_part:
                pp.append(p['+'])
        maxn=0; maxa=1; maxr=0
        #print "pp1=",pp
        for P in pp: #rincipal_part:
            for (r,n) in list(P.keys()):
                # Remember that the principal part (i,j):c means different things for scalar and vector-valued forms, i.e. i is the cusp in the first case and the component in the second
                a = P[(r,n)]
                if(a>maxa):
                    maxa=a
                if self.multiplier().ambient_rank()>1:
                    if r >=0 and r< self.multiplier().rank():
                        aln = n + self.alpha(0)[r][0]
                        if aln < maxn:
                            maxr = r
                            maxn = aln
                else:
                    aln = n + self.alpha(r)[0]
                    if aln < maxn:
                        maxr = r
                        maxn = aln

        if self._verbose > 1:
            print("maxa={0}".format(maxa))
            print("maxn={0}".format(maxn))
            print("digits={0}".format(digs))
        if not self._holomorphic or self._weak:
            maxa = maxa * len(pp)
            pp_max = {(maxr,maxn):maxa}
            if self._verbose > 1:
                print("pp_max={0}".format(pp_max))
            #[Y,M]=self.get_Y_and_M(prec,principal_part=pp)
            [Y,M]=get_Y_and_M_for_hwmf(self._group,pp_max,self._weight,digs)
        else:
            Y=self._group.minimal_height()*mpmath.mpf(95)/mpmath.mpf(100)
            M=get_M_for_holom(Y,self._weight,digs)
        return [Y,M]

    ## def get_element(self,principal_part={},C=None,prec=10,M0_in=None,Y_in=None):
    ##     r"""
    ##     Get an element of self given by either principal part
    ##     or normalization of coefficients.


    ##     """
    ##     print "self.type=",type(self)
    ##     F=AutomorphicFormElement(self,principal_part=principal_part)
    ##     if(Y_in<>None and M0_in<>None):
    ##         Y=Y_in
    ##         M=M0_in
    ##     elif(Y_in<>None):
    ##         Y=Y_in
    ##         if(not self._holomorphic):
    ##             M=get_M_for_hwmf(Y,self._weight,prec,principal_part)
    ##         else:
    ##             M=get_M_for_holom(Y,self._weight,prec)               
    ##     elif(M0_in<>None):
    ##         M=M0_in
    ##         Y=self._group.minimal_height()*mpmath.mpf(95)/mpmath.mpf(100)
    ##     else:
    ##         [Y,M]=self.get_Y_and_M(prec,principal_part)
    ##     #Y=Y*0.7
    ##     Q=M+30
    ##     Ymp=mpmath.mp.mpf(Y)
    ##     if(not (principal_part.has_key('+') or principal_part.has_key('-'))):
    ##         raise ValueError,"Need principal part with keys '+' and '-'!"
    ##     PP=principal_part
    ##     V=setup_matrix_for_harmonic_Maass_waveforms_sv(self,Ymp,M,Q,PP)
    ##     V['PP']=PP
    ##     print "Use M,Y=",M,Y
    ##     # recall that the zeroth coefficient is counted in the principal part
    ##     return V
    ##     if(C==None):
    ##         C=dict()
    ##     for j  in range(len(self._group.cusps())):
    ##         al=self.alpha(j)
    ##         print "al(",j,")=",al
    ##         if(al[1]==1):
    ##             if(PP.has_key((j,0))):
    ##                 for r in C.keys():
    ##                     C[r]=dict()
    ##                     C[r][(j,0)]=0
    ##     print "C=",C
    ##     N=self.set_normalization(C)
    ##     print "N=",N
    ##     D=solve_system_for_harmonic_weak_Maass_waveforms(V,N,deb=True)
    ##     F._coeffs=D
    ##     return F
    def _get_element(self,principal_part,digs=10,dbase_prec=None,SetC=None,SetY=None,SetM=None,SetQ=None,do_mpmath=0,get_mat=False,use_sym=1,get_c=False,gr=0,version=0,threads=1,**kwds):
        r"""
        INPUT:
        
        - `principal_part`   -- list of principal parts of the form:
                 RR = { '+' : {(j,n) : c^+(j,n)}     # j is a cusp and n>=0 an index
                        '-' : {(j,n) : c^-(j,n)}     # j is a cusp and n<=0 an index
                     }
                corresponding to principal parts (in notation of Bruinier-Funke):
                    \(\Sum_{n>0} c^+(j,n)q^{-n} +  \Sum_{n<0} c^-(j,n)H(n\tau)
        
                    PP[c,m]=a if the principal at cusp c contains a*q^m
        - `digs` -- integer (default 10): the number of requested digits
        - `dbase_prec` -- integer (default None): if set, use this number of digits for precision in all mpmath calculations
        - `SetC` -- dictionary containing fourier coefficients to keep fixed (and their values)
                      of the form SetC[n][i]=c_i(n)
        """
        from .vv_harmonic_weak_maass_forms import solve_system_for_vv_harmonic_weak_Maass_waveforms_new

        ## the principal part and the SetC should be lists if present
        if self._verbose>0:
            print("PP={0}".format(principal_part))
            print("gr={0}".format(gr))
        if(not isinstance(principal_part,list)):
            ppart = [principal_part]
        else:
            ppart = principal_part
        ppart1=list()
        for pp in ppart:
            d=dict()
            d['-']=pp.get('-',{}) # By default we have no princ. part
            d['+']=pp.get('+',{}) 
            #print "pp=",pp
            if isinstance(list(pp.keys())[0],(list,tuple)):
                d['+']=pp # If only one is given we assume it holomorphic
            # If self._holomorphic is True and we have a negative principal part we assume
            # that the only non-holomorphic part is the principal part
            #if self._holomorphic: # Make sure no non-holom. ppart is given for a holom. form
            #    d['-']={}
            ppart1.append(d)
        if self._verbose>0:
            print("PP1={0}".format(ppart1))
        ppart = ppart1 #principal_part

        ## Check whether the set coefficients are the same for all elements
        ## if thy are the same we only need to solve the system once (using LU decomposition).
        ## otherwise we need to rerun the system solving several times
        
        if SetC!=None and not isinstance(SetC,list):
            setc=list()
            for i in range(len(ppart)):
                setc.append(SetC)
        elif not SetC:
            setc = []
        else:
            setc=SetC
        if len(setc)>0 and len(setc)!=len(ppart):
            raise ValueError("Inconsistent lengths of principal part and set coefficients!")        
        # recall that we treat 0-coefficients in the principal part
        # as variables.
        if self._verbose>0:
            print("setc0={0}".format(setc))
            #print "group=",self.group()
        for i in range(len(ppart)):
            for j in range(self.group().ncusps()):
                if (j,0) in ppart[i]['+']:
                    for ii in range(len(setc),i+1):
                        setc.append({})
                    setc[i][(j,0)]=ppart[i]['+'][(j,0)]
        if self._verbose>0:
            print("setc1={0}".format(setc))

        solve_once = True
        if setc != None and len(setc) > 0:
            d = list(setc[0].keys())
            for c in setc:
                if list(c.keys()) != d:
                    solve_once=False
                    break
#                    raise ValueError," Inconsistent set coefficients. Need the same length in all entries! Got: %s" %(setc)
        if self._verbose>0:
            print("solve_once={0}".format(solve_once))
            print("ppart={0}".format(ppart))
            print("setc={0}".format(setc))
        #pos_part=list()
        #for pp in ppart:
        #   pos_part.append(pp['+'])
        if(not (SetY and SetM)):
            [Y,M]=self.get_Y_and_M(digs,ppart)
        # print "dps=",mpmath.mp.dps
        #Y=Y*0.5 #97.
        if SetY != None:
            Y=SetY
        if SetM != None:
            M=SetM        
        if SetQ != None and SetQ>M:
            Q = SetQ
        else:
            Q=M+10

        mpmath.mp.prec = self._prec
        
        Ymp=RealField(self._prec)(Y) #mpmath.mp.mpf(Y)
        dpold=mpmath.mp.dps
        #if dbase_prec<>None:
        #    mpmath.mp.dps=max(self._dprec,dbase_prec)
        #else:
        #    mpmath.mp.dps=self._dprec
        sv = 1
        d = self.multiplier().rank()

        if d>1 or  hasattr(self.multiplier(),"D"):
            sv=0
        if self._verbose>0:
            print("dps={0}".format(mpmath.mp.dps))
            print("setc={0}".format(setc))
            print("Y={0}".format(Ymp))
            print("M={0}".format(M))
            print("Q={0}".format(Q))
            print("PP={0}".format(ppart))
            print("do_mpmath={0}".format(do_mpmath))
            print("alphas={0}".format(self.alphas()))
            print("dim={0}".format(d))
            print("scalar={0}".format(sv))
        C = None

        if sv==1:
            if do_mpmath==1:
                Ymp = mpmath.mp.mpf(Ymp)
                V=setup_matrix_for_harmonic_Maass_waveforms_sv_bak(self,Ymp,M,Q,ppart)
            elif do_mpmath==2:
                Ymp = mpmath.mp.mpf(Ymp)
                V=setup_matrix_for_harmonic_Maass_waveforms_sv_bak_22(self,Ymp,M,Q,ppart)
            elif (version==0 or gr==1):
                V=setup_matrix_for_harmonic_Maass_waveforms(self,Ymp,M,Q,ppart,use_sym=use_sym,threads=threads)
            else:
                C = setup_and_solve_for_harmonic_Maass_waveforms(self,Ymp,M,Q,ppart,cset=setc)
        else:
            ## Only one principal part is implemented in vv-case
            if self._holomorphic:
                V=vv_holomorphic_setupV_mpc(self,Ymp,M,Q)
            else:
                V=vv_harmonic_wmwf_setupV_mpc2(self,ppart[0],Ymp,M,Q)
            V['space']=self
            V['PP']=ppart
        if gr==1:
            return V
        if solve_once and C==None:
            #N=set_norm_harmonic_weak_maass_forms(self,ppart,setc)
            #if isinstance(setc,list):
            if sv==1:                
                if do_mpmath==0:
                    N = self.set_normalization(setc)
                else:
                    N = self.set_norm(ppart,setc)
            else:
                N = self.set_normalization_vv(ppart,setc)
            #else:
            #    N = self.set_normalization(setc)
            V['PP']=ppart
            #return V,N
            if self._verbose>0:
                print("N={0}".format(N))
            if do_mpmath!=0:
                C=solve_system_for_harmonic_weak_Maass_waveforms_mpmath(V,N)
            else:
                if sv==1:
                    C=solve_system_for_harmonic_weak_Maass_waveforms(V,N)
                else:
                     C=solve_system_for_vv_harmonic_weak_Maass_waveforms_new(self,V,N)

        elif C==None:
            C=list()
            RHS=V['RHS']
            if RHS.cols != len(ppart):
                raise ValueError("Inconsistent lengths of principal part and right hand sides!")        
            for i in range(len(ppart)):
                pp=[ppart[i]]; cc=[setc[i]]
                V['PP']=pp
                V['RHS']=RHS.column(i)
                if self._verbose>1:
                    print("cc={0}".format(cc))
                    print("pp={0}".format(pp))
                N = self.set_norm(ppart,setc)
                #N=set_norm_harmonic_weak_maass_forms(self,pp,cc)
                if self._verbose>1:
                    print("N={0}".format(N))
                #return V,N
                try:
                    C.append(solve_system_for_harmonic_weak_Maass_waveforms(V,N)[0])
                except:
                    return C.append((V,N))
        #mpmath.mp.dps=dpold
        if self._verbose>0:
            print("C[0][-1]=",C.get(0,{}).get(0,{}).get(-1,None))

        res=list()
        if get_c:
            return C
        prec = mpmath.mp.dps
        if len(C)>0:
            for i in range(len(C)):
                #print "PP=",ppart[i]
                ppf=dict()
                ppf['+']=ppart[i]['+']
                ppf['-']=ppart[i]['-']
                if self._verbose>1:
                    print("type={0}".format(type(self)))
                if str(type(self)).find("HalfIntegralWeightForms")>0:
                    F=HalfIntegralWeightFormElement(self,C[i],principal_part=ppf)
                elif str(type(self)).find("HarmonicWeakMaassFormSpace")>0:
                    if self._verbose>1:
                        print("Constructing a Harmonic Weak Maassform")
                        print("pp={0}".format(ppf))
                    F=HarmonicWeakMaassFormElement(self,C[i],prec=prec,principal_part=ppf)
                else:
                    F=AutomorphicFormElement(self,C[i],prec=prec,principal_part=ppf)

                #print "M0=",M
                F._M0 = M
                res.append(F)
                #print "appended"
                #print "res=",res
        if len(res)==1:
            return res[0]
        else:
            return res

    def smallest_M0(self):
        r"""
        Smallest M0 which we can use if we want to test using Hecke relations.
        """
        if is_Hecke_triangle_group(self._group):
            if self._group.is_Gamma0():
                self._smallest_M0=int(12)
            else:
                self._smallest_M0=int(12*self._group._lambdaq)
        if self._smallest_M0>0:
            return self._smallest_M0
        a = self.get_primitive_p()
        b = self.get_primitive_p(a)
        c = a*b
        self._smallest_M0=c+3
        return self._smallest_M0


        
    def get_primitive_p(self,p0=0,notone=1):
        r"""
        Gives a prime p to use for Hecke operator on M
        p should be relative prime to the level of M._group
        and to the modulus of M._multiplier._character
        INPUT:
        - 'p0' -- return prime greater than p0
        - 'notone' -- if set to one we return a prime with chi(p)<>1
        """
        if not self._group._is_congruence:
            return next_prime(p0)
        m=self._multiplier
        x=m._character
        if hasattr(x,"modulus"):
            modulus=x.modulus()
        else:
            modulus=1
        prim_to=lcm(self.level(),modulus)
        p00 = next_prime(p0)
        p01 = p00 + prim_to
        if notone:
            if self.level() % 9 ==0 :
                pq=3
                # celif self._group._level % 4 ==0 :
                #    pq=4
            else:
                pq=1

        for p in prime_range(p00,p01+1):
            if notone==1 and p%pq==1:
                continue
            if gcd(p,prim_to)==1:
                return p
        raise ArithmeticError(" Could not find appropriate p rel. prime to {0}!".format(prim_to))