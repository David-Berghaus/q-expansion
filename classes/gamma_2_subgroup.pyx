r"""
A general class for subgroups of Gamma(2).
AUTHORS:
 - Sage Group: Base class for SL2Z subgroups
 - Fredrik Strömberg: Modifications for computing cusp forms of SL2Z subgroups
 - David Berghaus: Modifications for computing cusp forms of Gamma(2) subgroups
"""

#*****************************************************************************
#  Copyright (C) 2021 The Sage Group -- http://www.sagemath.org/,
#      Fredrik Strömberg <stroemberg@mathematik.tu-darmstadt.de>,
#                David Berghaus <berghaus@th.physik.uni-bonn.de>,
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
import math

from sage.all import Gamma, is_odd
from sage.rings.all import QQ, ZZ, RR
from sage.modular.arithgroup.congroup_sl2z import SL2Z
from sage.modular.arithgroup.arithgroup_perm import sl2z_word_problem
from sage.modular.modsym.p1list import lift_to_sl2z
from sage.modular.cusps import Cusp

from psage.groups.permutation_alg cimport MyPermutation
from psage.modform.arithgroup.mysubgroups_alg cimport SL2Z_elt
from psage.modform.arithgroup.mysubgroups_alg import SL2Z_elt, closest_vertex, apply_sl2z_map_dp

cdef extern from "complex.h":
    cdef double cimag(double complex)
    cdef double creal(double complex)
    cdef double cabs(double complex)

cpdef pullback_to_gamma2_mat(double x, double y):
    """
    Returns matrix that maps point into fundamental domain of Gamma(2).
    """
    cdef SL2Z_elt T = SL2Z_elt(1,2,0,1)
    cdef SL2Z_elt T_inv = SL2Z_elt(1,-2,0,1)
    cdef SL2Z_elt T_0 = SL2Z_elt(1,0,-2,1)
    cdef SL2Z_elt T_0_inv = SL2Z_elt(1,0,2,1)
    cdef SL2Z_elt pb_map = SL2Z_elt(1,0,0,1)
    cdef double complex z_pb = x + y*1j
    while True:
        if abs(creal(z_pb)) > 1:
            if creal(z_pb) < 0: #We need to move point to the right
                z_pb += 2
                pb_map = T*pb_map
            else: #We need to move point to the left
                z_pb -= 2
                pb_map = T_inv*pb_map
        elif cabs(z_pb+0.5) < 0.5:
            z_pb = z_pb/(2*z_pb+1)
            pb_map = T_0_inv*pb_map
        elif cabs(z_pb-0.5) < 0.5:
            z_pb = z_pb/(-2*z_pb+1)
            pb_map = T_0*pb_map
        else:
            break
    return pb_map

cpdef pullback_to_gamma2_mat_perm(double x, double y, MyPermutation perm_inf, MyPermutation perm_inf_inv, MyPermutation perm_0, MyPermutation perm_0_inv):
    """
    Returns matrix that maps point into fundamental domain of Gamma(2)
    as well as the corresponding permutation.
    """
    cdef SL2Z_elt T = SL2Z_elt(1,2,0,1)
    cdef SL2Z_elt T_inv = SL2Z_elt(1,-2,0,1)
    cdef SL2Z_elt T_0 = SL2Z_elt(1,0,-2,1)
    cdef SL2Z_elt T_0_inv = SL2Z_elt(1,0,2,1)
    cdef SL2Z_elt pb_map = SL2Z_elt(1,0,0,1)
    cdef MyPermutation pb_perm = MyPermutation(length=perm_inf.N()) #Identity
    cdef double complex z_pb = x + y*1j
    while True:
        if abs(creal(z_pb)) > 1:
            if creal(z_pb) < 0: #We need to move point to the right
                z_pb += 2
                pb_map = T*pb_map
                pb_perm = perm_inf*pb_perm
            else: #We need to move point to the left
                z_pb -= 2
                pb_map = T_inv*pb_map
                pb_perm = perm_inf_inv*pb_perm
        elif cabs(z_pb+0.5) < 0.5:
            z_pb = z_pb/(2*z_pb+1)
            pb_map = T_0_inv*pb_map
            pb_perm = perm_0_inv*pb_perm
        elif cabs(z_pb-0.5) < 0.5:
            z_pb = z_pb/(-2*z_pb+1)
            pb_map = T_0*pb_map
            pb_perm = perm_0*pb_perm
        else:
            break
    return pb_map, pb_perm

cpdef pullback_general_gamma2_subgroup_dp(G,double x,double y,int ret_mat=0): 
        r"""
        Pullback for a general subgroup of Gamma(2).
        """
        cdef double xpb,ypb
        cdef int a,b,c,d,found_coset_rep,j
        cdef MyPermutation p,pj
        cdef SL2Z_elt B
        A, p = pullback_to_gamma2_mat_perm(x,y,G.perm_inf,G.perm_inf_inv,G.perm_0,G.perm_0_inv)
        # p = G.permutation_action(A) #.inverse()
        found_coset_rep=0
        for j in range(1,G._index+1):
            pj = G.coset_reps_perm[j]
            if p(pj(1))==1:
                found_coset_rep=1
                B=G.coset_reps[j]*A
                break
        if found_coset_rep==0:
            raise ArithmeticError,"Did not find coset rep. for A^-1={0}, x,y={1},{2}, G={3}".format(A.inverse(),x,y,G)
        xpb,ypb=apply_sl2z_map_dp(x,y,B[0],B[1],B[2],B[3])
        a,b,c,d = B[0],B[1],B[2],B[3]
        if ret_mat==1:
            return xpb,ypb,a,b,c,d
        else:
            return xpb,ypb

class Gamma_2_Subgroup:
    def __init__(self,MyPermutation o_0,MyPermutation o_inf):
        self._verbose = False #This functionality is currently not supported
        self.perm_0 = o_0
        self.perm_0_inv = o_0.inverse()
        self.perm_inf = o_inf
        self.perm_inf_inv = o_inf.inverse()
        self.perm_1_inv = o_inf*o_0
        self.perm_1 = (self.perm_1_inv).inverse()
        self._index = self.perm_inf.N()
        self._gamma_2_farey = Gamma(2).farey_symbol()
        self.coset_reps = self._get_coset_reps_from_perms()
        self.coset_reps_perm = self._init_coset_reps_perm()
        self._vertices, self._vertex_data, self._cusps, self._cusp_data = None, None, None, None
        self._nvertices, self._ncusps = None, None
        self._vertex_widths = list()
        self._vertex_maps, self._cusp_maps= list(), list()
        self._get_data_from_group()
        self._genus = None

    def __contains__(self, A):
        if A.SL2Z() not in self._gamma_2_farey:
            return False
        p = self.permutation_action(A)
        return p(1) == 1
    
    def _get_coset_reps_from_perms(self):
        p0 = self.perm_0; p1=self.perm_1; pI=self.perm_inf
        T_inf = SL2Z_elt(1,2,0,1)
        T_0 = SL2Z_elt(1,0,-2,1)
        Id = SL2Z_elt(1,0,0,1)
        coset_reps = dict()
        coset_reps[1] = Id
        cycI = pI.cycle_tuples()
        if len(cycI) > 1:
            raise NotImplementedError("This functionality has not been added yet!")
        else:
            iterations = 0
            for i in cycI[0]:
                if i != 1:
                    coset_reps[i] = T_inf**iterations
                iterations += 1
        return coset_reps

    def _init_coset_reps_perm(self):
        coset_reps = self.coset_reps
        coset_reps_perm = dict()
        for i in range(1,len(coset_reps)+1):
            coset_reps_perm[i] = self.permutation_action(coset_reps[i])
        return coset_reps_perm

    def _get_data_from_group(self):  
        ## Get information about cusps and vertices
        l = self._get_all_cusp_data()
        self._vertices, self._vertex_data, self._cusps, self._cusp_data = l
        self._nvertices, self._ncusps = len(self._vertices), len(self._cusps)
        for i in range(len(self._vertices)):
            wi = self._cusp_data[self._vertex_data[i]['cusp']]['width']
            self._vertex_widths.append(wi)
            N = self._cusp_data[self._vertex_data[i]['cusp']]['normalizer']
            N = SL2Z_elt(N[0],N[1],N[2],N[3])
            U = self._vertex_data[i]['cusp_map']
            self._cusp_maps.append(U) #[U[0,0],U[0,1],U[1,0],U[1,1]])
            N = N.inverse()*U
            self._vertex_maps.append(N) #[N[0,0],N[0,1],N[1,0],N[1,1]])

    def _get_all_cusp_data(self):
        coset_reps = self.coset_reps
        Id = SL2Z_elt(1,0,0,1)
        cusps = list()
        vertices = list()
        vertex_data = dict()
        cusp_data = dict()
        #We begin with the vertex at infinity which we also choose to be a cusp representative
        vertex_data[0] = {'cusp':0,'cusp_map':Id,'coset':self.perm_inf.cycles_as_lists()[0]}
        vertices.append(Cusp(1,0))
        cusp_data[0] = {'width':self.cusp_width(Cusp(1,0)),'normalizer':Id}
        cusps.append(Cusp(1,0))
        vi = 1
        for i in range(1,len(coset_reps)+1):
            if coset_reps[i][2]==0:
                center = 2*(i-1) #Center of coset_rep
                for j in (-1,0,1): #These are the positions on the real line of the vertices of the standard fundamental domain
                    v = Cusp(center+j,1)
                    if v not in vertices: #Found a new vertex
                        is_cusp = True
                        for c in cusps:
                            cusp_map = self.are_equivalent(v, c, trans=True)
                            if cusp_map != False: #Found an equivalent cusp
                                ci = cusps.index(c)
                                vertex_data[vi] = {'cusp':ci,'cusp_map':cusp_map,'coset':[i]}
                                vertices.append(v)
                                is_cusp = False
                                break
                        if is_cusp: #Found a new cusp representative
                            ci = len(cusps)
                            vertex_data[vi] = {'cusp':ci,'cusp_map':Id,'coset':[i]}
                            vertices.append(v)
                            cusp_data[ci] = {'width':self.cusp_width(v),'normalizer':self.cusp_normalizer(v)}
                            cusps.append(v)
                        vi += 1
                    else:
                        vj = vertices.index(v)
                        vertex_data[vj]['coset'].append(i)
            else:
                raise NameError("This case has not been implemented yet")
        return vertices, vertex_data, cusps, cusp_data

    def permutation_action(self, A):
        """
        Returns permutation corresponding to A.
        """
        w = self._gamma_2_word_problem(A)
        p = MyPermutation(length=self._index) #identity
        for i in w:
            if i[0] == 0:
                p = p*self.perm_inf**i[1] #Note that we multiply from the right
            else:
                p = p*self.perm_1_inv**i[1] #Note that we multiply from the right
        return p

    def _gamma_2_word_problem(self, A):
        """
        Solves word problem of SL2Z_elt A in Gamma(2).
        We work with the implementation of the word problem of Gamma(2) in sage
        which uses T_inf and T_1_inv as first and second generators respectively.
        This function is quite slow and should not be used extensively.
        """
        gamma2_w = self._gamma_2_farey.word_problem(A.SL2Z(), output='syllables')
        w = list()
        for i in gamma2_w:
            if i[0] != 2: #We ignore the third generator which is [-1,0,0,-1] since it gets projected out
                w.append(i)
        return w

    def are_equivalent(self,x,y,trans=False):
        r"""
        Check whether two cusps are equivalent with respect to self
        The algorithm is the same as in sage except that we use SL2Z_elt to increase speed
        and modified the generator to the one of Gamma(2).
        """
        x = Cusp(x)
        y = Cusp(y)
        vx = lift_to_sl2z(x.numerator(),x.denominator(), 0)
        dx = SL2Z_elt(vx[2], -vx[0], vx[3], -vx[1])
        vy = lift_to_sl2z(y.numerator(),y.denominator(), 0)
        dy = SL2Z_elt(vy[2], -vy[0], vy[3], -vy[1])
        for i in range(self._index):
            # Note that the width of any cusp is bounded above by the index of self.
            t = dy * SL2Z_elt(1,2*i,0,1)*dx.inverse() #Note that T_inf is different to SL2Z here
            if t in self:
                if trans:
                    return t
                else:
                    return True
        return False

    def cusp_width(self,c):
        r""":
        Returns width of cusp c. 
        INPUT:
        - 'c' -- Integer or cusp
        """
        cusp = Cusp(c)
        if self._cusps != None and cusp in self._cusps:
            j = self._cusps.index(cusp)
            return self._cusp_data[j]['width']
        else:
            g = self.cusp_normalizer(c)
            g_inv = g.inverse()
            for d in range(1,1+self._index):
                t = g * SL2Z_elt(1,2*d,0,1) * g_inv #Note that T_inf is different to SL2Z here
                if t in self:
                    return 2*d
            raise ArithmeticError("Can't get here!")
    
    def cusp_normalizer(self,c):
        r"""
        Return the cusp normalizer of cusp c.
        """
        cusp = Cusp(c)
        if self._cusps != None and cusp in self._cusps:
            j = self._cusps.index(cusp)
            return self._cusp_data[j]['normalizer']
        else:
            w = lift_to_sl2z(cusp.denominator(), cusp.numerator(), 0)
            g = SL2Z_elt(w[3], w[1], w[2], w[0])
            return g

    def closest_vertex(self,x,y,as_integers=1):
        r"""
        The closest vertex to the point z=x+iy in the following sense:
        Let sigma_j be the normalized cusp normalizer of the vertex p_j, 
        i.e. sigma_j^-1(p_j)=Infinity and sigma_j*S_j*sigma_j^-1=T, where
        S_j is the generator of the stabiliser of p_j
        The closest vertex is then the one for which Im(sigma_j^-1(p_j))
        is maximal.
        INPUT:
         - ''x,y'' -- x+iy  in the upper half-plane
        OUTPUT:
        
         - ''v'' -- the closest vertex to x+iy
         
        
        EXAMPLES::
        sage: G=MySubgroup(Gamma0(5))
        sage: G.closest_vertex(-0.4,0.2)
        Infinity
        sage: G.closest_vertex(-0.1,0.1)
        0
        """
        ci=closest_vertex(self._vertex_maps,self._vertex_widths,self._nvertices,x,y,self._verbose)
        if as_integers:
            return ci
        else:
            return self._vertices[ci]

    def closest_cusp(self,x,y,vertex=0,as_integers=1):
        r"""
        The closest cusp to the point z=x+iy in the following sense:
        Let sigma_j be the normalized cusp normalizer of the vertex p_j, 
        i.e. sigma_j^-1(p_j)=Infinity and sigma_j*S_j*sigma_j^-1=T, where
        S_j is the generator of the stabiliser of p_j
        The closest vertex is then the one for which Im(sigma_j^-1(p_j))
        is maximal and the closest cusp is the cusp associated to this vertex
        INPUT:
         - ''x,y'' -- x+iy  in the upper half-plane
        OUTPUT: 
         - ''v'' -- the closest vertex to x+iy
        """
        vi = closest_vertex(self._vertex_maps,self._vertex_widths,self._nvertices,x,y,self._verbose)
        ci = self._vertex_data[vi]['cusp']
        if vertex==1:
            if as_integers:
                return ci, vi
            else:
                return self._cusps[ci], self._vertices[vi]
        else:
            if as_integers:
                return ci
            else:
                return self._cusps[ci]

    def minimal_height(self):
        r""" Computes the minimal (invariant) height of the fundamental region of self.
        Note: Not guaranteed.
        """
        # Locate the largest width
        maxw=0
        for i in range(self._ncusps):
            l=self._cusp_data[i]['width']
            if l>maxw:
                maxw=l
        return math.sqrt(3)/(2*maxw)
    
    def ncusps(self):
        return self._ncusps
    
    def genus(self):
        if self._genus == None:
            #Begin by computing the excesses of the generating permutations
            tmp = self.perm_0.cycle_lens()
            e_perm_0 = sum(tmp)-len(tmp)
            tmp = self.perm_1.cycle_lens()
            e_perm_1 = sum(tmp)-len(tmp)
            tmp = self.perm_inf.cycle_lens()
            e_perm_inf = sum(tmp)-len(tmp)
            self._genus = 1 - self._index + QQ(e_perm_0+e_perm_1+e_perm_inf)/QQ(2)
        return self._genus

    def dimension_cusp_forms(self,k):
        r"""
        Returns the dimension of the space of cuspforms on G of weight k
        where k is an even integer
        """
        ki = ZZ(k)
        if is_odd(ki):
            raise ValueError("Use only for even weight k! not k={0}".format(k))
        if ki < 2:
            dim=0 
        elif ki == 2:
            dim=self.genus()
        elif ki >= 4:
            kk = RR(k)
            dim=ZZ(kk-1)*(self._genus-1) #Note that Gamma(2) has no elliptic points
            dim+= ZZ(kk/2.0 - 1)*self._ncusps
        return dim
