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

from sage.all import Gamma
from sage.modular.arithgroup.arithgroup_perm import sl2z_word_problem
from sage.modular.modsym.p1list import lift_to_sl2z
from sage.modular.cusps import Cusp

from psage.groups.permutation_alg cimport MyPermutation
from psage.modform.arithgroup.mysubgroups_alg import SL2Z_elt

class Gamma_2_Subgroup:
    def __init__(self,MyPermutation o_0,MyPermutation o_inf):
        self.perm_0 = o_0
        self.perm_inf = o_inf
        self.perm_1_inv = o_inf*o_0 #This is one of the gens of Gamma(2).farey_symbol() which is why we need it
        self.perm_1 = (self.perm_1_inv).inverse()
        self._index = self.perm_inf.N()
        self._gamma_2_farey = Gamma(2).farey_symbol()

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

    def permutation_action(self, A):
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
    
    def cusp_data(self,c):
        r""":
        Returns cuspdata in the same format as for the generic Arithmetic subgroup, i.e. a tuple (A,h,s) where A is a generator of the stabiliser of c, h is the width of c and s is the orientation. 
        INPUT:
        - 'c' -- Integer or cusp
        """       
        cusp = Cusp(c)
        ## Then we compute everything using the same method as in sage but with SL2Z_elt
        ## to make it faster
        w = lift_to_sl2z(c.denominator(), c.numerator(), 0)
        g = SL2Z_elt([w[3], w[1], w[2],w[0]])
        for d in range(1,1+self._index):
            t = g * SL2Z_elt(1,2*d,0,1) * g.inverse()
            if t in self:
                return t, 2*d, 1
        raise ArithmeticError("Can't get here!")
