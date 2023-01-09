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

"""
NOTICE: This code is part of psage https://github.com/fredstro/psage/blob/master/psage and has been copied with permission of the license holder.
We copy it here to remove the dependency on the installation of psage which can be a bit tricky sometimes.
"""

from psage.modform.arithgroup.all import MySubgroup,MySubgroup_class
from sage.all import Gamma0
from sage.rings.integer import Integer
from sage.rings.rational import Rational

class AutomorphicFormSpace():
    def __init__(self,G,weight=0):
        r""" Initialize the space of automorphic forms.
        Note that this class is much smaller than the one provided by psage because we do not use all its functionalities.
        """
        if isinstance(G,MySubgroup_class):
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
        self._weight = weight
    
    def group(self):
        return self._group
    
    def weight(self):
        return self._weight

def is_int(q):
    r"""
    Find out if the rational number q is an integer.
    INPUT:
    -''q'' -- integer/rational/real
    OUTPUT:
    - logical -- True if q is an integer otherwise False
    EXAMPLES::
        sage: is_int(1)
        True
        sage: is_int(float(1.0))
        True
        sage: is_int(RR(1.0))   
        True
        sage: is_int(6/3) 
        True
        sage: is_int(6/4)
        False
        sage: is_int(1.5)    
        False
        sage: is_int(Gamma0(1))    
        False
        
    """
    if(isinstance(q,Integer) or isinstance(q,int)):
        return True
    if(isinstance(q,Rational)):
        n=q.denominator()
        if(n==1):
            return True
    if(isinstance(q,tuple)):
        return False
    try:
        if(floor(q)==ceil(q)):
            return True
    except:
        pass
    return False