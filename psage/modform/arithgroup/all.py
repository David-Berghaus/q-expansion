"""
NOTICE: This code is part of psage https://github.com/fredstro/psage/blob/master/psage and has been copied with permission of the license holder.
We copy it here to remove the dependency on the installation of psage which can be a bit tricky sometimes.
"""

from __future__ import absolute_import
from . import mysubgroups_alg

from .mysubgroups_alg import SL2Z_elt,factor_matrix_in_sl2z,ncf_to_SL2Z_element,GL2Z_elt
from .mysubgroup import MySubgroup,HeckeTriangleGroup,nearest_integer_continued_fraction,list_valid_signatures,MySubgroup_class