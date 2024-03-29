#sage setup.py develop
#python setup.py build_ext -> for getting html file

from setuptools import setup, Extension
from Cython.Build import cythonize
import Cython.Compiler.Options
import numpy

from sage.env import cython_aliases

Cython.Compiler.Options.annotate = False #Set to "True" if html should be created that highlights python parts

extensions = [
    Extension("classes.acb_mat_class", ["classes/acb_mat_class.pyx"]),
    Extension("classes.acb_dft_class", ["classes/acb_dft_class.pyx"]),
    Extension("classes.plu_class", ["classes/plu_class.pyx"]),
    Extension("classes.block_factored_mat_class", ["classes/block_factored_mat_class.pyx"]),
    Extension("classes.modular_splitting_polynomial", ["classes/modular_splitting_polynomial.pyx"]),
    Extension("classes.modform_class", ["classes/modform_class.pyx"]),
    Extension("classes.gamma_2_subgroup", ["classes/gamma_2_subgroup.pyx"]),
    Extension("classes.factored_polynomial", ["classes/factored_polynomial.pyx"]),
    Extension("classes.fourier_expansion", ["classes/fourier_expansion.pyx"]),
    Extension("classes.belyi_map", ["classes/belyi_map.pyx"]),
    Extension("pullback.my_pullback", ["pullback/my_pullback.pyx"]),
    Extension("point_matching.point_matching_dp", ["point_matching/point_matching_dp.pyx"]),
    Extension("point_matching.point_matching_arb_wrap", ["point_matching/point_matching_arb_wrap.pyx"]),
    Extension("iterative_solvers.gmres_dp", ["iterative_solvers/gmres_dp.pyx"]),
    Extension("iterative_solvers.gmres_arb_wrap", ["iterative_solvers/gmres_arb_wrap.pyx"]),
    Extension("iterative_solvers.iterative_refinement_arb_wrap", ["iterative_solvers/iterative_refinement_arb_wrap.pyx"]),
    Extension("belyi.newton_genus_zero", ["belyi/newton_genus_zero.pyx"]),
    Extension("belyi.number_fields", ["belyi/number_fields.pyx"]),
    Extension("belyi.expression_in_u_and_v", ["belyi/expression_in_u_and_v.pyx"]),
    Extension("belyi.elliptic_curve", ["belyi/elliptic_curve.pyx"]),
    Extension("eisenstein.haberland", ["eisenstein/haberland.pyx"]),
    Extension("eisenstein.eisenstein_computation", ["eisenstein/eisenstein_computation.pyx"]),
    Extension("psage.groups.permutation_alg", ["psage/groups/permutation_alg.pyx"]),
    Extension("psage.modform.arithgroup.mysubgroup", ["psage/modform/arithgroup/mysubgroup.pyx"]),
    Extension("psage.modform.arithgroup.mysubgroups_alg", ["psage/modform/arithgroup/mysubgroups_alg.pyx"]),
    Extension("psage.modform.maass.automorphic_forms_alg", ["psage/modform/maass/automorphic_forms_alg.pyx"]),
    Extension("psage.rings.mp_cimports", ["psage/rings/mp_cimports.pyx"]),
]

setup(
    ext_modules = cythonize(extensions,aliases=cython_aliases()),
    include_dirs=[numpy.get_include()],
    zip_safe = False,
)
