#sage setup.py develop
#python setup.py build_ext -> for getting html file

from setuptools import setup, Extension
from Cython.Build import cythonize
import Cython.Compiler.Options

from sage.env import cython_aliases

Cython.Compiler.Options.annotate = False #Set to "True" if html should be created that highlights python parts

extensions = [
    Extension("classes.acb_mat_class", ["classes/acb_mat_class.pyx"]),
    Extension("classes.acb_dft_class", ["classes/acb_dft_class.pyx"]),
    Extension("classes.plu_class", ["classes/plu_class.pyx"]),
    Extension("classes.block_factored_mat_class", ["classes/block_factored_mat_class.pyx"]),
    Extension("classes.modform_class", ["classes/modform_class.pyx"]),
    Extension("classes.gamma_2_subgroup", ["classes/gamma_2_subgroup.pyx"]),
    Extension("classes.factored_polynomial", ["classes/factored_polynomial.pyx"]),
    Extension("classes.fourier_expansion", ["classes/fourier_expansion.pyx"]),
    Extension("classes.belyi_map", ["classes/belyi_map.pyx"]),
    Extension("pullback.my_pullback", ["pullback/my_pullback.pyx"]),
    Extension("point_matching.point_matching_dp", ["point_matching/point_matching_dp.pyx"]),
    Extension("point_matching.point_matching_arb_wrap", ["point_matching/point_matching_arb_wrap.pyx"]),
    Extension("point_matching.fft_test", ["point_matching/fft_test.pyx"]),
    Extension("iterative_solvers.gmres_dp", ["iterative_solvers/gmres_dp.pyx"]),
    Extension("iterative_solvers.gmres_arb_wrap", ["iterative_solvers/gmres_arb_wrap.pyx"]),
    Extension("iterative_solvers.iterative_refinement_arb_wrap", ["iterative_solvers/iterative_refinement_arb_wrap.pyx"]),
    Extension("belyi.newton_genus_zero", ["belyi/newton_genus_zero.pyx"]),
    Extension("belyi.number_fields", ["belyi/number_fields.pyx"]),
    Extension("belyi.expression_in_u_and_v", ["belyi/expression_in_u_and_v.pyx"]),
    Extension("eisenstein.haberland", ["eisenstein/haberland.pyx"]),
    Extension("eisenstein.eisenstein_computation", ["eisenstein/eisenstein_computation.pyx"]),
]

setup(
    ext_modules = cythonize(extensions,aliases=cython_aliases()),
    zip_safe = False,
)
