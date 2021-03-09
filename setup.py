#sage setup.py develop
#python setup.py build_ext -> for getting html file

from setuptools import setup
from Cython.Build import cythonize
import Cython.Compiler.Options
Cython.Compiler.Options.annotate=True

setup(
    ext_modules=cythonize(["my_pullback.pyx","point_matching.pyx"],annotate=True),
    zip_safe=False,
)