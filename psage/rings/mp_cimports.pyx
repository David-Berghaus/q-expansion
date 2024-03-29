r"""
Standard cimports for working with multiprecision numbers.

NOTICE: This code is part of psage https://github.com/fredstro/psage/blob/master/psage and has been copied with permission of the license holder.
We copy it here to remove the dependency on the installation of psage which can be a bit tricky sometimes.
"""
from cysignals.memory cimport sig_free,sig_malloc
from cysignals.signals cimport sig_on,sig_off
from sage.ext.stdsage cimport PY_NEW
from sage.libs.gmp.all cimport *
from sage.libs.mpfr cimport *

from sage.rings.complex_mpc cimport * 
from sage.rings.complex_mpc cimport MPComplexNumber
from sage.rings.real_mpfr cimport RealNumber,RealField_class
from sage.rings.integer cimport Integer
from sage.rings.complex_mpc cimport MPComplexNumber
from sage.rings.complex_mpc cimport MPComplexField_class
from sage.rings.real_mpfr cimport RealNumber

## Define the rounding modes
cdef mpc_rnd_t rnd
cdef mpfr_rnd_t rnd_re
rnd = MPC_RNDNN
rnd_re = MPFR_RNDN