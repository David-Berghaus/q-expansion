from sage.libs.arb.arb cimport *
from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_mat cimport *

from arblib_helpers.acb_approx cimport *

def test_fft(int N, bit_prec):
    acb_test_fft(N, bit_prec)
