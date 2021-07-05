from cysignals.signals cimport sig_on, sig_str, sig_off

from arblib_helpers.acb_approx cimport *
from classes.acb_mat_class cimport Acb_Mat, Acb_Mat_Win

cdef class Acb_DFT():
    def __cinit__(self, int two_Q, int bit_prec):
        sig_str("Arb exception")
        acb_dft_precomp_init(self.value, two_Q, bit_prec)
        sig_off()

    def __dealloc__(self):
        acb_dft_precomp_clear(self.value)

    def perform_fft(self, Acb_Mat b, Acb_Mat x, int bit_prec):
        sig_on()
        acb_dft_column_precomp(b.value, self.value, x.value, bit_prec)
        sig_off() 