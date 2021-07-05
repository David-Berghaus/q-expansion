#include "acb.h"
#include "acb_mat.h"

void acb_dft_column_precomp(acb_mat_t b, acb_dft_pre_t pre_comp, acb_mat_t x, int prec);

void acb_compute_dft_matrix(acb_mat_t A, int N, int bit_prec);

void acb_test_fft(int N, int bit_prec);