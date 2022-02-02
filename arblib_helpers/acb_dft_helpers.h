#include "acb.h"
#include "acb_mat.h"

//This is an ugly include of acb_dft, it would be cleaner if Sage would include this in its arb wrapper...
#define ulong mp_limb_t
#define slong mp_limb_signed_t
#include "acb_dft.h"
#undef ulong
#undef slong

void acb_dft_column_precomp(acb_mat_t b, acb_dft_pre_t pre_comp, acb_mat_t x, int prec);

void acb_compute_dft_matrix(acb_mat_t A, int N, int bit_prec);

void acb_test_fft(int N, int bit_prec);
