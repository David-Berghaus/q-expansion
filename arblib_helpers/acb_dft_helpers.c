#include "acb.h"
#include "acb_mat.h"
#include "acb_approx_helpers.h"

void //Perform DFT on matrix column
acb_dft_column_precomp(acb_mat_t b, acb_dft_pre_t pre_comp, acb_mat_t x, int prec)
{
    int len = acb_mat_nrows(x);
    acb_ptr tmp_vec;

    TMP_INIT;

    TMP_START;
    tmp_vec = TMP_ALLOC(sizeof(acb_struct) * len);

    int i;
    for (i = 0; i < len; i++)
    {
        tmp_vec[i] = *acb_mat_entry(x, i, 0);
    }

    acb_dft_precomp(tmp_vec, tmp_vec, pre_comp, prec);
    for (i = 0; i < len; i++)
    {
        acb_swap(acb_mat_entry(b, i, 0), &tmp_vec[i]);
    }

    TMP_END;
}

void acb_compute_dft_matrix(acb_mat_t A, int N, int bit_prec){
    int i, j;
    acb_t two_pi_i_over_N, acb_tmp;
    acb_init(acb_tmp);
    acb_init(two_pi_i_over_N);
    acb_const_pi(two_pi_i_over_N, bit_prec);
    acb_mul_si(two_pi_i_over_N, two_pi_i_over_N, 2, bit_prec);
    acb_div_si(two_pi_i_over_N, two_pi_i_over_N, N, bit_prec);
    arb_swap(acb_realref(two_pi_i_over_N), acb_imagref(two_pi_i_over_N));

    acb_mat_ones(A);
    for (i = 1; i < N; i++){
        acb_mul_si(acb_tmp, two_pi_i_over_N, i, bit_prec);
        acb_exp(acb_tmp, acb_tmp, bit_prec);
        for (j = 1; j < N; j++){
            acb_approx_mul(acb_mat_entry(A,i,j), acb_mat_entry(A,i,j-1), acb_tmp, bit_prec);
        }
    }
    acb_clear(acb_tmp);
    acb_clear(two_pi_i_over_N);
}

void acb_test_fft(int N, int bit_prec){
    int i;
    acb_mat_t F, dft_res, x;
    acb_ptr x_vec = _acb_vec_init(N);
    acb_ptr fft_res_vec = _acb_vec_init(N);
    acb_mat_init(x, N, 1);
    acb_t acb_tmp;
    arb_t arb_tmp;
    acb_init(acb_tmp);
    arb_init(arb_tmp);
    arb_const_pi(acb_realref(acb_tmp), bit_prec);
    arb_const_e(acb_imagref(acb_tmp), bit_prec);
    for (i = 0; i < N; i++){
        // arb_set_ui(arb_tmp, i+1);
        // arb_sqrt(arb_tmp, arb_tmp, bit_prec);
        // acb_approx_mul_arb(acb_tmp, acb_tmp, arb_tmp, bit_prec); //Just create some (seemingly) random entries
        acb_set(&x_vec[i], acb_tmp);
        acb_set(acb_mat_entry(x,i,0), acb_tmp);
    }
    acb_mat_init(F, N, N);
    acb_mat_init(dft_res, N, 1);
    acb_dft_pre_t pre_comp;
    clock_t start, end;
    double elapsed_time;

    // start = clock();
    // acb_compute_dft_matrix(F, N, bit_prec);
    // end = clock();
    // elapsed_time = (double)(end - start)/(double)CLOCKS_PER_SEC;
    // printf("DFT Matrix computation took: %f.\n", elapsed_time);
    // start = clock();
    // acb_mat_approx_mul(dft_res, F, x, bit_prec);
    // end = clock();
    // elapsed_time = (double)(end - start)/(double)CLOCKS_PER_SEC;
    // printf("Classic DFT took: %f.\n", elapsed_time);
    // acb_mat_printd(dft_res, 10);
    // printf("\n");

    start = clock();
    acb_dft_precomp_init(pre_comp, N, bit_prec);
    end = clock();
    elapsed_time = (double)(end - start)/(double)CLOCKS_PER_SEC;
    printf("FFT pre-computation took: %f.\n", elapsed_time);
    start = clock();
    acb_dft_precomp(fft_res_vec, x_vec, pre_comp, bit_prec);
    end = clock();
    elapsed_time = (double)(end - start)/(double)CLOCKS_PER_SEC;
    printf("FFT took: %f.\n", elapsed_time);
    acb_printd(&fft_res_vec[0], 10);
    // for (i = 0; i < N; i++){
    //     acb_printd(&fft_res_vec[i], 10);
    //     printf("\n");
    // }

    acb_mat_clear(F);
    acb_mat_clear(dft_res);
    acb_mat_clear(x);
    _acb_vec_clear(x_vec, N);
    _acb_vec_clear(fft_res_vec, N);
    acb_dft_precomp_clear(pre_comp);
    acb_clear(acb_tmp);
    arb_clear(arb_tmp);
}