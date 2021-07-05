from sage.libs.arb.acb_mat cimport *

cdef _get_J_block_matrix_arb_wrap(acb_mat_t J,int Ms,int Mf,int weight,int Q,coordinates,int bit_prec)

cdef _get_W_block_matrix_arb_wrap(acb_mat_t W,int Ms,int Mf,int weight,coordinates,int bit_prec)