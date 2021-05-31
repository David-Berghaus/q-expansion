void
acb_mat_approx_add(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, int prec);

void
acb_mat_approx_sub(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, int prec);

void
acb_mat_approx_scalar_addmul(acb_mat_t res, acb_mat_t B, const acb_mat_t A, const acb_t c, int prec);

void //Computes c*A where c is a complex scalar and A is a matrix
acb_mat_approx_scalar_mul(acb_mat_t res, const acb_mat_t A, const acb_t c, int prec);

void //Computes c*A where c is a real scalar and A is a matrix
acb_mat_approx_scalar_mul_arb(acb_mat_t res, const acb_mat_t A, const arb_t c, int prec);

void //Computes conjugate(x) dot y where x&y are len*1 matrices.
acb_mat_approx_dotc(acb_t res, acb_mat_t x, acb_mat_t y, int prec);

void //Computes x dot y where x&y are len*1 matrices.
acb_mat_approx_dot(acb_t res, acb_mat_t x, acb_mat_t y, int prec)

void acb_mat_approx_norm(arb_t res, acb_mat_t x, int prec);

void acb_approx_complex_sign(acb_t res, acb_t z, arb_t z_abs, int prec);

void //See Algorithm 1 of https://www.netlib.org/lapack/lawnspdf/lawn148.pdf
lartg(acb_t c, acb_t s, acb_t r, acb_t f, acb_t g, int prec);

void acb_mat_change_prec(acb_mat_t res, acb_mat_t A, int prec);

void //Computes D*A where D is a diagonal matrix which is stored as a Nx1 vector
acb_mat_approx_left_mul_diag(acb_mat_t res, const acb_mat_t D, const acb_mat_t A, int prec);

void //Computes A*D where D is a diagonal matrix which is stored as a Nx1 vector
acb_mat_approx_right_mul_diag(acb_mat_t res, const acb_mat_t A, const acb_mat_t D, int prec);