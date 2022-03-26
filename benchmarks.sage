"""
Some benchmarks to compute a single cuspform using different approaches
"""
import resource

from point_matching.point_matching_arb_wrap import get_coefficients_cuspform_arb_wrap, get_coefficients_cuspform_ir_arb_wrap, digits_to_bits, get_horo_height_arb_wrap, get_V_tilde_matrix_factored_b_cuspform_arb_wrap, get_M_0, get_Q
from iterative_solvers.gmres_arb_wrap import gmres_mgs_arb_wrap

#resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
def run_benchmark(S,digit_prec,version):
    """
    Run benchmark with a version of choice.
    """
    import time
    t = time.time()
    if version == 1:
        res = classical_benchmark(S,digit_prec)
    if version == 2:
        res = non_precond_gmres_benchmark(S,digit_prec)
    if version == 3:
        res = optimized_ir_benchmark(S,digit_prec)
    print("Elapsed time in seconds: ", time.time()-t)
    max_mem_usage_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1e6) #CAREFUL, ONLY RUN ONE BENCHMARK PER PROCESS!!!
    print("Max memory usage in Gb: ", max_mem_usage_gb)

def classical_benchmark(S,digit_prec):
    """
    Compute cuspform classically. We first construct J and W and afterwards form V_tilde through matrix multiplication.
    Then we use arb's direct solving routine to compute the coefficients.
    """
    bit_prec = digits_to_bits(digit_prec)
    res = get_coefficients_cuspform_arb_wrap(S,digit_prec)._get_mcbd(bit_prec)
    print("res[:5]", res[:5].str())
    return res

def non_precond_gmres_benchmark(S,digit_prec):
    """
    Compute cuspform using GMRES without preconditioner. We form V_tilde only in factored form and use GMRES for V_tilde_factored_sc
    without preconditioner. The actions of J and W are computed using classical matrix-vector multiplication.
    """
    bit_prec = digits_to_bits(digit_prec)
    RBF = RealBallField(bit_prec)
    CBF = ComplexBallField(bit_prec)
    M_0 = get_M_0(S,digit_prec)
    Y = get_horo_height_arb_wrap(S,RBF,M_0)
    Q = get_Q(Y,S.weight(),digit_prec,is_cuspform=True)
    print("Y = ", Y)
    print("M_0 = ", M_0)
    print("Q = ", Q)
    print("ncusps = ", S.group().ncusps())
    V, b_vecs = get_V_tilde_matrix_factored_b_cuspform_arb_wrap(S,M_0,Q,Y,bit_prec,use_FFT=False,use_splitting=False,labels=[0])
    b = b_vecs[0]

    tol = RBF(10.0)**(-digit_prec)
    x_gmres_arb_wrap = gmres_mgs_arb_wrap(V, b, bit_prec, tol, maxiter=10000)

    res = x_gmres_arb_wrap[0]
    V.diag_inv_scale_vec(res, res, bit_prec)
   
    res = x_gmres_arb_wrap[0]._get_mcbd(bit_prec)
    print("res[:5]", res[:5].str())
    return res

def optimized_ir_benchmark(S,digit_prec):
    bit_prec = digits_to_bits(digit_prec)
    res = get_coefficients_cuspform_ir_arb_wrap(S,digit_prec,use_FFT=True,use_splitting=True,use_scipy_lu=True)._get_mcbd(bit_prec)
    print("res[:5]", res[:5].str())
    return res