from cysignals.signals cimport sig_on, sig_off

from math import sqrt

from sage.rings.complex_arb cimport *
from sage.libs.arb.acb_mat cimport *
from sage.rings.real_arb import RealBallField
from sage.rings.complex_arb import ComplexBallField

from arblib_helpers.acb_approx cimport *

cdef class Modular_splitting_polynomial():
    """
    Class for storing and evaluating rows of W using modular splitting.
    """
    def __init__(self, int Ms, int Mf, ComplexBall q, RealBall y_fund_fact, int bit_prec):
        cdef CBF = ComplexBallField(bit_prec)
        cdef ComplexBall tmp = CBF(0)
        self.y_fund_fact = y_fund_fact
        coeff_len = Mf+1
        n_sqrt = int(sqrt(coeff_len))
        k = n_sqrt
        if k**2 != coeff_len: #coeff_len is not a square number
            k += 1
        j = n_sqrt
        self.xs, self.ys = Acb_Mat(1,j), Acb_Mat(1,k)
        acb_mat_set_powers_approx(self.xs.value,q.value,bit_prec)
        acb_approx_mul(tmp.value,acb_mat_entry(self.xs.value,0,j-1),q.value,bit_prec) #This is y = x^j
        acb_mat_set_powers_approx(self.ys.value,tmp.value,bit_prec)
        self.j, self.k, self.Ms, self.Mf = j, k, Ms, Mf

    cdef evaluate(self, acb_t res, Acb_Mat coeffs, int bit_prec):
        sig_on()
        evaluate_modular_splitting_polynomial(res,coeffs.value,self.xs.value,self.ys.value,self.j,self.k,self.Ms,self.Mf,bit_prec)
        acb_approx_mul_arb(res,res,self.y_fund_fact.value,bit_prec)
        sig_off()

def test_evaluate(z, two_pi_i, Ms, Mf):
    import time
    CBF = z.parent()
    cdef ComplexBall q = (two_pi_i*z).exp()
    print("q: ", q)
    cdef Acb_Mat coeffs = Acb_Mat(Mf+1,1)
    for i in range(Mf+1):
        #acb_one(acb_mat_entry(coeffs.value,i,0))
        acb_set_ui(acb_mat_entry(coeffs.value,i,0),i-Ms+1)
    bit_prec = CBF.precision()
    RBF = RealBallField(bit_prec)
    cdef RealBall y_fund_fact = RBF(1)/2
    t = time.time()
    p_split = Modular_splitting_polynomial(Ms, Mf, q, y_fund_fact, bit_prec)
    print("Modular splitting setup: ", time.time()-t)
    cdef ComplexBall res = CBF(0)
    t = time.time()
    p_split.evaluate(res.value,coeffs,bit_prec)
    print("Modular splitting evaluate: ", time.time()-t)
    print("res (modular splitting): ", res)
    cdef Acb_Mat q_arr = Acb_Mat(1,Mf+1)
    cdef Acb_Mat res_arr = Acb_Mat(1,1)

    # t = time.time()
    # acb_one(acb_mat_entry(q_arr.value,0,0))
    # for i in range(1,Mf+1):
    #     acb_approx_mul(acb_mat_entry(q_arr.value,0,i),acb_mat_entry(q_arr.value,0,i-1),q.value,bit_prec)
    # print("q_arr setup: ", time.time()-t)
    # t = time.time()
    # acb_mat_approx_mul(res_arr.value,q_arr.value,coeffs.value,bit_prec)
    # acb_approx_mul_arb(acb_mat_entry(res_arr.value,0,0),acb_mat_entry(res_arr.value,0,0),y_fund_fact.value,bit_prec)
    # print("q_arr evaluate: ", time.time()-t)
    # print("res (q_arr): ", res_arr._get_mcbd(bit_prec)) #We are obviously not considering Ms here
    
    cdef ComplexBall tmp = y_fund_fact*((two_pi_i*Ms*z).exp())
    acb_set(acb_mat_entry(q_arr.value, 0, 0), tmp.value) #q_arr actually contains y_fund_fact here
    for l in range(Ms+1,Mf+1):
        acb_approx_mul(acb_mat_entry(q_arr.value, 0, l-Ms), acb_mat_entry(q_arr.value, 0, l-Ms-1), q.value, bit_prec)
    acb_mat_approx_mul(res_arr.value,q_arr.value,coeffs.value,bit_prec)
    print("q_arr evaluate: ", time.time()-t)
    print("res (q_arr): ", res_arr._get_mcbd(bit_prec)) #We are obviously not considering Ms here