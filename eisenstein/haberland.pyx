import functools

from sage.functions.other import binomial
from sage.misc.misc_c import prod
from sage.all import cached_function
from sage.rings.rational_field import QQ
from sage.symbolic.constants import pi

@cached_function
def memoized_binomial(k, n):
    """
    Returns the binomial coefficient and caches the result.
    """
    return binomial(k,n)

@cached_function
def memoized_prod(a, b):
    """
    Returns the product a*...*(b-1) and caches the result.
    """
    return prod(range(a,b)) #We could probably also try to improve this by computing it recursively but have not done this yet

@cached_function
def memoized_pow(a, n):
    """
    Compute a**n recursively to cache intermediate results.
    """
    if n == 0:
        return 1
    if n == 1:
        return a
    return a*memoized_pow(a,n-1)

@cached_function
def memoized_exp_two_pi_i_a_div_width(two_pi_i_a, width): #We want to avoid computing this exponential for all m so we precompute this result
    return (two_pi_i_a/width).exp()

def get_minus_one_pow(n):
    """
    Return (-1)**n without actually performing the multiplications.
    """
    if n%2 == 0:
        return 1
    else:
        return -1

def get_cusp_width_from_var_name(var_name):
    """
    Given a string of the form 'q_N' where N is an integer, return N.
    """
    return int(var_name[2:])

def get_m(cusp_normalizer,coset_representative):
    """
    Determine the shift m s.t. \gamma_j = \gamma_c*T**m
    """
    if cusp_normalizer[0] != 0:
        m = (coset_representative[1]-cusp_normalizer[1])/cusp_normalizer[0]
    else:
        m = (coset_representative[3]-cusp_normalizer[3])/cusp_normalizer[2]
    if isinstance(m,int) or m.is_integer==True:
        m = int(m)
    else:
        raise ArithmeticError("Cannot represent m as integer!")
    return m

def get_coset_expansions(F):
    """
    Compute expansion of F for all coset representatives.
    """
    coset_expansions = dict()
    G = F.G
    cosets = G.coset_reps()
    weight = F.weight
    for ci in range(G.ncusps()):
        c = G.cusps()[ci]
        cusp_expansion = F.get_cusp_expansion(c)
        CF = cusp_expansion[0].parent()
        width = G._vertex_data[ci]['width']
        R = cusp_expansion.parent()
        q = R.gen()
        if ci == 0:
            cusp_expansion *= width**(-weight/2)
        roots_of_unity = [(2*CF(0,pi)*i/width).exp() for i in range(width)]
        cusp_normalizer = G.cusp_normalizer(c)
        for coset_i in G._vertex_data[ci]['coset']:
            m = get_m(cusp_normalizer,cosets[coset_i])
            coset_expansion = 0
            for n in range(cusp_expansion.degree()+1):
                coeff = cusp_expansion[n]*roots_of_unity[(n*m)%width]
                coset_expansion += coeff*q**n
            coset_expansions[coset_i] = coset_expansion.O(cusp_expansion.degree()+1)
    
    return coset_expansions

def compute_petersson_product_haberland(F,G):
    """
    Compute Petersson product by using the Haberland-type formula as described in https://arxiv.org/abs/1809.10908
    It is important that F is a cuspform (G can be either cuspidal or non-cuspidal).
    """
    weight = F.weight
    index = F.G.index()
    F_coset_exp, G_coset_exp = get_coset_expansions(F), get_coset_expansions(G)
    CC = F_coset_exp[0][0].parent()
    two_pi_i = CC(0,2*pi)
    rho = (two_pi_i/3).exp()
    scale_fact = 1/(index*(CC(0,2))**(weight-1))
    res = 0
    for j in range(index):
        f_j, g_j = F_coset_exp[j], G_coset_exp[j]
        width = get_cusp_width_from_var_name(f_j.variable())
        for n in range(weight-2+1):
            # term = (-1)**n*binomial(weight-2,n)*I(weight-2-n,rho+1,"ioo",f_j,width,two_pi_i)*conjugate(I(n,CC(0,1),CC(1,1),g_j,width,two_pi_i))
            term = get_minus_one_pow(n)*memoized_binomial(weight-2,n)*I(weight-2-n,rho+1,"ioo",f_j,width,two_pi_i)*(I(n,CC(0,1),CC(1,1),g_j,width,two_pi_i).conjugate())
            res += term
    res *= scale_fact
    return res

def I(n,a,b,coset_expansion,width,two_pi_i):
    """
    Compute I as defined in https://arxiv.org/abs/1809.10908
    We assume that 'a' is a point in the upper half plane and b is either i*infinity or also in H.
    """
    if b == "ioo": #This corresponds to infinity
        if coset_expansion[0] != 0:
            raise ArithmeticError("We have only implemented the case where F is a cuspform!")
        return compute_period_integral(n,a,coset_expansion,width,two_pi_i)
    else:
        m_zero_term = -compute_m_zero_period_integral_summand(n,a,b,coset_expansion) #WHY DO WE NEED THE MINUS?!
        m_larger_zero_sum = compute_period_integral(n,a,coset_expansion,width,two_pi_i)-compute_period_integral(n,b,coset_expansion,width,two_pi_i)
        return m_zero_term+m_larger_zero_sum

def compute_period_integral(n,a,coset_expansion,width,two_pi_i):
    """
    Evaluate int_a^(i*infinity) tau^n F(tau) dtau. We only conside the terms m>0 and need to treat the m=0 case separately.
    """
    res = 0
    for i in range(1,coset_expansion.degree()+1):
        m = QQ(i)/width
        res += coset_expansion[i]*compute_period_integral_term(n,a,m,two_pi_i)
    return res

@cached_function #This function can be used for all cosets in a cusp so it makes sense to re-use it
def compute_period_integral_term(n,a,m,two_pi_i):
    """
    Evaluate int_a^(i*infinity) tau^n exp(2*pi*I*m*tau) dtau.
    We use the formula given in https://wstein.org/books/modform/modform/periods.html#approximating-period-integrals
    """
    res = 0
    a_pow_n = memoized_pow(a,n) #a**(n-s) below now comes almost for free
    two_pi_i_m = two_pi_i*m
    for s in range(n+1):
        #res += ((-1)**s*a**(n-s)/(two_pi_i*m)**(s+1))*prod(range(n+1-s,n+1))
        res += (get_minus_one_pow(s)*memoized_pow(a,n-s)/memoized_pow(two_pi_i_m,s+1))*memoized_prod(n+1-s,n+1)
    res *= memoized_pow(memoized_exp_two_pi_i_a_div_width(two_pi_i*a,m.denominator()),m.numerator())
    return res

def compute_m_zero_period_integral_summand(n,a,b,coset_expansion):
    """
    Compute a_0 * int_a^b t^n dt where a,b are in H.
    """
    return coset_expansion[0]*(b**(n+1)-a**(n+1))/(n+1)
