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
    return prod(range(a,b)) #We could probably also try to optimize this by computing it recursively but have not done this yet

@cached_function
def memoized_pow(a, n):
    """
    Compute a**n recursively to cache intermediate results.
    Warning: Note that the default maximum recursion depth in python is 1000 (and by increasing it one risks stackoverflows).
    This function should therefore only be used for relatively small powers.
    """
    if n == 0:
        return 1
    if n == 1:
        return a
    return a*memoized_pow(a,n-1)

@cached_function
def memoized_exp_two_pi_i_a_div_width(two_pi_i_a, width): #We want to avoid computing this exponential for all m so we precompute this result
    return (two_pi_i_a/width).exp()

def clear_memoized_caches(): #Clear caches to decrease memory usage and avoid interference with additional computations
    memoized_binomial.clear_cache()
    memoized_prod.clear_cache()
    memoized_pow.clear_cache()
    memoized_exp_two_pi_i_a_div_width.clear_cache()
    compute_period_integral_term.clear_cache()

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

def get_m(cusp_normalizer, coset_representative):
    """
    Determine the shift m s.t. \gamma_j = \gamma_c*T**m
    """
    if cusp_normalizer[0] != 0:
        m = (coset_representative[1]-cusp_normalizer[1])/cusp_normalizer[0]
    else:
        m = (coset_representative[3]-cusp_normalizer[3])/cusp_normalizer[2]
    if isinstance(m,int) or m.is_integer() == True:
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
        cusp_expansion /= width**(weight//2) #We follow the convention of Cohen's paper who uses different cusp-normalizers
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

class Hashabledict(dict): #See: https://stackoverflow.com/a/16162138
    """
    Note that after using this class the dictionaries should not be mutated!
    """
    def __hash__(self):
        return hash(frozenset(self))

def get_exp_two_pi_i_a_m_dict(a_values, f, g, CC):
    """
    Precompute exp(2*pi*I*a*m) where m = i/cusp_width for i in range(trunc_order) efficiently through recursive multiplications.
    We do not use "memoized_pow" for this task because it might exceed the max_recursion_depth.
    Additionally, storing the values like this should require less memory usage.
    """
    G = f.G #The subgroup
    exp_two_pi_i_a_m_dict = dict()
    for a in a_values:
        exp_two_pi_i_a_m_dict[a] = dict()
        exp_two_pi_i_a_m_dict[a][QQ(0)] = 1 #The zeroth power is trivial
        for ci in range(G.ncusps()):
            cusp_width = G._vertex_data[ci]['width']
            c = G.cusps()[ci]
            f_deg, g_deg = f.get_cusp_expansion(c).degree(), g.get_cusp_expansion(c).degree()
            exp_two_pi_i_a_div_cusp_width = (CC(0,2*pi)*(a/cusp_width)).exp()
            tmp = 1
            for i in range(1,max(f_deg,g_deg)+1):
                m = QQ(i)/cusp_width
                tmp *= exp_two_pi_i_a_div_cusp_width
                exp_two_pi_i_a_m_dict[a][m] = tmp #Some of these will already be defined but it does not seem worth it to check whether they exist...
    return Hashabledict(exp_two_pi_i_a_m_dict) #We need our dict to be hashable in order to use it as input for memoized functions

def get_exp_two_pi_i_a_m(a, m, exp_two_pi_i_a_m_dict):
    try:
        res = exp_two_pi_i_a_m_dict[a][m]
    except KeyError: #It sometimes happens that the modular forms have slightly different amounts of coefficients
        CC = a.parent()
        res = (CC(0,2*pi)*a*m).exp()
        exp_two_pi_i_a_m_dict[a][m] = res
    return res

def compute_petersson_product_haberland(f, g, clear_memoized_caches_bool=True, exp_two_pi_i_a_m_dict=None):
    """
    Compute Petersson product by using the Haberland-type formula as described in https://arxiv.org/abs/1809.10908
    It is important that f is a cuspform (g can be either cuspidal or non-cuspidal).
    """
    weight = f.weight
    index = f.G.index()
    f, g = f._convert_to_CC(), g._convert_to_CC() #Working with CBFs seems to cause problems with the hashing
    f_coset_exp, g_coset_exp = get_coset_expansions(f), get_coset_expansions(g)
    CC = f_coset_exp[0][0].parent()
    two_pi_i = CC(0,2*pi)
    rho = (two_pi_i/3).exp()
    a_values = [rho+1,CC(0,1),CC(1,1)]
    if exp_two_pi_i_a_m_dict == None:
        exp_two_pi_i_a_m_dict = get_exp_two_pi_i_a_m_dict(a_values,f,g,CC)
    scale_fact = 1/(index*(CC(0,2))**(weight-1))
    res = 0
    for j in range(index):
        f_j, g_j = f_coset_exp[j], g_coset_exp[j]
        width = get_cusp_width_from_var_name(f_j.variable())
        for n in range(weight-2+1):
            # term = (-1)**n*binomial(weight-2,n)*I(weight-2-n,rho+1,"ioo",f_j,width,two_pi_i)*conjugate(I(n,CC(0,1),CC(1,1),g_j,width,two_pi_i))
            term = get_minus_one_pow(n)*memoized_binomial(weight-2,n)*I(weight-2-n,rho+1,"ioo",f_j,width,two_pi_i,exp_two_pi_i_a_m_dict)*(I(n,CC(0,1),CC(1,1),g_j,width,two_pi_i,exp_two_pi_i_a_m_dict).conjugate())
            res += term
    res *= scale_fact
    if clear_memoized_caches_bool == True:
        clear_memoized_caches()
    return res

def I(n, a, b, coset_expansion, width, two_pi_i, exp_two_pi_i_a_m_dict):
    """
    Compute I as defined in https://arxiv.org/abs/1809.10908
    We assume that 'a' is a point in the upper half plane and b is either i*infinity or also in H.
    """
    if b == "ioo": #This corresponds to infinity
        if coset_expansion[0] != 0:
            raise ArithmeticError("We have only implemented the case where F is a cuspform!")
        return compute_period_integral(n,a,coset_expansion,width,two_pi_i,exp_two_pi_i_a_m_dict)
    else:
        m_zero_term = -compute_m_zero_period_integral_summand(n,a,b,coset_expansion) #WHY DO WE NEED THE MINUS?!
        m_larger_zero_sum = compute_period_integral(n,a,coset_expansion,width,two_pi_i,exp_two_pi_i_a_m_dict)-compute_period_integral(n,b,coset_expansion,width,two_pi_i,exp_two_pi_i_a_m_dict)
        return m_zero_term+m_larger_zero_sum

def compute_period_integral(n, a, coset_expansion, width, two_pi_i, exp_two_pi_i_a_m_dict):
    """
    Evaluate int_a^(i*infinity) tau^n F(tau) dtau. We only conside the terms m>0 and need to treat the m=0 case separately.
    """
    res = 0
    for i in range(1,coset_expansion.degree()+1):
        m = QQ(i)/width
        res += coset_expansion[i]*compute_period_integral_term(n,a,m,two_pi_i,exp_two_pi_i_a_m_dict)
    return res

@cached_function
def compute_period_integral_term(n, a, m, two_pi_i, exp_two_pi_i_a_m_dict):
    """
    Evaluate int_a^(i*infinity) tau^n exp(2*pi*I*m*tau) dtau.
    We use the formula given in https://wstein.org/books/modform/modform/periods.html#approximating-period-integrals
    """
    res = 0
    two_pi_i_m = two_pi_i*m
    for s in range(n+1):
        #res += ((-1)**s*a**(n-s)/(two_pi_i*m)**(s+1))*prod(range(n+1-s,n+1))
        res += (get_minus_one_pow(s)*memoized_pow(a,n-s)/memoized_pow(two_pi_i_m,s+1))*memoized_prod(n+1-s,n+1)
    #res *= exp(two_pi_i*m*a)
    res *= get_exp_two_pi_i_a_m(a,m,exp_two_pi_i_a_m_dict)
    return res

def compute_m_zero_period_integral_summand(n, a, b, coset_expansion):
    """
    Compute a_0 * int_a^b t^n dt where a,b are in H.
    """
    return coset_expansion[0]*(b**(n+1)-a**(n+1))/(n+1)
