from sage.rings.power_series_poly import PowerSeries_poly

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
        width_int = G._vertex_data[ci]['width']
        width = CF(width_int,0)
        R = cusp_expansion.parent()
        q = R.gen()
        cusp_expansion *= width**(-weight/2) #We assume that the cusp-expansions of F have the width absorbed so we need to rescale
        roots_of_unity = [exp(2*CF(0,pi)*i/width) for i in range(width_int)]
        cusp_normalizer = G.cusp_normalizer(c)
        for coset_i in G._vertex_data[ci]['coset']:
            m = get_m(cusp_normalizer,cosets[coset_i])
            coset_expansion = PowerSeries_poly(R,prec=cusp_expansion.degree()+1) #This corresponds to O(q**M_0) in sage syntax
            for n in range(cusp_expansion.degree()+1):
                coeff = cusp_expansion[n]*roots_of_unity[(n*m)%width_int]
                coset_expansion += coeff*q**n
            coset_expansions[coset_i] = coset_expansion
    
    return coset_expansions

def compute_petersson_product_haberland(F,G):
    """
    Compute Petersson product by using the Haberland-type formula as described in https://arxiv.org/abs/1809.10908
    It is important that F is a cuspform (G can be either cuspidal or non-cuspidal).
    """
    weight = F.weight
    index = F.G.index()
    F_coset_exp, G_coset_exp = get_coset_expansions(F), get_coset_expansions(G)
    CF = F_coset_exp[0][0].parent()
    two_pi_i = CF(0,2*pi)
    rho = exp(two_pi_i/3)
    scale_fact = 1/(index*(CF(0,2))**(weight-1))
    res = 0
    for j in range(index):
        f_j, g_j = F_coset_exp[j], G_coset_exp[j]
        width = CF(get_cusp_width_from_var_name(f_j.variable()),0)
        for n in range(weight-2+1):
            term = (-1)**n*binomial(weight-2,n)*I(weight-2-n,rho+1,"ioo",f_j,width,two_pi_i)*conjugate(I(n,CF(0,1),CF(1,1),g_j,width,two_pi_i))
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
        m = i/width
        res += coset_expansion[i]*compute_period_integral_term(n,a,m,two_pi_i)
    return res

def compute_period_integral_term(n,a,m,two_pi_i):
    """
    Evaluate int_a^(i*infinity) tau^n exp(2*pi*I*m*tau) dtau.
    We use the formula given in https://wstein.org/books/modform/modform/periods.html#approximating-period-integrals
    """
    res = 0
    for s in range(n+1):
        res += ((-1)**s*a**(n-s)/(two_pi_i*m)**(s+1))*prod([j for j in range(n+1-s,n+1)])
    res *= exp(two_pi_i*m*a)
    return res

#This formula diverges!
# def compute_period_integral_term(n,a,m,two_pi_i):
#     """
#     Evaluate int_a^(i*infinity) tau^n exp(2*pi*I*m*tau) dtau.
#     """
#     return (-1)**n*exp(-two_pi_i*m*a)*(factorial(n)/(two_pi_i*m)**(n+1))*R(n,-two_pi_i*m*a)

# def R(n,x):
#     """
#     Compute R as defined in https://arxiv.org/abs/1809.10908
#     """
#     res = 0
#     for i in range(n+1):
#         res += x**i/factorial(i) #This is a very inefficient way of computing these sums...
#     print(x, res)
#     return res

def compute_m_zero_period_integral_summand(n,a,b,coset_expansion):
    """
    Compute a_0 * int_a^b t^n dt where a,b are in H.
    """
    return coset_expansion[0]*(b**(n+1)-a**(n+1))/(n+1)