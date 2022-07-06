from belyi.number_fields import is_effectively_zero, lindep, to_K

def example():
    from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
    from psage.modform.arithgroup.mysubgroup import MySubgroup
    from classes.fourier_expansion import get_cuspform_q_expansion_approx
    #G, Kv = MySubgroup(Gamma0(11)), NumberField(x,'v')
    #G = MySubgroup(o2='(2 3)(1 4)(5 7)(6 10)(8 13)(9 16)(11 19)(12 22)(14 25)(15 28)(17 31)(18 20)(21 26)(23 29)(24 32)(27 33)(30 34)(35 36)',o3='(2 4 3)(1 10 7)(5 16 13)(6 22 19)(8 28 25)(9 20 31)(11 26 18)(12 32 29)(14 33 21)(15 34 23)(17 27 24)(30 35 36)')
    #G, Kv = MySubgroup(o2='(1)(2 5)(3 7)(4 8)(6 9)',o3='(1 2 6)(3 8 5)(4 9 7)'), NumberField(x,'v')
    G, Kv = MySubgroup(o2='(1)(2)(3 12)(4 7)(5 9)(6 10)(8 11)',o3='(1 2 3)(4 8 12)(5 10 7)(6 11 9)'), NumberField(x**2-3,'v',embedding=-1.732)
    #G, Kv = MySubgroup(o2='(10)(2)(3 12)(4 7)(5 9)(6 1)(8 11)',o3='(10 2 3)(4 8 12)(5 1 7)(6 11 9)'), NumberField(x**2-3,'v',embedding=-1.732)
    digit_prec = 100
    f = get_cuspform_q_expansion_approx(AutomorphicFormSpace(G,2),digit_prec) #Weight 2 cuspform to 50 digits precision
    period_lattice_els = compute_period_lattice_els(G,f, digit_prec)
    tmp = get_w_1_w_2(period_lattice_els, digit_prec)
    if tmp == None:
        return None
    w_1, w_2 = tmp
    tau = w_1/w_2
    g_2, g_3 = elliptic_invariants = get_elliptic_invariants(tau)
    print("elliptic_invariants: ", elliptic_invariants)
    j_fl = 1728*g_2**3/(g_2**3-27*g_3**2)
    print("j_fl: ", j_fl)
    #j = QQ(pari(j_fl.real()).bestappr(int(digit_prec/2)))
    j = recognize_expr(j_fl,Kv)
    if j == None:
        return None
    print("j: ", j)
    print("j_fl-j: ", abs(j_fl-j_fl.parent(j)))
    E = EllipticCurve_from_j(j).global_minimal_model()
    print("E: ", E)
    return E

# def example2():
#     from psage.modform.arithgroup.mysubgroup import MySubgroup
#     G, Kv = MySubgroup(o2='(1)(2 5)(3 8)(4 9)(6 12)(7 10)(11)',o3='(1 2 6)(3 9 5)(4 10 8)(7 11 12)'), NumberField(x**2-3,'v',embedding=-1.732)
#     v = Kv.gen()
#     j = 14621235235888115443/16708992677662604*v - 20728089694692551503/16708992677662604
#     E = EllipticCurve_from_j(j).global_minimal_model()
#     return E

def example3():
    from psage.modform.arithgroup.mysubgroup import MySubgroup
    G = MySubgroup(o2='(1)(2 5)(3 8)(4 9)(6 12)(7 10)(11)',o3='(1 2 6)(3 9 5)(4 10 8)(7 11 12)')
    C = CyclotomicField(12)
    z_12 = C.gen()
    v = -z_12**3+2*z_12
    j = 14621235235888115443/16708992677662604*v - 20728089694692551503/16708992677662604
    E = EllipticCurve_from_j(j).global_minimal_model()
    return E

def recognize_expr(expr_fl, K):
    expr = to_K(expr_fl,K)
    if expr == None:
        print("Not enough precision to recognize expression ", expr_fl)
        return None
    if K.degree() == 1:
        expr = QQ(expr) #Change j to QQ because sage is able to reduce the curve then
    return expr

def I(f, tau_0, cusp):
    """
    Evaluate int_{tau_0}^{i*infty} 2*pi*i f(tau) dtau.
    """
    CC = f.get_cusp_expansion(cusp).base_ring()
    RR = RealField(CC.precision())
    cusp_expansion = f.get_cusp_expansion(cusp)
    q = exp(CC(0,2*pi)*tau_0/f.G.cusp_width(cusp))
    res = 0
    for n in range(1,cusp_expansion.prec()):
        res += -cusp_expansion[n]/n * q**n
    return res

def get_tau_0(c, d, CC):
    """
    Get a good choice of the base point to evaluate the period integral.
    """
    return CC(-d,1)/c

def get_best_cusp_and_gamma(gamma, G):
    """
    Given gamma, select the best choice of c to compute P(f,gamma) together with the corresponding gamma-map.
    """
    heights = {}
    gamma_primes = {}
    for cusp in G.cusps():
        gamma_prime = G.cusp_normalizer(cusp).inverse().matrix()*gamma*G.cusp_normalizer(cusp).matrix()
        c = gamma_prime[1][0]
        if c == 0:
            continue
        if c < 0:
            gamma_prime = -gamma_prime
            c = gamma_prime[1][0]
        heights[cusp] = 1/(c*G.cusp_width(cusp))
        gamma_primes[cusp] = gamma_prime
    max_vals = list(heights.values())
    keys = list(heights.keys())
    cusp = keys[max_vals.index(max(max_vals))]
    gamma_prime = gamma_primes[cusp]
    return cusp, gamma_prime

def P(f, gamma):
    cusp, gamma = get_best_cusp_and_gamma(gamma,f.G)
    CC = f.get_cusp_expansion(cusp).base_ring()
    a, b, c, d = gamma[0][0], gamma[0][1], gamma[1][0], gamma[1][1]
    tau_0 = get_tau_0(c,d,CC)
    gamma_of_tau_0 = (a*tau_0+b)/(c*tau_0+d)
    return I(f,tau_0,cusp)-I(f,gamma_of_tau_0,cusp)

def compute_period_lattice_els(G, f, digit_prec):
    """
    Given a group G and the weight two cuspform f, compute the numerical approximations of the period lattice.
    """
    gamma_list = get_gamma_list(G)
    period_lattice_els = [P(f,gamma) for gamma in gamma_list]
    non_zero_unique_period_lattice_els = []
    for lattice_element in period_lattice_els:
        if is_effectively_zero(lattice_element,0.9*digit_prec) == True: #Filter out zero elements
            continue
        is_unique = True
        for non_zero_el in non_zero_unique_period_lattice_els:
            if is_effectively_zero(abs(lattice_element-non_zero_el),0.9*digit_prec) == True: #We don't want duplicates
                is_unique = False
                break
        if is_unique == True:
            non_zero_unique_period_lattice_els.append(lattice_element)
    if len(non_zero_unique_period_lattice_els) == 0:
        raise ArithmeticError("Period lattice elements are either empty or zero!")
    return non_zero_unique_period_lattice_els

def get_gamma_list(G):
    gamma_list = []
    for gen in G.generators():
        c = gen[1][0]
        if c == 0:
            continue #P(f) is zero for this so we can ignore it
        elif c < 0:
            minus_id = SL2Z([-1, 0, 0, -1])
            gamma_list.append(minus_id*gen) #Note that we work in PSL(2,ZZ)
        else:
            gamma_list.append(gen)
    return gamma_list

def get_w_1_w_2(non_zero_unique_period_lattice_els, digit_prec):
    if len(non_zero_unique_period_lattice_els) < 2:
        print("Not enough elements to form lattice basis.")
        return None
    elif len(non_zero_unique_period_lattice_els) == 2: #We don't need to use LLL to reduce the basis
        w_1, w_2 = non_zero_unique_period_lattice_els[0], non_zero_unique_period_lattice_els[1]
    else:
        w_1 = non_zero_unique_period_lattice_els[0]
        for i in range(1,len(non_zero_unique_period_lattice_els)):
            if lindep([w_1,non_zero_unique_period_lattice_els[i]],check_if_result_is_invalid=True) == None:
                w_2 = non_zero_unique_period_lattice_els[i] #w_1 and w_2 are linearly independent
                break
        lattice_factors = get_lattice_factors(w_1, w_2, non_zero_unique_period_lattice_els, digit_prec)
        if lattice_factors == None:
            return None
        w_1, w_2 = improve_w_1_w_2(w_1, w_2, lattice_factors)
    if imag(w_1/w_2) > 0:
        return w_1, w_2
    else:
        return w_2, w_1

def get_lattice_factors(w_1, w_2, non_zero_unique_period_lattice_els, digit_prec):
    """
    Compute the lattice factors n,m in QQ such that lattice_el = n*w_1+m*w_2. 
    """
    K = NumberField(x,'v')
    digit_prec_half = int(digit_prec/2)
    M = MatrixSpace(w_1.real().parent(),2,2)
    A = M([[w_1.real(),w_2.real()],[w_1.imag(),w_2.imag()]])
    n_m_list = [(QQ(1),QQ(0)),(QQ(0),QQ(1))]
    for lattice_el in non_zero_unique_period_lattice_els:
        if lattice_el != w_1 and lattice_el != w_2:
            n_fl, m_fl = A.solve_right(vector([lattice_el.real(),lattice_el.imag()]))
            #n, m = QQ(pari(n_fl).bestappr(digit_prec_half)), QQ(pari(m_fl).bestappr(digit_prec_half)) #Use heuristic upper bound in denominator to avoid overfitting due to incorrect digits
            n, m = recognize_expr(n_fl,K), recognize_expr(m_fl,K)
            if n == None or m == None:
                return None
            n_m_list.append((n,m))
    return n_m_list

def improve_w_1_w_2(w_1, w_2, lattice_factors):
    """
    Try to improve lattice so that all elements can be written as ZZ*w_1+ZZ*w_2. 
    """
    if len(lattice_factors) == 0: #Nothing to improve here
        return w_1, w_2
    largest_common_denominator = 1
    for lattice_factor in lattice_factors:
        tmp_n = lattice_factor[0]*largest_common_denominator
        if tmp_n.denominator() != 1:
            largest_common_denominator *= tmp_n.denominator()
        tmp_m = lattice_factor[1]*largest_common_denominator
        if tmp_m.denominator() != 1:
            largest_common_denominator *= tmp_m.denominator()
    integer_lattice_factors = [[largest_common_denominator*lattice_factor[0],largest_common_denominator*lattice_factor[1]] for lattice_factor in lattice_factors]
    V = ZZ^2
    smallest_lattice_factor = V.submodule(integer_lattice_factors).basis()
    w_1, w_2 = (smallest_lattice_factor[0][0]*w_1+smallest_lattice_factor[0][1]*w_2)/largest_common_denominator, (smallest_lattice_factor[1][0]*w_1+smallest_lattice_factor[1][1]*w_2)/largest_common_denominator
    return w_1, w_2

def reduce_tau_to_psl2z(w_1, w_2): #Actually we don't need this function because arb does it internally
    tau = w_1/w_2
    if real(tau) < -1/2:
        w_1, w_2 = w_1+w_2, w_2
        return reduce_tau_to_psl2z(w_1,w_2)
    if real(tau) > 1/2:
        w_1, w_2 = w_1-w_2, w_2
        return reduce_tau_to_psl2z(w_1,w_2)
    if abs(tau) > 1:
        w_1, w_2 = -w_2, w_1
        return reduce_tau_to_psl2z(w_1,w_2)
    return tau

def get_elliptic_invariants(tau):
    """
    Use Arb to compute g_2 and g_3.
    """
    CC = tau.parent()
    CBF = ComplexBallField(CC.precision())
    res = CBF(tau).elliptic_invariants()
    return (CC(res[0]),CC(res[1]))