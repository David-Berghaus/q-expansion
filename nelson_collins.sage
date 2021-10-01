def gp_petersson():
    gp('default(realprecision, 120)')
    #Gamma0(1) weight 12 cusp form
    gp('D=mfDelta(); mf=mfinit(D); DS=mfsymbol(mf,D)')
    print(gp('mfpetersson(DS)'))
    # # gp('mf=mfinit([1,12],3); L=mfbasis(mf); E=L[1]; ES=mfsymbol(mf,E)') #Flag '3' indicates that we are considering the Eisenstein space
    # gp('mf=mfinit([1,12],4); L=mfbasis(mf); E=L[1]; ES=mfsymbol(mf,E)')
    # print(gp('mfcoefs(E,8)'))
    # print(gp('mfpetersson(DS,ES)'))

    # #Gamma0(2) weight 8 cusp form
    # gp('mf=mfinit([2,8],1); L=mfbasis(mf); S=mfsymbol(mf,L[1])') #Flag '1' indicates that we are considering the complete cuspform space
    # print(gp('mfcoefs(S,8)'))
    # print(gp('mfpetersson(S)'))

    # #Gamma0(3) weight 6 cusp form
    # gp('mf=mfinit([3,6],0); L=mfbasis(mf); f=L[1]; fS=mfsymbol(mf,f)') #Flag '0' indicates that we are considering the new cuspform space
    # print(gp('mfcoefs(f,8)'))
    # print(gp('mfpetersson(fS)'))
    # print(gp('mfslashexpansion(mf,f,[0,-1;1,0],8,0)'))

def get_cusp_width_from_var_name(var_name):
    """
    Given a string of the form 'q_N' where N is an integer, return N.
    """
    return int(var_name[2:])

def compute_petersson_product_nelson_collins(F, G):
    """
    We use the notation of https://arxiv.org/pdf/1802.09740.pdf theorem 4.2
    But fixed the arguments of the Bessel functions...
    """
    bit_prec = F['bit_prec']
    epsilon = N(2**(-bit_prec),10)
    weight = F['weight']
    RF = RealField(bit_prec)
    res = 0
    PI = RF(pi)
    for c in F['cusp_expansions'].keys():
        cusp_res = 0
        f, g = F['cusp_expansions'][c], G['cusp_expansions'][c]
        cusp_width = get_cusp_width_from_var_name(f.variable())
        for n in range(1,min(f.prec(),g.prec())): #We always assume that at least one of f, g is a cuspform (i.e. has c_0 = 0)
            x = 4*PI*sqrt(RF(n)/cusp_width)
            cusp_res += f[n] * conjugate(g[n]) * W(weight, x, epsilon) / n**(weight-1)
        res += cusp_res #note that the width is already inside the normalizer
    return 4*(8*PI)**(-(weight-1))*res/F['index']

def W(weight, x, epsilon): #There are faster ways to compute this function, see https://arxiv.org/pdf/1809.10908.pdf section 4.3
    res = 0
    summand = 10.0 #Dummy initialization
    m = 1
    #The convergence criterion is non-rigorous but since (f[n]*conjugate(g[n])/n**(weight-1)) grows like n,
    #it should work decently in practice.
    while abs(summand) > epsilon:
        y = m*x
        #We perform the besselk computations in arb because sage's implementation is incredibly slow
        summand = pow(y,weight-1)*(y*arb_besselk(weight-2,y,epsilon)-arb_besselk(weight-1,y,epsilon))
        res += summand
        m += 1
    return res

def arb_besselk(nu, z, epsilon):
    """
    For a ComplexFloat z, evaluate bessel_k(nu,z) using Arb.
    Since the bessel functions in arb are sometimes unstable, we increase the precision until sufficient precision is reached.
    """
    init_prec = z.prec()
    arb_prec = int(1.5*init_prec)
    CF = ComplexField(init_prec)
    CBF = ComplexBallField(arb_prec)
    res = CBF(z).bessel_K(nu)
    while res.rad() > epsilon:
        arb_prec = int(1.5*arb_prec)
        CBF = ComplexBallField(arb_prec)
        res = CBF(z).bessel_K(nu)
    return CF(res)