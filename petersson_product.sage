def gp_delta():
    #Gamma0(1) weight 12 cusp form
    gp('D=mfDelta(); mf=mfinit(D); DS=mfsymbol(mf,D)')
    print(gp('mfpetersson(DS)'))

    # #Gamma0(4) weight 8 cusp form
    # gp('mf=mfinit([2,8],1); L=mfbasis(mf); S=mfsymbol(mf,L[1])') #Flag '1' indicates that we are considering the complete cuspform space
    # print(gp('mfcoefs(S,8)'))
    # print(gp('mfpetersson(S)'))

    #Gamma0(3) weight 6 cusp form
    gp('mf=mfinit([3,6],0); L=mfbasis(mf); S=mfsymbol(mf,L[1])')
    print(gp('mfcoefs(S,8)'))
    print(gp('mfpetersson(S)'))

def petersson_product_horocycle(F, G, Y):
    """
    We use the approach suggested in https://mathoverflow.net/a/115134
    """
    bit_prec = F['bit_prec']
    weight = F['weight']
    RF = RealField(bit_prec)
    res = 0
    for c in F['cusp_expansions'].keys(): #Iterate over cusp classes
        cusp_res = 0
        f, g = F['cusp_expansions'][c], G['cusp_expansions'][c]
        cusp_width = get_cusp_width_from_var_name(f.variable())
        for n in range(1,f.prec()): #We always assume that at least one of f, g is a cuspform (i.e. has c_0 = 0)
            cusp_res += f[n] * conjugate(g[n]) * exp(RF(-4*pi*n*Y/cusp_width))
        cusp_res *= cusp_width * Y**weight
        res += cusp_res
    return res

def petersson_product_nelson_collins(F, G):
    """
    We use the notation of https://arxiv.org/pdf/1809.10908.pdf section 4.2
    """
    bit_prec = F['bit_prec']
    weight = F['weight']
    RF = RealField(bit_prec)
    res = 0
    Q = 200 #The order up to which we expand Phi
    PI = RF(pi)
    for c in F['cusp_expansions'].keys():
        cusp_res = 0
        f, g = F['cusp_expansions'][c], G['cusp_expansions'][c]
        cusp_width = get_cusp_width_from_var_name(f.variable())
        for n in range(1,f.prec()): #We always assume that at least one of f, g is a cuspform (i.e. has c_0 = 0)
            n_w = RF(n)/cusp_width
            print(n_w)
            cusp_res += f[n] * conjugate(g[n]) * W(4*PI*sqrt(n_w),weight,Q,RF) / (n_w**(weight-1))
        cusp_res *= cusp_width
        res += cusp_res
    res *= 4*(8*PI)**(-(weight-1))
    return res

def W(x, weight, Q, RF): #Note that this is an inefficient way of computing W
    res = 0
    for m in range(1,Q+1):
        y = m*x
        res += y**(weight-1)*RF(y*bessel_K(weight-2,y)-bessel_K(weight-1,y))
    return res

def get_cusp_width_from_var_name(var_name):
    """
    Given a string of the form 'q_N' where N is an integer, return N.
    """
    return int(var_name[2:])
