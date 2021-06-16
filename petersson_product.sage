def gp_petersson():
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

def petersson_product_nelson_collins(F, G):
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
        for n in range(1,f.prec()): #We always assume that at least one of f, g is a cuspform (i.e. has c_0 = 0)
            n_w = RF(n)/cusp_width
            cusp_res += f[n] * conjugate(g[n]) * Phi(weight, epsilon, PI, n_w) / n**(weight-1)
        res += cusp_res #note that the width is already inside the normalizer
    return 4*res/F['index']

def Phi(weight, epsilon, PI, n_w): #There are faster ways to compute this function, see https://arxiv.org/pdf/1809.10908.pdf section 4.3
    res = 0
    summand = 10.0 #Dummy initialization
    m = 1
    #The convergence criterion is non-rigorous but since (f[n]*conjugate(g[n])/n**(weight-1)) grows like n,
    #it should work decently in practice.
    while abs(n_w*summand) > epsilon:
        x = 4*PI*m*sqrt(n_w)
        summand = (x/(8*PI))**(weight-1)*(x*bessel_K(weight-2,x)-bessel_K(weight-1,x))
        res += summand
        m += 1
    return res


#Now to the bruteforce approach

two_pi = RR(2*pi)
two_pi_i = two_pi*I

@CachedFunction
def e_x(x, w, m):
    return exp(two_pi_i*x*m/w)

@CachedFunction
def e_y(y, w, m, weight):
    return y**(weight-2)*exp(-two_pi*y*m/w)

def get_psl2z_grid(Q, dr, max_r):
    """
    Return a grid of points in a fundamental domain of PSL(2,Z). The points are sampled along arcs parallel to the x^2+y^2 = 1.
    We choose Q values of -0.5<x<0.5 and afterwards sample values of y until sqrt(x^2+y^2) < max_r.
    """
    r = 1
    dx = RR(1/(Q+1))
    x_vals = [i*dx-0.5 for i in range(1,Q+1)]
    y_rows = []
    while r < max_r:
        y_row = [sqrt(r**2-x**2) for x in x_vals]
        y_rows.append(y_row)
        r += dr
    return x_vals, y_rows

def bruteforce_integration(f, g, weight):
    """
    Compute Petersson product by numerically solving the double integral.
    This function is so far only implemented for SL(2,ZZ).
    """
    w = 1
    res = 0
    Q = 50
    dr = 0.01
    max_r = 1.4
    x_vals, y_rows = get_psl2z_grid(Q, dr, max_r)
    amount_of_sampling_points = len(x_vals)**2*len(y_rows)
    for n in range(1,f.prec()): #We always assume that at least one of f, g is a cuspform (i.e. has c_0 = 0)
        for m in range(1,g.prec()): #We always assume that at least one of f, g is a cuspform (i.e. has c_0 = 0)
            #Compute the double integrals
            summand = 0
            for x in x_vals:
                for y_row in y_rows:
                    for y in y_row:
                        summand += e_y(y, w, n+m, weight)*e_x(x, w, n-m)
            summand /= amount_of_sampling_points

            #Multiply cuspform coefficients
            summand *= f[n] * conjugate(g[m])

            res += summand
    return res/RR(pi/3) #scale by volume
            