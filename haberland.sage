from sage.rings.power_series_poly import PowerSeries_poly

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
        width = CF(G._vertex_data[ci]['width'],0)
        R = cusp_expansion.parent()
        q = R.gen()
        if ci != 0: #We assume that the cusp-expansions of F have the width absorbed so we need to rescale
            fact = width**(-weight/2)
            cusp_expansion *= fact
        zeta_w = exp(2*CF(0,pi)/width)
        cusp_normalizer = G.cusp_normalizer(c)
        for coset_i in G._vertex_data[ci]['coset']:
            m = get_m(cusp_normalizer,cosets[coset_i])
            coset_expansion = PowerSeries_poly(R,prec=cusp_expansion.degree()+1) #This corresponds to O(q**M_0) in sage syntax
            for n in range(cusp_expansion.degree()+1):
                coeff = cusp_expansion[n]*zeta_w**(n*m) #It would be faster to work with mod here instead of generic powers
                coset_expansion += coeff*q**n
            coset_expansions[coset_i] = coset_expansion
    
    return coset_expansions