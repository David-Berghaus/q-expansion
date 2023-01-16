def get_derivative(f):
    """
    Compute 1/(2 pi i) d/dq f(q).
    """
    q = f.parent().gen()
    f_expansion_prime = 0
    for n in range(f.valuation(),f.prec()):
        f_expansion_prime += n*f[n]*q**n
    return f_expansion_prime

def get_monomials(X, Y):
    monomials = []
    for i in range(7):
        monomials.append(X**i)
    # for i in range(4):
    #     monomials.append(X**i*Y)
    monomials.append(Y**2)
    return monomials

def get_coefficient_matrix(monomials):
    min_prec = min([m.prec() for m in monomials])
    max_val = max([m.valuation() for m in monomials])
    L = monomials[0].base_ring()
    M = MatrixSpace(L, len(monomials), min_prec-max_val)
    A = []
    for monomial in monomials:
        coeffs = [monomial[n] for n in range(max_val, min_prec)]
        A.append(coeffs)
    A = M(A)
    return A.transpose()

def get_hyperelliptic_curve_genus_two(f_1, f_2):
    """
    Given a basis of S_2 defined as power series over L, identify the hyperelliptic curve (if it exists).
    """
    X = f_1/f_2
    Y = get_derivative(X)/f_2
    monomials = get_monomials(X, Y)
    A = get_coefficient_matrix(monomials)
    x = var("x")
    y = var("y")
    k = A.right_kernel_matrix()
    L = f_1.base_ring()
    P.<x> = L[]
    res = 0
    for i in range(7):
        res += k[0,i]*x**i
    if A[0,7] == 0:
        return None
    res = res/(-A[0,7])
    C = HyperellipticCurve(res)
    return C

