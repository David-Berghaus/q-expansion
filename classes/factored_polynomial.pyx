from sage.rings.polynomial.polynomial_complex_arb import Polynomial_complex_arb

cpdef construct_poly_from_root_tuple(x, root_tuple): #We might want to use native arb here later
    p = 1
    (roots, order) = root_tuple
    for root in roots:
        p *= (x-root)
    return [p,order]

cpdef construct_poly_from_coeff_tuple(x, coeff_tuple):
    (coeffs, order) = coeff_tuple
    p = Polynomial_complex_arb(x.parent(), coeffs)
    return [p,order]

class Factored_Polynomial():
    """
    Class for working with polynomials that are factored like:
    p = p_1^(n_1)*...*p_N^(n_N)
    """
    def __init__(self, polygen, root_tuples=None, coeff_tuples=None):
        """
        root_tuples contains a list of tuples. 
        Each tuple is given by a list of ComplexBalls reflecting the roots as well as the order of each of the roots.
        Alternatively we can construct the polynomials directly from their coefficients.
        """
        self.polygen = polygen
        factors = []
        if (root_tuples == None and coeff_tuples == None) or (root_tuples != None and coeff_tuples != None):
            raise ArithmeticError("Invalid construction. Construct either through root_tuples or coeff_tuples!")
        if root_tuples != None: #Construct polynomial from root tuples. This usually happens during the first iteration where the haupt-values are passed
            for root_tuple in root_tuples:
                factors.append(construct_poly_from_root_tuple(polygen,root_tuple))
        if coeff_tuples != None:
            for coeff_tuple in coeff_tuples:
                factors.append(construct_poly_from_coeff_tuple(polygen,coeff_tuple))
        self.factors = factors

    def construct(self):
        """
        Explicitly constructs self as a polynomial.
        """
        factors = self.factors
        res = 1
        for (factor, order) in factors:
            res *= factor**order
        return res

    def derivative(self, poly_index, coeff_index):
        """
        Returns a polynomial that corresponds to the derivative with respect to the 'coeff_index'th coefficient of polynomial at 'poly_index'.
        """
        factors = self.factors
        outer_derivative = 1
        for i in range(len(factors)):
            p, order = factors[i]
            if i == poly_index:
                order -= 1
            outer_derivative *= p**order
        inner_derivative = (factors[poly_index][1])*(self.polygen)**coeff_index
        derivative = inner_derivative*outer_derivative #this multiplication by a monomial can certainly be optimized
        return derivative


