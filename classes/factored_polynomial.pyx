import copy

cpdef construct_poly_from_root_tuple(x, root_tuple): #We might want to use native arb here later
    p = 1
    (roots, order) = root_tuple
    for root in roots:
        p *= (x-root)
    return [p,order]

class Factored_Polynomial():
    """
    Class for working with polynomials that are factored like:
    p = p_1^(n_1)*...*p_N^(n_N)
    """
    def __init__(self, polygen, root_tuples):
        """
        root_tuples contains a list of tuples. 
        Each tuple is given by a list of ComplexBalls reflecting the roots as well as the order of each of the roots.
        """
        self.polygen = polygen
        factors = []
        for root_tuple in root_tuples:
            factors.append(construct_poly_from_root_tuple(polygen,root_tuple))
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


