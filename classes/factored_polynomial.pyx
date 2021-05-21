from copy import copy

cpdef construct_poly_from_root_tuple(x, root_tuple): #We might want to use native arb here later
    p = 1
    (roots, order) = root_tuple
    for root in roots:
        p *= (x-root)
    return [p,order]

cpdef construct_poly_from_sorted_factors(factors):
    reduced_factors = [] #By reduced factors we mean that all factors of the same order get multiplied
    reduced_factors.append(factors[0])
    for i in range(1,len(factors)):
        (poly_fact, order) = factors[i]
        if reduced_factors[-1][1] == order:
            reduced_factors[-1][0] *= poly_fact
        else:
            reduced_factors.append(factors[i])
    res = 1
    for (poly_fact, order) in reduced_factors:
        res *= poly_fact**order
    return res

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
    
    def derivative(self, poly_index, coeff_index):
        """
        Returns a polynomial that corresponds to the derivative with respect to the 'coeff_index'th coefficient of polynomial at 'poly_index'.
        """
        outer_derivative = copy(self.factors)
        outer_derivative[poly_index][1] -= 1
        outer_derivative.sort(key=lambda tup: tup[1])
        outer_derivative = construct_poly_from_sorted_factors(outer_derivative)
        inner_derivative = (self.factors[poly_index][1])*(self.polygen)**coeff_index
        derivative = inner_derivative*outer_derivative #this multiplication by a monomial can certainly be optimized


