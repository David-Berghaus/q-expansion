from sage.rings.qqbar import QQbar
from sage.modular.cusps import Cusp

from belyi.newton_genus_zero import run_newton

class BelyiMap():
    """
    Class for working with Belyi maps that are expressed in terms of Factored_Polynomials with algebraic coefficients.
    """
    def __init__(self, S, starting_digit_prec=42, target_digit_prec=50000, max_extension_field_degree=None):
        """
        Compute Belyi map from scratch. This is done by first using Hejhal's method to get approximate values at the elliptic points and cusps.
        Afterwards, Newton iterations are used to refine the solution. In the last step, the result is verified.
        """
        G = S.group()
        if G.genus() != 0:
            raise ArithmeticError("This function only works for genus zero subgroups!")
        if S.weight() != 0:
            raise ArithmeticError("This function only works for weight zero!")
        if max_extension_field_degree == None:
            max_extension_field_degree = G.index() #For the groups that we are considering the conjugacy class size is <= the index
        
        (p3, p2, pc) = run_newton(S, starting_digit_prec, target_digit_prec, stop_when_coeffs_are_recognized=True)
        self.G = G
        self.p3, self.p2, self.pc = p3, p2, pc
        self.p3_constructed, self.p2_constructed, self.pc_constructed = p3.construct(), p2.construct(), pc.construct()

        #If these 4 steps succeed, then the result is verified
        self.ell_2_point_evaluations = self.get_elliptic_two_point_evaluations()
        self.ell_3_point_evaluations = self.get_elliptic_three_point_evaluations()
        self.cusp_evaluations = self.get_cusp_evaluations()
        self.verify_polynomial_equation()

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return self.p3.__str__() + " / " + self.pc.__str__()
    
    def get_elliptic_two_point_evaluations(self):
        """
        Returns algebraic evaluations of hauptmodul at elliptic points of order two.
        """
        ell_2_point_evaluations = []
        for (p,multiplicity) in self.p2.factors:
            roots = p.roots(ring=QQbar)
            for (root, order) in roots:
                if order != 1:
                    raise ArithmeticError("Something went wrong, we need distinct roots here!")
                else:
                    ell_2_point_evaluations.append([root,multiplicity])
        return ell_2_point_evaluations

    def get_elliptic_three_point_evaluations(self):
        """
        Returns algebraic evaluations of hauptmodul at elliptic points of order three.
        """
        ell_3_point_evaluations = []
        for (p,multiplicity) in self.p3.factors:
            roots = p.roots(ring=QQbar)
            for (root, order) in roots:
                if order != 1:
                    raise ArithmeticError("Something went wrong, we need distinct roots here!")
                else:
                    ell_3_point_evaluations.append([root,multiplicity])
        return ell_3_point_evaluations
    
    def get_cusp_evaluations(self):
        """
        Returns algebraic evaluations of hauptmodul at cusps.
        The cusp at infinity gets treated separately.
        """
        cusp_evaluations = []
        for (p,multiplicity) in self.pc.factors:
            roots = p.roots(ring=QQbar)
            for (root, order) in roots:
                if order != 1:
                    raise ArithmeticError("Something went wrong, we need distinct roots here!")
                else:
                    cusp_evaluations.append([root,multiplicity])
        return cusp_evaluations

    def verify_polynomial_equation(self):
        p3_constructed, p2_constructed, pc_constructed = self.p3_constructed, self.p2_constructed, self.pc_constructed
        if p3_constructed-p2_constructed-1728*pc_constructed != 0: #Verify result
            raise ArithmeticError("Verification of polynomial_equation failed!")
        if pc_constructed.degree()+self.G.cusp_width(Cusp(1,0)) != self.G.index():
            raise ArithmeticError("Wrong behavior at infinity!")

    # def get_x_expansion(self, trunc_order):
    #     """
    #     Get expansion of rational function in terms of x (the hauptmodul).
    #     """
    #     R = self.p3_constructed/self.pc_constructed

    def get_q_expansion(self, trunc_order):
        raise NotImplementedError("This functionality has not been implemented yet!")
    
    def get_approx_q_expansion(self, trunc_order, digit_prec):
        raise NotImplementedError("This functionality has not been implemented yet!")
