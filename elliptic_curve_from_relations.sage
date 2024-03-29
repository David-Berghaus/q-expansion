# ****************************************************************************
#       Copyright (C) 2022 David Berghaus <berghaus@th.physik.uni-bonn.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

load("run.sage")
load("passports.sage")

def example_relations():
    #Kv = NumberField(x,'v')
    G = MySubgroup(o2='(1)(2 5)(3 7)(4 8)(6 9)',o3='(1 2 6)(3 8 5)(4 9 7)')
    pp = compute_passport_data_higher_genera([G,None],100,100,6,compute_embeddings=False,return_floating_expansions=False)
    [X, Y, Z] = pp["q_expansions"][6]["cuspforms_raw"] #X,Y,Z now correspond to weight 6 cuspforms
    Kw = X.base_ring()
    w = Kw.gen()
    monomial_list = get_monomial_list(X,Y,Z)
    A = get_monomial_matrix(monomial_list)
    ker = get_kernel(A)
    verify_equation(ker,X,Y,Z)
    R = PolynomialRing(Kw,"x, y, z")
    x, y, z = R.gens()
    monomial_equation = get_monomial_equation(R,ker)
    monomial_equation = monomial_equation.subs(x=w*x,z=z/w)
    monomial_equation /= w
    monomial_equation = monomial_equation.change_ring(QQ)
    E = EllipticCurve_from_cubic(monomial_equation,[1,0,0],morphism=False) 
    E = E.global_minimal_model()
    return E

def get_monomial_list(X, Y, Z):
    """
    Form all degree 3 combinations of XYZ and return terms as a list.
    """
    return [X**3, X**2*Y, X**2*Z, X*Y**2, X*Y*Z, X*Z**2, Y**3, Y**2*Z, Y*Z**2, Z**3]

def get_monomial_matrix(monomial_list):
    """
    Return monomial matrix in a way that its kernel corresponds to the monomial coefficients.
    """
    min_val = min([monomial.valuation() for monomial in monomial_list])
    M = MatrixSpace(monomial_list[0].base_ring(),len(monomial_list),len(monomial_list))
    A = list(M.zero())
    for monomial_index in range(len(monomial_list)):
        for coeff_index in range(min_val,len(monomial_list)+min_val):
            A[monomial_index][coeff_index-min_val] = monomial_list[monomial_index][coeff_index]
    return M(A)

def get_kernel(A):
    return A.kernel() #Is there a better way?

def get_monomial_equation(R, kernel):
    """
    Write the monomial equation in terms of the generators of A3.
    """
    ker_basis = kernel.basis()[0]
    x, y, z = R.gens()
    monomial_equation = ker_basis[0]*x**3 + ker_basis[1]*x**2*y + ker_basis[2]*x**2*z + ker_basis[3]*x*y**2 + ker_basis[4]*x*y*z + ker_basis[5]*x*z**2 + ker_basis[6]*y**3 + ker_basis[7]*y**2*z + ker_basis[8]*y*z**2 + ker_basis[9]*z**3
    return monomial_equation

def verify_equation(kernel, x, y, z):
    """
    Confirm the monomial equation by plugging in the corresponding forms.
    """
    ker_basis = kernel.basis()[0]
    equation = ker_basis[0]*x**3 + ker_basis[1]*x**2*y + ker_basis[2]*x**2*z + ker_basis[3]*x*y**2 + ker_basis[4]*x*y*z + ker_basis[5]*x*z**2 + ker_basis[6]*y**3 + ker_basis[7]*y**2*z + ker_basis[8]*y*z**2 + ker_basis[9]*z**3
    if equation != 0:
        raise ArithmeticError("Monomial equation does not hold!")