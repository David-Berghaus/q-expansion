load("run.sage")

#These notations are part of the Fiori, Franc paper. Note that we do not get the same permT...
G1 = MySubgroup(o2='(1 2)(3 4)(5 6)(7)',o3='(1)(2 3 5)(4 6 7)')
u_G1 = (-7)**(1/4)/7**2
H1 = MySubgroup(o2='(1 2)(3 4)(5 6)(7)',o3='(1)(2 3 5)(7 6 4)')
u_H1 = (-7**3)**(1/5)/7**2
z3 = exp(2*pi*I/3)
u_U1 = ((1763*z3 + 1255)*2**2*3/7**7)**(1/6)
U1 = MySubgroup(o2='(1 2)(3 4)(5)(6 7)',o3='(1)(2 3 5)(4 6 7)')

def get_coeffs():
    digit_prec=400
    bit_prec = digits_to_bits(digit_prec-10)
    tmp = get_coefficients_eisenstein_ir_arb_wrap(AutomorphicFormSpace(U1,4),digit_prec)
    c = tmp._get_mcbd(bit_prec)
    a = []
    CF = ComplexField(bit_prec)
    for i in range(15):
        a.append(N(CF(c[i][0])/u_U1**(i+1),digits=digit_prec-10))
    for i in range(len(a)):
        print(a[i].algebraic_dependency(10))
        print('')
    return a