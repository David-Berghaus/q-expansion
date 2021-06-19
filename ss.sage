#
# Transformation of cusp forms: Berghaus -> Selander, Stroembergson
#
# cusp enumeration
#
# SS  cusp     Berghaus
#
# 1 ->   0  -> 0
# 2 ->  -1  -> 2
# 3 ->   7  -> 3
# 4 ->  oo  -> 1

f_0 = load('F_0')
f_1 = load('F_1')

ss = dict()

ss[Cusp(0)] = 1
ss[Cusp(-1)] = 2
ss[Cusp(7)] = 3
ss[Cusp(1, 0)] = 4

# coefficients of the a_1 for each cusp

v_0 = vector([f[1] for f in f_0.values()])
v_1 = vector([f[1] for f in f_1.values()])

# F_0 = a*v_0 + b*v_1
#
# a*v_0[0] + b*v_1[0] = 0
# a*v_0[2] + b*v_1[2] = 1

det = v_0[0]*v_1[2] - v_0[2]*v_1[0]

a = -v_1[0] / det
b =  v_0[0] / det

F_0 = dict()

for k in f_0.keys():
    F_0[k] = a*f_0[k] + b*f_1[k]

# check that the a_1 coefficients of F_0 have the correct value

print("Coefficient a_1 of F_0 at cusps:")

for k in F_0.keys():
    print(ss[k], CC(F_0[k][1]))
#
# F_1
#
# c*v0[0] + d*v_1[0] = 1
# c*v0[1] + d*v_1[1] = 0

det = v_0[0]*v_1[1] - v_0[1]*v_1[0]

c =  v_1[1] / det
d = -v_0[1] / det

F_1 = dict()

for k in f_0.keys():
    F_1[k] = c*f_0[k] + d*f_1[k]

print("Coefficient a_1 of F_1 at cusps:")

for k in F_0.keys():
    print(ss[k], CC(F_1[k][1]))
