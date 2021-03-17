import math
import numpy as np

from sage.modular.arithgroup.congroup_sl2z import SL2Z_class
from sage.modular.arithgroup.congroup_gamma import Gamma_class
from sage.rings.complex_arb cimport *
from sage.rings.real_arb cimport *
from sage.libs.arb.arb cimport *
from sage.rings.complex_arb import ComplexBallField
from sage.rings.real_arb import RealBallField

from psage.modform.arithgroup.mysubgroup import MySubgroup
from psage.modform.arithgroup.mysubgroups_alg import apply_sl2z_map_dp, pullback_general_group_dp, normalize_point_to_cusp_dp, SL2Z_elt
from psage.modform.maass.automorphic_forms import AutomorphicFormSpace

cpdef my_pullback_pts_dp(S,int Qs,int Qf,double Y): #Slow version
    G = S.group()
    if G.minimal_height() <= Y:
        raise ArithmeticError, "Need smaller value of Y!"
    Xm = dict()
    twoQl = 2*(Qf-Qs+1)
    Pb = dict()
    for ci in G._cusps:
        cii = G._cusps.index(ci)
        Pb[cii] = dict()
        for cj in G._cusps:
            cjj = G._cusps.index(cj)
            Pb[cii][cjj] = []
    for m in range(Qs,Qf+1):
        Xm[m] = (2.0*m-1)/twoQl

    for ci in G._cusps:
        cii = G._cusps.index(ci)
        for j in range(Qs,Qf+1):
            z_horo = Xm[j]+Y*1j
            x,y = normalize_point_to_cusp_dp(G,ci,Xm[j],Y)
            x1,y1,T_a,T_b,T_c,T_d = pullback_general_group_dp(G,x,y,ret_mat=1)
            vjj = G.closest_vertex(x1,y1,as_integers=1)
            cjj = G._vertex_data[vjj]['cusp'] #closest cusp
            U_w = G._vertex_data[vjj]['cusp_map']
            x2, y2 = apply_sl2z_map_dp(x1,y1,U_w[0],U_w[1],U_w[2],U_w[3]) #apply U_w
            cj = G._cusps[cjj]
            x3, y3 = normalize_point_to_cusp_dp(G,cj,x2,y2,inv=1)
            z3 = x3+y3*1j
            tmp = G.cusp_normalizer(cj).inverse()*U_w*SL2Z_elt(T_a,T_b,T_c,T_d)*G.cusp_normalizer(ci)
            swi = math.sqrt(G.cusp_width(ci))
            swj = math.sqrt(G.cusp_width(cj))
            tmp = simple_two_by_two_matmul(tmp,[swi, 0, 0, 1/swi])
            tmp = simple_two_by_two_matmul([1/swj, 0, 0, swj],tmp)
            a,b,c,d = tmp[0],tmp[1],tmp[2],tmp[3]
            Pb[cii][cjj].append((z_horo,a,b,c,d,z3))
    return Pb

cdef apply_moebius_transformation_dp(z,a,b,c,d): #z is a double complex and a,b,c,d are doubles
    return (a*z+b)/(c*z+d)

cpdef my_pullback_pts_arb_wrap(S,int Qs,int Qf,Y,prec): #Slow version
    RR = RealBallField(prec)
    CC = ComplexBallField(prec)
    G = S.group()
    if G.minimal_height() <= float(Y):
        raise ArithmeticError, "Need smaller value of Y!"
    Xm = dict()
    twoQl = 2*(Qf-Qs+1)
    Pb = dict()
    for ci in G._cusps:
        cii = G._cusps.index(ci)
        Pb[cii] = dict()
        for cj in G._cusps:
            cjj = G._cusps.index(cj)
            Pb[cii][cjj] = []
    for m in range(Qs,Qf+1):
        Xm[m] = RR(2.0*m-1)/twoQl
    for ci in G._cusps:
        cii = G._cusps.index(ci)
        if is_int(G.cusp_width(ci)) == False:
            raise NameError('Cannot convert cusp width to Arb without sacrifizing precision!')
        for j in range(Qs,Qf+1):
            z_horo = CC(Xm[j],Y)
            x,y = normalize_point_to_cusp_dp(G,ci,float(Xm[j]),float(Y)) #We don't need to perform this computation with arbs...
            x1,y1,T_a,T_b,T_c,T_d = pullback_general_group_dp(G,x,y,ret_mat=1)
            vjj = G.closest_vertex(x1,y1,as_integers=1)
            cjj = G._vertex_data[vjj]['cusp'] #closest cusp
            U_w = G._vertex_data[vjj]['cusp_map']
            x2, y2 = apply_sl2z_map_dp(x1,y1,U_w[0],U_w[1],U_w[2],U_w[3]) #apply U_w
            cj = G._cusps[cjj]
            x3, y3 = normalize_point_to_cusp_dp(G,cj,x2,y2,inv=1)
            z3 = x3+y3*1j
            tmp = G.cusp_normalizer(cj).inverse()*U_w*SL2Z_elt(T_a,T_b,T_c,T_d)*G.cusp_normalizer(ci) #This computation involves only integers
            swi = RR(G.cusp_width(ci)).sqrt()
            swj = RR(G.cusp_width(cj)).sqrt()
            tmp = simple_two_by_two_matmul(tmp,[swi, 0, 0, 1/swi])
            tmp = simple_two_by_two_matmul([1/swj, 0, 0, swj],tmp)
            a,b,c,d = tmp[0],tmp[1],tmp[2],tmp[3]
            Pb[cii][cjj].append((z_horo,a,b,c,d,z3)) #These are all arb (acb) types despite z3
    return Pb

cdef apply_moebius_transformation_arb_wrap(z,a,b,c,d): #z is a ComplexBall and a,b,c,d are RealBalls
    return (a*z+b)/(c*z+d)

cpdef simple_two_by_two_matmul(A1,A2): #Currently slow because of arbitrary types and python-arrays
    r"""
        Performs a two-by-two matrix multiplication of two two-dimensional arrays
        Works with Arrays of type A = [a,b,c,d] that represent two-by-two matrices [[a,b],[c,d]]
    """
    return [A1[0]*A2[0]+A1[1]*A2[2], A1[0]*A2[1]+A1[1]*A2[3], A1[2]*A2[0]+A1[3]*A2[2], A1[2]*A2[1]+A1[3]*A2[3]]

cdef is_int(val):
    if type(val) == int:
        return True
    else:
        if val.is_integer():
            return True
        else:
            return False

def compare_pb_to_psage(S,Qs,Qf,Y):
    my_pb = my_pullback_pts_dp(S,Qs,Qf,Y)
    from psage.modform.maass.pullback_algorithms import pullback_pts_dp
    tmp = pullback_pts_dp(S,Qs,Qf,Y)
    G = S.group()
    eps = 1e-14

    psage_pb = dict()
    for ci in G._cusps:
        cii = G._cusps.index(ci)
        psage_pb[cii] = dict()
        for cj in G._cusps:
            cjj = G._cusps.index(cj)
            psage_pb[cii][cjj] = []
            ypb = tmp["ypb"][cii][cjj]
            for i in range(len(ypb)):
                if ypb[i] != 0:
                    psage_pb[cii][cjj].append(ypb[i])
    
    nc = G.ncusps()
    for cii in range(nc):
        for cjj in range(nc):
            if len(my_pb[cii][cjj]) != len(psage_pb[cii][cjj]):
                print("Error: ", len(my_pb[cii][cjj]), ", ", len(psage_pb[cii][cjj]))
            for i in range(len(my_pb[cii][cjj])):
                if abs(my_pb[cii][cjj][i][5].imag-psage_pb[cii][cjj][i]) > eps:
                    print(my_pb[cii][cjj][i][5].imag, ", ", psage_pb[cii][cjj][i])

    return psage_pb



#print(my_pullback_pts_dp(AutomorphicFormSpace(Gamma_class(2)),-3,3,0.1))