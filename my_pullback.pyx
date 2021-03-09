import math
import numpy as np

from sage.modular.arithgroup.congroup_sl2z import SL2Z_class
from sage.modular.arithgroup.congroup_gamma import Gamma_class

from psage.modform.arithgroup.mysubgroup import MySubgroup
from psage.modform.arithgroup.mysubgroups_alg import apply_sl2z_map_dp, pullback_general_group_dp, normalize_point_to_cusp_dp, SL2Z_elt
from psage.modform.maass.automorphic_forms import AutomorphicFormSpace

SL2Z = SL2Z_class()

cpdef my_pullback_pts_dp(S,Qs,Qf,Y): #Slow version
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
            #print(ci, ", ", j)
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
            # print(x3, ", ", y3)
            tmp = G.cusp_normalizer(cj).inverse()*U_w*SL2Z_elt(T_a,T_b,T_c,T_d)*G.cusp_normalizer(ci)
            tmp = np.matrix([[tmp[0],tmp[1]],[tmp[2],tmp[3]]])
            swi = math.sqrt(G.cusp_width(ci))
            swj = math.sqrt(G.cusp_width(cj))
            tmp = np.matmul(tmp,np.matrix([[swi, 0], [0, 1/swi]]))
            tmp = np.matmul(np.matrix([[1/swj, 0], [0, swj]]),tmp)
            a,b,c,d = tmp[0,0],tmp[0,1],tmp[1,0],tmp[1,1]
            Pb[cii][cjj].append((z_horo,a,b,c,d,z3))
            # zm = Xm[j]+Y*1j
            # print((tmp[0,0]*zm+tmp[0,1])/(tmp[1,0]*zm+tmp[1,1]))
            # print("")
    return Pb

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