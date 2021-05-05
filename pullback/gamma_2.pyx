cdef extern from "complex.h":
    cdef double cimag(double complex)
    cdef double creal(double complex)
    cdef double cabs(double complex)

from sage.modular.arithgroup.congroup_sl2z import SL2Z
from sage.modular.arithgroup.arithgroup_perm import sl2z_word_problem

from psage.modform.arithgroup.mysubgroups_alg import SL2Z_elt

cpdef gamma_2_pb_dp(double x, double y):
    """
    Maps point into fundamental domain of Gamma(2)
    """
    T = SL2Z_elt(1,2,0,1)
    T_inv = SL2Z_elt(1,-2,0,1)
    T_0 = SL2Z_elt(1,0,-2,1)
    T_0_inv = SL2Z_elt(1,0,2,1)
    pb_map = SL2Z_elt(1,0,0,1)
    cdef double complex z_pb = x + y*1j

    while True:
        if abs(creal(z_pb)) > 1:
            if creal(z_pb) < 0: #We need to move point to the right
                z_pb += 2
                pb_map = T*pb_map
            else: #We need to move point to the left
                z_pb -= 2
                pb_map = T_inv*pb_map
        elif cabs(z_pb+0.5) < 0.5:
            z_pb = z_pb/(2*z_pb+1)
            pb_map = T_0_inv*pb_map
        elif cabs(z_pb-0.5) < 0.5:
            z_pb = z_pb/(-2*z_pb+1)
            pb_map = T_0*pb_map
        else:
            break

    return creal(z_pb), cimag(z_pb), pb_map[0], pb_map[1], pb_map[2], pb_map[3]
