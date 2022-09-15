from sage.all_cmdline import *   # import sage library

import os
import sys
from pathlib import Path

load("run.sage")

storage_path = "/cephfs/user/s6ddberg/H_5"

def get_M_2_H_5_single_label(digit_prec, maxiter, label=0):
    from classes.fourier_expansion import FourierExpansion, c_vec_to_cusp_expansions
    from point_matching.point_matching_arb_wrap import _get_normalization_modforms
    
    S = AutomorphicFormSpace(H_5,2)

    if digit_prec < 60:
        raise ArithmeticError("Digit precision must be at least 60.")
    if label > 9:
        raise ArithmeticError("Label too high.")
    
    starting_order = 0
    bit_prec = digits_to_bits(digit_prec)
    basis = []
    if label < 9:
        c_vecs, M_0, labels = get_modform_basis_gmres_arb_wrap(S,digit_prec,multiplicity=9,labels=[label],normalization_zeros=[1],return_M_and_labels=True,maxiter=maxiter)
        normalization = _get_normalization_modforms(S,9,label=labels[0])
    else:
        c_vecs, M_0, _ = get_modform_basis_gmres_arb_wrap(S,digit_prec,multiplicity=11,labels=[10],return_M_and_labels=True,maxiter=maxiter)
        normalization = _get_normalization_modforms(S,11,label=10)
    c_vec_mcbd = c_vecs[0]._get_mcbd(bit_prec)
    cusp_expansions = c_vec_to_cusp_expansions(c_vec_mcbd,S,starting_order,normalization,M_0)
    base_ring = cusp_expansions[Cusp(1,0)].base_ring()
    basis.append(FourierExpansion(S.group(),S.weight(),cusp_expansions,"ModForm",base_ring))

    return basis[0]

def get_M_4_H_5_single_label(digit_prec, maxiter, label=0):
    from classes.fourier_expansion import FourierExpansion, c_vec_to_cusp_expansions
    from point_matching.point_matching_arb_wrap import _get_normalization_modforms
    
    S = AutomorphicFormSpace(H_5,4)

    if digit_prec < 150:
        raise ArithmeticError("Digit precision must be at least 150.")
    if label > 19:
        raise ArithmeticError("Label too high.")
    
    starting_order = 0
    bit_prec = digits_to_bits(digit_prec)
    basis = []
    if label < 19:
        c_vecs, M_0, labels = get_modform_basis_gmres_arb_wrap(S,digit_prec,multiplicity=19,labels=[label],normalization_zeros=[2],return_M_and_labels=True,maxiter=maxiter)
        normalization = _get_normalization_modforms(S,19,label=labels[0])
    else:
        c_vecs, M_0, _ = get_modform_basis_gmres_arb_wrap(S,digit_prec,multiplicity=21,labels=[20],return_M_and_labels=True,maxiter=maxiter)
        normalization = _get_normalization_modforms(S,21,label=20)
    c_vec_mcbd = c_vecs[0]._get_mcbd(bit_prec)
    cusp_expansions = c_vec_to_cusp_expansions(c_vec_mcbd,S,starting_order,normalization,M_0)
    base_ring = cusp_expansions[Cusp(1,0)].base_ring()
    basis.append(FourierExpansion(S.group(),S.weight(),cusp_expansions,"ModForm",base_ring))

    return basis[0]

def get_S_4_H_5_single_label(digit_prec, maxiter, label=0):
    from classes.fourier_expansion import FourierExpansion, c_vec_to_cusp_expansions
    from point_matching.point_matching_arb_wrap import _get_normalization_modforms
    
    S = AutomorphicFormSpace(H_5,4)

    if digit_prec < 150:
        raise ArithmeticError("Digit precision must be at least 150.")
    if label > 9:
        raise ArithmeticError("Label too high.")
    
    starting_order = 1
    bit_prec = digits_to_bits(digit_prec)
    basis = []
    if label < 9:
        c_vecs, M_0, labels = get_cuspform_basis_gmres_arb_wrap(S,digit_prec,multiplicity=9,labels=[label],normalization_zeros=[1],return_M_and_labels=True,maxiter=maxiter)
        normalization = _get_normalization_modforms(S,9,label=labels[0])
    else:
        c_vecs, M_0, _ = get_cuspform_basis_gmres_arb_wrap(S,digit_prec,multiplicity=11,labels=[10],return_M_and_labels=True,maxiter=maxiter)
        normalization = _get_normalization_modforms(S,11,label=10)
    c_vec_mcbd = c_vecs[0]._get_mcbd(bit_prec)
    cusp_expansions = c_vec_to_cusp_expansions(c_vec_mcbd,S,starting_order,normalization,M_0)
    base_ring = cusp_expansions[Cusp(1,0)].base_ring()
    basis.append(FourierExpansion(S.group(),S.weight(),cusp_expansions,"CuspForm",base_ring))

    return basis[0]

def compute_M_k_H_5(weight, digit_prec, maxiter, label):
    if weight == 2:
        res = get_M_2_H_5_single_label(digit_prec, maxiter, label)
    elif weight == 4:
        res = get_M_4_H_5_single_label(digit_prec, maxiter, label)
    else:
        raise ArithmeticError("Weight must be 2 or 4.")
    save(res,os.path.join(storage_path,"m_{}_{}.sobj".format(weight,label)))

def compute_S_k_H_5(weight, digit_prec, maxiter, label):
    if weight == 2:
        res = get_cuspform_basis_approx(AutomorphicFormSpace(H_5,2),digit_prec,labels=[label])
    elif weight == 4:
        res = get_S_4_H_5_single_label(digit_prec, maxiter, label)
    else:
        raise ArithmeticError("Weight must be 2 or 4.")
    save(res,os.path.join(storage_path,"s_{}_{}.sobj".format(weight,label)))

label = int(sys.argv[1])
weight = 2
digit_prec = 60
maxiter = 30

compute_M_k_H_5(weight, digit_prec, maxiter, label)
#compute_S_k_H_5(weight, digit_prec, maxiter, label)
