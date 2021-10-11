"""
Routines for loading noncong_reps.sobj
"""

from psage.modform.arithgroup.mysubgroup import MySubgroup

def get_list_of_all_subgroups(max_index, only_genus_zero=False):
    if max_index > 17:
        raise ArithmeticError("The database of subgroups only goes up to and including order 17!")
    data = load("noncong_reps.sobj")
    res = []
    sorted_signatures = sorted(data.keys())
    for signature in sorted_signatures:
        if signature[0] > max_index:
            break
        if only_genus_zero == True and signature[1] != 0:
            continue
        signature_reps = data[signature]
        for rep in signature_reps:
            (permS,permR,_) = rep
            G = MySubgroup(o2=permS,o3=permR)
            res.append(G)
    return res
