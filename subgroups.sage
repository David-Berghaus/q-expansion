"""
Routines for loading data/noncong_reps.sobj
"""

from psage.modform.arithgroup.mysubgroup import MySubgroup
from psage.groups.permutation_alg import MyPermutation

def get_conjugator_permutation(o_start, o_end):
    """
    Return permutation p, s.t. p**(-1)*o_start*p == o_end.
    Note that p is not unique!
    """
    o_s = o_start.cycles_ordered_as_list()
    o_e = o_end.cycles_ordered_as_list()
    p_as_list = [o_e[o_s.index(i)] for i in range(1,len(o_s)+1)]
    p = MyPermutation(p_as_list)
    assert p.inverse()*o_start*p == o_end
    return p

def get_permT_conjugator(permT, reverse_permT):
    """
    Given permT, return conjugator p, such that p**(-1)*permT*p has cycle lens ordered by size and labels sorted by size
    unless there are several cusps with the same smallest width in which case we try to set a unique cusp at infinity.
    If reverse_permT is True, we place the largest cusp at infinity (unless it has higher multiplicity).
    This is useful for higher genera cases because it typically leads to smaller coefficients.
    For genus zero cases it is not very useful because it leads to slower rigorous arithmetic.
    """
    if reverse_permT == False:
        ct = get_desired_permT_cycle_type(permT)
    else:
        ct = get_desired_permT_cycle_type_reversed(permT)
    i = 1
    normalized_perm_T_str = ''
    for cycle_len in ct:
        normalized_perm_T_str += '('
        for j in range(i,i+cycle_len):
            normalized_perm_T_str += str(j) + ","
        normalized_perm_T_str += ')' #We have an empty komma at the end which MyPermutation however seems to ignore
        i += cycle_len
    normalized_perm_T = MyPermutation(normalized_perm_T_str)
    p = get_conjugator_permutation(permT,normalized_perm_T)
    return p

def get_desired_permT_cycle_type(permT):
    """
    Returns the desired cycle type given permT.
    If the cusp-width of the smallest cusp-width is unique, this is given by the traditional cycle type
    (i.e., the cycle lens sorted by their size).
    Otherwise we try to place a cusp at infinity that has a unique cusp-width.
    """
    ct = permT.cycle_type()
    for i in range(len(ct)):
        if has_equal_list_entry(ct,i) == False: #Found a unique cusp-width
            ct[0], ct[i] = ct[i], ct[0] #Swap cusps to put the unique one at infinity
            break
    #In case no unique cusp has been found, we could also check to place the one with the lowest multiplicity at infinity.
    #This case however never seems to happen for the subgroups that we consider so we dont implement it...
    return ct

def get_desired_permT_cycle_type_reversed(permT):
    """
    Returns the desired cycle type given permT.
    For this function we try to place the largest cusp at infinity (if it is unique).
    """
    ct = permT.cycle_type()
    ct.reverse() #Note that this is done inplace!
    for i in range(len(ct)):
        if has_equal_list_entry(ct,i) == False: #Found a unique cusp-width
            ct[0], ct[i] = ct[i], ct[0] #Swap cusps to put the unique one at infinity
            break
    return ct

def has_equal_list_entry(list, index):
    """
    If there exists an element in list (outside index) that is equal to list[index], return True, otherwise return False.
    """
    for i in range(len(list)):
        if i != index and list[i] == list[index]:
            return True
    return False

def get_list_of_all_passports(max_index, genus=None, reverse_permT=True):
    """
    Return a list of all passports up to specified index and with specified genus.
    """
    if max_index > 17:
        raise ArithmeticError("The database of subgroups only goes up to and including order 17!")
    data = load("data/noncong_reps.sobj")
    res = []
    sorted_signatures = sorted(data.keys())
    for signature in sorted_signatures:
        if signature[0] > max_index:
            break
        if genus != None and signature[1] != genus: #Only search for groups with specified genus
            continue
        signature_reps = data[signature]
        sorted_by_cycle_types = dict()
        for rep in signature_reps:
            (permS_str,permR_str,_) = rep
            permS, permR = MyPermutation(permS_str), MyPermutation(permR_str)
            permT = permS*permR
            p = get_permT_conjugator(permT,reverse_permT)
            permS_new, permR_new, permT_new = p.inverse()*permS*p, p.inverse()*permR*p, p.inverse()*permT*p
            cycle_types = (tuple(permS_new.cycle_type()),tuple(permR_new.cycle_type()),tuple(permT_new.cycle_type())) #We need tuples for dict keys
            G = MySubgroup(o2=permS_new,o3=permR_new)
            monodromy_group_order = G.perm_group().order() #It is probably better to recognize groups by order instead of some string...
            if cycle_types in sorted_by_cycle_types: #Dict entry already exists
                sorted_by_monodromy = sorted_by_cycle_types[cycle_types]
                if monodromy_group_order in sorted_by_monodromy: #Dict entry already exists
                    sorted_by_monodromy[monodromy_group_order].append(G)
                else:
                    sorted_by_monodromy[monodromy_group_order] = [G]
            else:
                sorted_by_cycle_types[cycle_types] = dict()
                sorted_by_monodromy = sorted_by_cycle_types[cycle_types]
                sorted_by_monodromy[monodromy_group_order] = [G]
        #We have now sorted all representatives into passports
        for sorted_by_monodromy in sorted_by_cycle_types.values():
            for passport in sorted_by_monodromy.values():
                res.append(passport)
    return res
