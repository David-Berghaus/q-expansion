"""
Routines for loading data/noncong_reps.sobj
"""
from json.encoder import INFINITY
import time

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

def get_list_of_all_passports(max_index, genus=None, reverse_permT=True, optimize_gen_1_for_period_int=True):
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
            permS_new, permR_new, permT_new = get_conjugate_permutation(permS,p), get_conjugate_permutation(permR,p), get_conjugate_permutation(permT,p)
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
    if optimize_gen_1_for_period_int == True: #Try to find a better choice of G that optimizes the period integral evaluation
        for (i,pp) in enumerate(res):
            if pp[0].genus() == 1:
                pp_conj = optimize_passport_for_period_integral(pp) #Conjugate all subgroups by the best conjugation
                res[i] = pp_conj #Update res
    return res

def get_conjugate_permutation(sigma, conjugator):
    """
    Conjugate sigma by conjugator, i.e., return conjugator^(-1)*sigma*conjugator.
    """
    return conjugator.inverse()*sigma*conjugator

def get_fixed_one_permutation_iterators(d):
    """
    Return an iterator of all permutations in S_d that keep 1 fixed.
    The output is an iterator over lists.
    """
    return Permutations(list(range(2,d+1)))

def optimize_passport_for_period_integral(passport):
    """
    Conjugate passport so that its first element is optimized for period integral evaluations.
    """
    pp = passport #For the sake of simpler notation
    conj_and_heights = [get_optimal_conjugator_and_max_height_for_period_integral(G) for G in pp]
    best_index = 0 #Index of the G in pp with optimal properties for period integral evaluation
    best_conj, best_height = conj_and_heights[best_index]
    for j in range(1,len(conj_and_heights)):
        new_conj, new_height = conj_and_heights[j]
        if new_height > best_height:
            best_index = j
            best_conj, best_height = new_conj, new_height
    pp.insert(0,pp.pop(best_index)) #Push best choice of G to front
    #Conjugate all subgroups by the best conjugation
    pp_conj = [MySubgroup(o2=get_conjugate_permutation(G.permS,best_conj),o3=get_conjugate_permutation(G.permR,best_conj)) for G in pp]
    return pp_conj

def get_optimal_conjugator_and_max_height_for_period_integral(G, preserve_cusp_width_at_inf=True):
    """
    Given a subgroup G, return permutation sigma, such that sigma^(-1)*G*sigma is optimal as well as the
    height under this choice.
    Optimal in this setting means that the period integrals can be as high inside H as possible
    If preserve_cusp_width_at_inf = True, we make sure that the cusp width at infinity remains the same.
    Otherwise we look for a potentially even better configuration.
    """
    max_height = get_min_height_for_period_integral(G)
    max_conj = MyPermutation(list(range(1,G.index()+1))) #The one element in S_d
    o2, o3 = G.permS, G.permR
    if preserve_cusp_width_at_inf == True:
        pot_swaps = G.permT.cycles()[0][1:] #All labels that we could potentially swap with 1
    else:
        pot_swaps = list(range(2,G.index()+1))
    for i in pot_swaps:
        l = list(range(1,G.index()+1))
        l[1-1], l[i-1] = l[i-1], l[1-1] #Swap 1 and i
        o = MyPermutation(l)
        G_p = MySubgroup(get_conjugate_permutation(o2,o),get_conjugate_permutation(o3,o))
        min_height_p = get_min_height_for_period_integral(G_p)
        if min_height_p > max_height:
            max_height, max_conj = min_height_p, o
    return max_conj, max_height

def get_min_height_for_period_integral(G):
    min_height = Infinity
    for cusp in G.cusps():
        N, N_inv = G.cusp_normalizer(cusp).matrix(), G.cusp_normalizer(cusp).inverse().matrix()
        cusp_width = G.cusp_width(cusp)
        for gen in G.generators():
            gamma = N_inv*gen*N
            c = gamma[1][0]
            if c == 0:
                continue
            height = 1/(abs(c)*cusp_width)
            if height < min_height:
                min_height = height
    if min_height == Infinity:
        raise ArithmeticError("Min height is infinity which should not happen.")
    return min_height