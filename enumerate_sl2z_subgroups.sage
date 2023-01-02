# ****************************************************************************
#       Copyright (C) 2023 Danylo Radchenko <danradchenko@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

#This code uses the GAP package grape, but only for the function libgap.SmallestImageSet (written by Steve Linton).
#It also uses bliss for graph.canonical_label (though sage's built-in version of canonical_label will likely also work).
#Install those with:
#sage -i gap_packages bliss

import itertools
from collections import Counter
assert libgap.LoadPackage("grape")

def edge_subsets(m, l, w):
    #Brute-force implementation, can be made much better with a recursive algorithm
    outl = []
    for p in subsets(l):
        c = Counter(flatten(p))
        f = True
        for i in range(m):
            if i in c and c[i] > w[i]:
                f = False
                break
        if f:
            outl.append(p)
    return outl

def sl2_multigraphs(k1, k3):
    #k1 - number of fixed points of pi3; k3 - number of orbits of size 3 in pi3
    #we always assume that k3 >= 1
    n = k1 + 3*k3
    if k3>1:
        graphs_list = list(graphs.nauty_geng("{} -c -d1 -D3".format(k3)))
    else:
        graphs_list = [Graph(1)]
    #set partition to use to restrict permutations allowed for canonical labels
    part = [list(range(k3)), list(range(k3,k3+k1))]
    graphs_set = set()
    #set containing fingerprints of graphs:
    #since graphs are mutable by default, we use tuple(sorted((a,b) for a,b,e in G.canonical_label().edges(sort=True)))
    graphs_output = []
    for G in graphs_list:
        G.allow_loops(True)
        G.allow_multiple_edges(True)
        l0 = [i for i in G for j in range(3-G.degree(i))]
        G.add_vertices(part[1])
        for p1 in itertools.combinations(l0, k1):
            #Add k1 vertices of degree 1
            H1 = G.copy()
            for i,j in enumerate(p1):
                H1.add_edge(j, part[1][i])
            w1 = [3-H1.degree(i) for i in range(k3)]
            l1 = [(i,j) for (i,j,e) in G.edges(sort=True) if w1[i] >= 1 and w1[j] >= 1]
            for p2 in edge_subsets(k3, l1, w1):
                #Add double edges in all possible ways
                H2 = H1.copy()
                for i,j in p2:
                    H2.add_edge(i,j)
                lloops = [i for i in range(k3) if H2.degree(i) == 1]
                for p3 in subsets(lloops):
                    #Add loops in all possible ways
                    H3 = H2.copy()
                    for i in p3:
                        H3.add_edge(i,i)
                    #Only collect non-isomorophic graphs
                    H4 = H3.canonical_label(partition=part)
                    fp = tuple(sorted((a,b) for a,b,e in H4.edges(sort=True)))
                    if not fp in graphs_set:
                        graphs_set.add(fp)
                        graphs_output.append(H4)
    return graphs_output

def pi2s_by_graph(G, k3):
    n = G.num_verts()
    el = [(i,j) for i,j,e in G.edges(sort=True)]
    ed = {(i,j):k for k,(i,j) in enumerate(el)}
    ml = [ed[(i,j)] for i,j,e in G.matching()]
    dallow = {i:({0,1,2} if i<k3 else {0}) for i in range(n)}
    labels = [(i,j) for i in dallow.keys() for j in dallow[i]]
    labeld = {(i,j):k for k,(i,j) in enumerate(labels)}
    #Fixing the values on edges from a matching and all loops
    #this gets rid of most of the symmetries
    sol_d0 = {i:(0,0) for i in ml}
    indl = list(range(len(el)))
    for i in ml:
        dallow[el[i][0]].remove(0)
        dallow[el[i][1]].remove(0)
        indl.remove(i)
    for i,j,e in G.loops():
        sol_d0[ed[i,j]] = (1,2)
        dallow[i] -= {1,2}
        indl.remove(ed[i,j])
    el0 = [el[i] for i in indl]
    if indl:
        outl = []
        #print(indl, el0, dallow)
        for d1 in submatching_problem(indl, el0, dallow):
            d1.update(sol_d0)
            pi2t = tuple(sorted(tuple(sorted([labeld[(el[i][0], t[0])]+1, labeld[(el[i][1], t[1])]+1])) for i,t in d1.items()))
            outl.append(pi2t)
        return outl
    else:
        pi2t = tuple(sorted(tuple(sorted([labeld[(el[i][0], t[0])]+1, labeld[(el[i][1], t[1])]+1])) for i,t in sol_d0.items()))
        return [pi2t]

def submatching_problem(indl, el, dallow):
    if len(indl) == 1:
        return [{indl[0] : (a,b)} for a in dallow[el[0][0]] for b in dallow[el[0][1]]]
    else:
        outl = []
        indl1 = indl[1:]
        el1 = el[1:]
        for a in dallow[el[0][0]]:
            for b in dallow[el[0][1]]:
                dallow1 = deepcopy(dallow)
                dallow1[el[0][0]].remove(a)
                dallow1[el[0][1]].remove(b)
                for dtmp in submatching_problem(indl1, el1, dallow1):
                    dtmp[indl[0]] = (a,b)
                    outl.append(dtmp)
        return outl

def generate_index_n_subgroups(n, PGL2=False):
    if n < 7: #n>6 is necessary with the current assumptions
        raise NotImplementedError("Indexes smaller than 7 are not supported yet.")
    outl = []
    outs = {}
    for k3 in range(1,n//3+1):
        k1 = n-3*k3
        pi3t = [(3*i+1,3*i+2,3*i+3) for i in range(k3)] + [(i+1,) for i in range(3*k3, n)]
        (Aut1, g1, d, pi3a, omega) = canonical_form(n, [], pi3t, precompute=True, PGL2=PGL2)
        graph_list = sl2_multigraphs(k1, k3)
        for G in graph_list:
            pi2l = pi2s_by_graph(G, k3)
            outl += [(pi2, pi3) for pi2, pi3 in canonicalize(Aut1, g1, d, pi3a, omega, pi2l)]
    return outl


def canonical_form(n, pi2, pi3, precompute=False, PGL2=False):
    #Canonical form for a pair of permutations (pi2, pi3)
    perm = Permutation(flatten(sorted(pi3, key=lambda t:-len(t)))).inverse()
    pi3a = sorted([tuple(perm(i) for i in t) for t in pi3])
    Sn = SymmetricGroup(n)
    g3 = Sn("".join(str(t) for t in pi3a if len(t)>1))
    if PGL2:
        #The group generated by pi2 and pi3 is conjugate to the one generated by pi2 and pi3^(-1) in PGL_2(Z)
        Aut = Sn.normalizer(Sn.subgroup([g3]))
    else:
        Aut = Sn.centralizer(g3)
    #omega = [(i,) for i in range(1,n+1)] + [(i,j) for i in range(1,n+1) for j in range(i+1,n+1)]
    #don't include fixed points into the definition of pi2
    omega = [(i,j) for i in range(1,n+1) for j in range(i+1,n+1)]
    omega_d = {a:i+1 for i,a in enumerate(omega)}
    gens_list = Aut.gens()
    gens_list1 = [Permutation([omega_d[tuple(sorted(g(i) for i in t))] for t in omega]) for g in gens_list]
    Aut1 = PermutationGroup(gens_list1)
    if precompute:
        return (Aut1, perm, omega_d, pi3a, omega)
    pi2a = [tuple(sorted(perm(i) for i in t)) for t in pi2]
    l = libgap.SmallestImageSet(Aut1, sorted([omega_d[t] for t in pi2a])).sage()
    return (tuple(omega[i-1] for i in l), tuple(pi3a))


def canonicalize(Aut1, perm, omega_d, pi3a, omega, pi2l):
    #Reduce a list of pi2's with the same pi3; input data precomputed given by canonical_form(..., precompute=True)
    s0 = set()
    pi3t = tuple(pi3a)
    for pi2 in pi2l:
        pi2a = [tuple(sorted(perm(i) for i in t)) for t in pi2]
        l = libgap.SmallestImageSet(Aut1, sorted([omega_d[t] for t in pi2a])).sage()
        s0.add((tuple(omega[i-1] for i in l), pi3t))
    return sorted(s0)


def get_signature(G):
    return (G.index(),G.genus(),G.ncusps(),G.nu2(),G.nu3())


def tuple_perm_notation_to_str_notation(perm):
    """
    Example: ((1, 4), (2, 7)) -> "(1, 4) (2, 7)"
    """
    return str(perm).replace("),", ")").replace(",)",")")[1:][:-1]


def convert_list_to_signature_dict(subgroup_list):
    """
    Sort list of subgroups into dictionary with signatures as keys.
    """
    sig_sorted = {}
    for G in subgroup_list:
        sig = get_signature(G)
        if sig in sig_sorted:
            sig_sorted[sig].append(G)
        else:
            sig_sorted[sig] = [G]
    return sig_sorted


def tuples_to_subgroups(perms):
    """
    Convert list of permutations in tuple notation to ArithSubgroups.
    """
    perms_str = [(tuple_perm_notation_to_str_notation(o2),tuple_perm_notation_to_str_notation(o3)) for (o2,o3) in perms]
    subgroups = [ArithmeticSubgroup_Permutation(S2=o2,S3=o3) for (o2,o3) in perms_str] 
    return subgroups