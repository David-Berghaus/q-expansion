from sage.all_cmdline import *   # import sage library

import itertools
import os
import sys
from pathlib import Path
from string import ascii_lowercase

database_path = "/mnt/c/Users/machs/Desktop/noncong database"
#database_path = "/cephfs/user/s6ddberg/noncong_database"

def get_belyi_friends():
    from lmfdb import db
    belyi_galmaps = db.belyi_galmaps

    paths = get_paths_list()
    belyi_friends = {}
    for path in paths:
        label = path_to_label(path)
        with open(path+".sage", "r") as f:
            data = f.readlines()
            with open("tmp.sage", "w") as f:
                f.write("G = " + data[36][11:] + "\n")
                f.write(data[40] + "\n")
                f.write(data[41])
        load("tmp.sage")

        P = G.perm_group()
        mu = int(G.index())
        [a_s, b_s, c_s] = sorted([int(G.S2().order()),int(G.S3().order()),int((G.S2()*G.S3()).order())])
        if K.degree() == 1:
            base_field = list(map(lambda x: int(x),[-1,1]))
        else:
            base_field = list(map(lambda x: int(x),list(K.polynomial())))
        for [a,b,c] in itertools.permutations([a_s, b_s, c_s]): #BelyiDB doesn't have this canonical atm so we have to check all permutations
            for belyi_map in belyi_galmaps.search({'a_s': a, 'b_s': b, 'c_s': c, 'base_field': base_field, 'deg': mu}, projection=['label','triples']):
                for triple in belyi_map['triples']:
                    P_belyi = PermutationGroup([Permutation(triple[2]), Permutation(triple[1])])
                    conj = libgap.RepresentativeAction(gap(f"SymmetricGroup({mu})"), gap(P), gap(P_belyi)) #Check for conjugacy, see: https://ask.sagemath.org/question/44357/determining-if-two-subgroups-of-a-symmetric-group-are-conjugate/?answer=47983#post-id-47983
                    if conj != gap("fail"):
                        belyi_friends[path] = belyi_map['label']
        if path not in belyi_friends:
            belyi_friends[path] = "\\N"
    return belyi_friends

def get_paths_list():
    return load("data/entry_paths.sobj")

def path_to_label(path):
    return path.split("/")[-1][:-5]