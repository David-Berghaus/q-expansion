"""
Routines for computing a passport and storing the result.
Note that this file is supposed to be called from the commandline like this:
sage create_database_entry.py <passport_index> <genus> <eisenstein_digit_prec> <max_weight>
"""

from sage.all_cmdline import *   # import sage library

import os
import sys
from pathlib import Path
from string import ascii_lowercase

load("subgroups.sage")
load("passports.sage")
load("database.sage")

# database_path = "/mnt/c/Users/machs/Desktop/database_test"
database_path = "/cephfs/user/s6ddberg/noncong_database"

def get_passport_list(genus):
    if genus == 0:
        return load("data/genus_zero_passport_list.sobj")
    elif genus == 1:
        return load("data/genus_one_passport_list.sobj")
    else:
        raise NotImplementedError("This case has not been implemented yet!")

def compute_database_entry(passport_index, genus, eisenstein_digit_prec, max_weight, rigorous_trunc_order=None):
    print(f"Computing passport of index {passport_index} and genus {genus} to {eisenstein_digit_prec} digits precision up to weight {max_weight}.")
    passport_list = get_passport_list(genus)
    passport = passport_list[passport_index]
    G = passport[0]
    storage_path = os.path.join(database_path,str(G.index()),str(G.genus()))
    Path(storage_path).mkdir(parents=True, exist_ok=True) #Check if path exists, if not, create it
    if genus == 0 and rigorous_trunc_order == None:
        rigorous_trunc_order = int(2000/(len(passport)+G.cusp_width(Cusp(1,0)))+50) #Heuristic choice that runs in reasonable amount of CPU time
    entry_name = get_signature_to_underscore(G.signature()) + "_" + str(get_signature_pos(passport_index,passport_list))
    unresolved_passport_path = os.path.join(storage_path,entry_name+"_unresolved_passport.sobj") #Path were we store unresolved passport elements in case the passport decays
    if os.path.exists(unresolved_passport_path) == True:
        passport = load(unresolved_passport_path)
        G = passport[0]
        for letter in ascii_lowercase[1:]:
            if os.path.exists(os.path.join(storage_path,entry_name+"_"+letter+".sobj")) == False:
                entry_name += "_"+letter
                break
    else:
        entry_name += "_a"
        if os.path.exists(os.path.join(storage_path,entry_name+".sobj")) == True: #Computation is finished
            return
    state_file_path = os.path.join(storage_path,entry_name+"_state.sobj") #Path where we store intermediate results to restart the computation
    if genus == 0:
        res, floating_expansions = compute_passport_data_genus_zero(passport,rigorous_trunc_order,eisenstein_digit_prec,max_weight,state_file_path=state_file_path)
    elif genus == 1:
        res, floating_expansions = compute_passport_data_higher_genera(passport,rigorous_trunc_order,eisenstein_digit_prec,max_weight,state_file_path=state_file_path)
    unresolved_passport_elements = get_unresolved_passport_elements(passport,res["embeddings"])
    if len(unresolved_passport_elements) != 0: #We still need to compute the remaining passport elements!
        if G.genus() == 1:
            unresolved_passport_elements = optimize_passport_for_period_integral(unresolved_passport_elements) #Update for better period integral evaluation
        save(unresolved_passport_elements,unresolved_passport_path)
    elif os.path.exists(unresolved_passport_path) == True:
        os.remove(unresolved_passport_path)
    save(res,os.path.join(storage_path,entry_name))
    save(floating_expansions,os.path.join(storage_path,entry_name+"_floating_expansions"))

passport_index = int(sys.argv[1])
genus = int(sys.argv[2])
eisenstein_digit_prec = int(sys.argv[3])
max_weight = int(sys.argv[4])
compute_database_entry(passport_index,genus,eisenstein_digit_prec,max_weight)
