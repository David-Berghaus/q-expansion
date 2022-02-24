"""
Routines for computing a passport and storing the result.
Note that this file is supposed to be called from the commandline like this:
sage create_database_entry.py <passport_index> <genus> <rigorous_trunc_order> <eisenstein_digit_prec> <max_weight>
"""

from sage.all_cmdline import *   # import sage library

import os
import sys
from pathlib import Path

load("subgroups.sage")
load("passports.sage")
load("database.sage")

# database_path = "/mnt/c/Users/machs/Desktop/database_test"
database_path = "/cephfs/user/s6ddberg/noncong_database"

def compute_database_entry(passport_index, genus, rigorous_trunc_order, eisenstein_digit_prec, max_weight):
    if genus == 0:
        passport_list = load("data/genus_zero_passport_list.sobj")
        passport = passport_list[passport_index]
        G = passport[0]
        entry_name = get_signature_to_underscore(G.signature()) + "_" + str(get_signature_pos(passport_index,passport_list))
        res, newton_res = compute_passport_data_genus_zero(passport,rigorous_trunc_order,eisenstein_digit_prec,max_weight,return_newton_res=True)
        storage_path = os.path.join(database_path,str(G.index()),str(G.genus()))
        Path(storage_path).mkdir(parents=True, exist_ok=True) #Check if path exists, if not, create it
        save(res,os.path.join(storage_path,entry_name))
        save(newton_res,os.path.join(storage_path,entry_name+"_newton_res"))
    elif genus == 1: #We could use less copy-paste here...
        passport_list = load("data/genus_one_passport_list.sobj")
        passport = passport_list[passport_index]
        G = passport[0]
        entry_name = get_signature_to_underscore(G.signature()) + "_" + str(get_signature_pos(passport_index,passport_list))
        res = compute_passport_data_higher_genera(passport,rigorous_trunc_order,eisenstein_digit_prec,max_weight)
        storage_path = os.path.join(database_path,str(G.index()),str(G.genus()))
        Path(storage_path).mkdir(parents=True, exist_ok=True) #Check if path exists, if not, create it
        save(res,os.path.join(storage_path,entry_name))
    else:
        raise NotImplementedError("This case has not been implemented yet!")

passport_index = int(sys.argv[1])
genus = int(sys.argv[2])
rigorous_trunc_order = int(sys.argv[3])
eisenstein_digit_prec = int(sys.argv[4])
max_weight = int(sys.argv[5])
compute_database_entry(passport_index,genus,rigorous_trunc_order,eisenstein_digit_prec,max_weight)
