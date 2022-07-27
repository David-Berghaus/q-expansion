"""
Routines for loading and working database entries
"""

import os
import sys
from string import ascii_lowercase

def load_database(database_path, index=None, genus=None, Kv_degree=None):
    """
    Loads entries of database and returns them as a list.
    We allow to filter for index, genus and Kv_degree by specifying them as ints or lists of ints.
    """
    if genus != None and genus > 1:
        raise NotImplementedError("This case has not been implemented yet!")
    if isinstance(index, list) == False and index != None: #Index is an int which we convert to a list
        index = [index]
    if isinstance(genus, list) == False and genus != None: #genus is an int which we convert to a list
        genus = [genus]
    if isinstance(Kv_degree, list) == False and Kv_degree != None: #Kv_degree is an int which we convert to a list
        Kv_degree = [Kv_degree]

    res = {}
    passport_list = []
    if genus == None:
        passport_list_g_zero = load("data/genus_zero_passport_list.sobj")
        passport_list_g_one = load("data/genus_one_passport_list.sobj")
        passport_list = passport_list_g_zero+passport_list_g_one
    else:
        if 0 in genus:
            passport_list += load("data/genus_zero_passport_list.sobj")
        if 1 in genus:
            passport_list += load("data/genus_one_passport_list.sobj")
    for i in range(len(passport_list)):
        G = passport_list[i][0]
        if index == None or G.index() in index:
            if genus == None or G.genus() in genus:
                for letter in ascii_lowercase:
                    entry_name = get_signature_to_underscore(G.signature()) + "_" + str(get_signature_pos(i,passport_list)) + "_" + letter
                    storage_path = os.path.join(database_path,str(G.index()),str(G.genus()))
                    file_path = storage_path + "/" + entry_name + ".sobj"
                    if os.path.isfile(file_path) == True:
                        entry = load(file_path)
                        if Kv_degree == None or entry["Kv"].degree() in Kv_degree:
                            res[entry_name] = entry
                    else:
                        break
    return res

def get_signature_pos(passport_index, passport_list):
    """
    In case passport_list contains several signatures of the same time, return the number n such that passport is the n-th passport of given signature.
    """
    passport_signature = passport_list[passport_index][0].signature()
    signatures = [passport_list[i][0].signature() for i in range(len(passport_list))]
    for i in range(len(signatures)):
        if passport_signature == signatures[i]:
            return passport_index-i

def get_signature_to_underscore(signature):
    """
    Given a signature (represented as a tuple), return string with underscores as separators.
    """
    res = ""
    for i in signature:
        res += str(i) + "_"
    res = res[:-1] #Remove last underscore
    return res

def print_missing_passports(genus, database_path, max_passport_index=None):
    """
    Prints all passports up to max_passport_index that have not been computed yet.
    """
    if genus > 1:
        raise NotImplementedError("This case has not been implemented yet!")
    if genus == 0:
        passport_list = load("data/genus_zero_passport_list.sobj")
    elif genus == 1:
        passport_list = load("data/genus_one_passport_list.sobj")
    if max_passport_index == None:
        max_passport_index = len(passport_list)
    for i in range(max_passport_index):
        G = passport_list[i][0]
        entry_name = get_signature_to_underscore(G.signature()) + "_" + str(get_signature_pos(i,passport_list))
        storage_path = os.path.join(database_path,str(G.index()),str(G.genus()))
        if os.path.isfile(storage_path + "/" + entry_name + ".sobj") == False:
            print(i)

def create_list_of_passport_labels(genus, max_passport_index=None):
    """
    Create a list of all passport labels (represented as strings).
    """
    passport_label_list = []
    if genus == 0:
        passport_list = load("data/genus_zero_passport_list.sobj")
    elif genus == 1:
        passport_list = load("data/genus_one_passport_list.sobj")
    else:
        raise NotImplementedError("This case has not been implemented yet!")
    if max_passport_index == None:
        max_passport_index = len(passport_list)
    for i in range(max_passport_index):
        G = passport_list[i][0]
        passport_label = get_signature_to_underscore(G.signature()) + "_" + str(get_signature_pos(i,passport_list))
        passport_label_list.append(passport_label)
    return passport_label_list
