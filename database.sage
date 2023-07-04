# ****************************************************************************
#       Copyright (C) 2022 David Berghaus <berghaus@th.physik.uni-bonn.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

"""
Routines for loading and working database entries
"""

import itertools
import os
from os.path import isfile
import re
import requests
import sys
from string import ascii_lowercase

def load_entry(database_path, entry_name, load_floating_expansions=False):
    """
    Given an entry_name (as a string), load data attached to this entry.
    """
    index = ""
    for i in entry_name:
        if i == "T":
            break
        index += i

    #Determine genus (we could in principle also compute it using the formula)
    #Search for all folders in database_path/index and see where the entry exists
    genus = None
    for folder in os.listdir(os.path.join(database_path,index)):
        if isfile(os.path.join(database_path,index,folder,entry_name+".sobj")):
            genus = folder
            break
    if genus == None:
        raise ValueError("Could not find entry in database!")

    storage_path = os.path.join(database_path,index,genus)
    file_path = storage_path + "/" + entry_name + ".sobj"
    passport = load(file_path)
    if load_floating_expansions == False:
        return passport
    else:
        file_path = storage_path + "/" + entry_name + "_floating_expansions.sobj"
        return passport, load(file_path)

def get_entry_name(G): #We have to copy this function here to avoid circular imports
    gap_label = str(G.index()) + "T" + str(libgap.TransitiveIdentification(G.perm_group()))
    return gap_label + "-" + get_lambda(G.permT) + "_" + get_lambda(G.permR) + "_" + get_lambda(G.permS)

def get_lambda(perm):
    cycle_type = perm.cycle_type()
    cycle_type.sort(reverse=True)
    res = ""
    for i in cycle_type:
        res += str(i) + "."
    return res[:-1]

def load_database(database_path, index=None, genus=None, Kv_degree=None, load_floating_expansions=False):
    """
    Loads entries of database and returns them as a dict.
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
    res_fl = {}
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
                entry_name_no_orbit = get_entry_name(G)
                for letter in ascii_lowercase:
                    entry_name = entry_name_no_orbit + "-" + letter
                    storage_path = os.path.join(database_path,str(G.index()),str(G.genus()))
                    file_path = storage_path + "/" + entry_name + ".sobj"
                    if isfile(file_path) == True:
                        entry = load(file_path)
                        if Kv_degree == None or entry["Kv"].degree() in Kv_degree:
                            res[entry_name] = entry
                            if load_floating_expansions == True:
                                file_path = storage_path + "/" + entry_name + "_floating_expansions.sobj"
                                res_fl[entry_name] = load(file_path)
                    else:
                        break
    if load_floating_expansions == False:
        return res
    else:
        return res, res_fl

def data_to_LMFDB_txt(passport_data, label, curves_file="curves.txt", spaces_file="spaces.txt", forms_file="forms.txt", lookup_belyi_friend=True):
    """
    Store passport_data into a txt file that can be easily inserted into the LMFDB.
    """
    curve_to_LMFDB_txt(passport_data, label, file=curves_file, lookup_belyi_friend=lookup_belyi_friend)
    G = passport_data["G"]
    if G.genus() == 0:
        space_to_LMFDB_txt(passport_data,label,0,"H",file=spaces_file)
        form_to_LMFDB_txt(passport_data,label,0,"H",0, file=forms_file)
    for weight in range(2,7,2):
        if G.dimension_modular_forms(weight) > 0:
            space_to_LMFDB_txt(passport_data,label,weight,"M",file=spaces_file)
            for modform_pos in range(G.dimension_modular_forms(weight)):
                form_to_LMFDB_txt(passport_data,label,weight,"M",modform_pos,file=forms_file)
            if G.dimension_cusp_forms(weight) > 0:
                space_to_LMFDB_txt(passport_data,label,weight,"C",file=spaces_file)
                for modform_pos in range(G.dimension_cusp_forms(weight)):
                    form_to_LMFDB_txt(passport_data,label,weight,"C",modform_pos,file=forms_file)

def curve_to_LMFDB_txt(passport_data, label, file="curves.txt", lookup_belyi_friend=True):
    """
    Given passport_data (as a dict), create a txt for the curve that can be easily inserted into the LMFDB.
    """
    if not os.path.exists(file):
        header1 = "genus|psl2z_index|n_c|n_e2|n_e3|cusp_widths|permS|permR|permT|label|monodromy_group|K|u_str|v_L|belyi_map|elliptic_curve|hyperelliptic_curve|friends|passport_reps|L|passport_embeddings|is_congruence"
        header2 = "integer|integer|integer|integer|integer|integer[]|text|text|text|text|text|text|text|text|text|text|text|text[]|text[]|numeric[]|double precision[]|boolean"
        with open(file, "a") as f:
            f.write(header1 + "\n")
            f.write(header2 + "\n")
            f.write("\n")

    G = passport_data["G"]
    with open(file, "a") as f:
        f.write(str(G.genus()) + "|")
        f.write(str(G.index()) + "|")
        f.write(str(G.ncusps()) + "|")
        f.write(str(G.nu2()) + "|")
        f.write(str(G.nu3()) + "|")
        f.write("{" + str(G.cusp_widths())[1:-1] + "}|")
        f.write(str(G.S2()) + "|")
        f.write(str(G.S3()) + "|")
        f.write(str(G.S2()*G.S3()) + "|")
        f.write(label + "|")
        monodromy_group_label = get_monodromy_group_label(G)
        f.write(monodromy_group_label + "|")
        f.write(str(get_number_field_LMFDB_label(passport_data["K"])) + "|")
        f.write(str(passport_data["u_str"]) + "|")
        f.write(str(passport_data["v"]) + "|")
        if G.genus() == 0:
            f.write(str(passport_data["curve"]) + "|") #To Do: Print in pretty form
            f.write('\\N' + "|")
            f.write('\\N' + "|")
        elif G.genus() == 1:
            f.write('\\N' + "|")
            f.write(str(passport_data["curve"]) + "|")
            f.write('\\N' + "|")
        else:
            f.write('\\N' + "|")
            f.write('\\N' + "|")
            f.write(str(passport_data["curve"]) + "|")
        friends = []
        belyi_friend = get_belyi_friend(passport_data,label,lookup=lookup_belyi_friend)
        if belyi_friend != "\\N":
            friends.append(belyi_friend)
        if G.genus() == 1:
            ec_friend = get_elliptic_curve_friend(passport_data["curve"],passport_data["K"])
            if ec_friend != "\\N":
                friends.append(ec_friend)
        f.write("{" + str(friends)[1:-1] + "}|")
        f.write("{{" + ','.join(["{" + str(embedding)[1:-1] + "}" for embedding in passport_data["embeddings"].keys()])[1:-1] + "}}|") #We use ".join" in order not to print the quotes
        f.write("{" + str(list(passport_data["L"].polynomial()))[1:-1] + "}|")
        f.write("{" + ','.join(["{" + str(complex_number_to_doubles(embedding))[1:-1] + "}" for embedding in passport_data["embeddings"].values()])[1:-1] + "}|")
        f.write("t" if passport_data["is_congruence"] else "f")
        f.write("\n")
        
def get_monodromy_group_label(G):
    i = libgap.TransitiveIdentification(G.perm_group())
    return str(G.index()) + "T" + str(i)

def get_number_field_LMFDB_label(K):
    """
    Call the lmfdb API to get the label of the number field K.
    """
    if K.degree() == 1:
        return "1.1.1.1"
    query = "https://www.lmfdb.org/api/nf_fields/?coeffs=li"
    for coeff in list(K.polynomial()):
        query += str(coeff) + ","
    query = query[:-1] + "&_format=json"

    response = requests.get(query)
    if response.status_code == 200:
        data = response.json()
        try:
            return data["data"][0]["label"]
        except:
            print("Number field not found in the LMFDB!")
    return list(K.polynomial()) #If the number field is not in the LMFDB, we return the polynomial

def get_elliptic_curve_friend(E, K):
    if K.degree() != 1:
        return "\\N" #We only search for elliptic curves over Q for now
    query = "https://www.lmfdb.org/api/ec_curvedata/?ainvs=li"
    for a_inv in E.a_invariants():
        query += str(a_inv) + ","
    query = query[:-1] + "&_format=json"

    response = requests.get(query)
    if response.status_code == 200:
        data = response.json()
        try:
            ec_friend = data["data"][0]["lmfdb_label"]
            return "EllipticCurve/Q/" + ec_friend
        except:
            print("Elliptic curve not found in the LMFDB!")
    return "\\N" #If the elliptic curve is not in the LMFDB, we return NULL

#Dont forget to run 'export PYTHONPATH=${PYTHONPATH}:<path to lmfdb>' in the terminal in advance
def get_belyi_friend(passport_data, label, lookup=True):
    if not lookup: #Try to find the Belyi friend in a precomputed database, for example if a local copy of the LMFDB is not installed
        try:
            d = load("data/belyi_friends.sobj")
            return d[label]
        except:
            return "\\N"

    from lmfdb import db
    belyi_galmaps = db.belyi_galmaps

    G = passport_data["G"]
    P = G.perm_group()
    mu = int(G.index())
    [a_s, b_s, c_s] = sorted([int(G.S2().order()),int(G.S3().order()),int((G.S2()*G.S3()).order())])
    if K.degree() == 1:
        base_field = list(map(lambda x: int(x),[-1,1]))
    else:
        base_field = list(map(lambda x: int(x),list(passport_data["K"].polynomial())))
    for [a,b,c] in itertools.permutations([a_s, b_s, c_s]): #BelyiDB doesn't have this canonical atm so we have to check all permutations
        for belyi_map in belyi_galmaps.search({'a_s': a, 'b_s': b, 'c_s': c, 'base_field': base_field, 'deg': mu}, projection=['label','triples']):
            for triple in belyi_map['triples']:
                P_belyi = PermutationGroup([Permutation(triple[2]), Permutation(triple[1])])
                conj = libgap.RepresentativeAction(gap(f"SymmetricGroup({mu})"), gap(P), gap(P_belyi)) #Check for conjugacy, see: https://ask.sagemath.org/question/44357/determining-if-two-subgroups-of-a-symmetric-group-are-conjugate/?answer=47983#post-id-47983
                if conj != gap("fail"):
                    belyi_friend_label = belyi_map['label']
                    return "Belyi/" + belyi_friend_label
    return "\\N" 

def complex_number_to_doubles(z):
    return [float(z.real()), float(z.imag())]

def space_to_LMFDB_txt(passport_data, label, weight, modform_type, file="spaces.txt"):
    """
    Given passport_data (as a dict), create a txt for the mf space that can be easily inserted into the LMFDB.
    """
    if not os.path.exists(file):
        header1 = "dim|weight|label|mf_curve"
        header2 = "integer|integer|text|text"
        with open(file, "a") as f:
            f.write(header1 + "\n")
            f.write(header2 + "\n")
            f.write("\n")
    if modform_type == "M":
        mf_space = passport_data["q_expansions"][weight]["modforms_pretty"]
    elif modform_type == "C":
        mf_space = passport_data["q_expansions"][weight]["cuspforms_pretty"]
    else:
        mf_space = [passport_data["q_expansions"][weight]["hauptmodul_pretty"]]

    with open(file, "a") as f:
        f.write(str(len(mf_space)) + "|")
        f.write(str(weight) + "|")
        f.write(label + "." + str(weight) + "." + modform_type + "|")
        f.write(label)
        f.write("\n")

def form_to_LMFDB_txt(passport_data, label, weight, modform_type, modform_pos, file="forms.txt"):
    """
    Create a txt for the form f_{modform_pos} that can be easily inserted into the LMFDB.
    """
    if not os.path.exists(file):
        header1 = "weight|cusp_width|valuation|label|K|u_str|v_L|mf_space|u|v|L|coefficient_numerators|coefficient_denominators"
        header2 = "integer|integer|integer|text|text|text|text|text|double precision[]|double precision[]|numeric[]|numeric[]|numeric[]"
        with open(file, "a") as f:
            f.write(header1 + "\n")
            f.write(header2 + "\n")
            f.write("\n")
    
    if modform_type == "M":
        form = passport_data["q_expansions"][weight]["modforms_pretty"][modform_pos]
    elif modform_type == "C":
        form = passport_data["q_expansions"][weight]["cuspforms_pretty"][modform_pos]
    else:
        form = passport_data["q_expansions"][weight]["hauptmodul_pretty"]
    
    with open(file, "a") as f:
        f.write(str(weight) + "|")
        f.write(str(passport_data["G"].cusp_width(Cusp(1,0))) + "|")
        f.write(str(form.valuation()) + "|")
        f.write(label + "." + str(weight) + "." + modform_type + "." + ascii_lowercase[modform_pos] + "|")
        f.write(get_number_field_LMFDB_label(passport_data["K"]) + "|")
        f.write(passport_data["u_str"] + "|")
        f.write(str(passport_data["v"]) + "|")
        f.write(label + "." + str(weight) + "." + modform_type + "|")
        f.write("{" + str(complex_number_to_doubles(passport_data["u"]))[1:-1] + "}|")
        f.write("{" + str(complex_number_to_doubles(passport_data["v"]))[1:-1] + "}|")
        f.write("{" + str(list(passport_data["L"].polynomial()))[1:-1] + "}|")
        numerators, denominators = get_numerators_and_denominators(form)
        f.write("{{" + ','.join(["{" + str(numerator)[1:-1] + "}" for numerator in numerators])[1:-1] + "}}|")
        f.write("{{" + ','.join(["{" + str(denominators)[1:-1] + "}" for denominator in denominators])[1:-1] + "}}|")
        f.write("\n")

def get_numerators_and_denominators_of_expression_in_K(x):
    """
    Given an expression x in a number field K, transform the result to a common denominator and return the numerators as a list.
    """
    denominator = 1
    for i in range(x.polynomial().degree()+1):
        denominator *= x[i].denominator()
        x *= x[i].denominator()
    return denominator, list(x)

def get_numerators_and_denominators(form):
    numerators, denominators = [], []
    for (i,coeff) in enumerate(list(form)):
        if coeff[i] == 0:
            coeff_denominator, coeff_numerators = get_numerators_and_denominators_of_expression_in_K(coeff[0]) #For u=1, we only have powers of u^0
        else:
            coeff_denominator, coeff_numerators = get_numerators_and_denominators_of_expression_in_K(coeff[i]) #We need to take the i-th index to get rid of "u"
        numerators.append(coeff_numerators)
        denominators.append(coeff_denominator)
    return numerators, denominators

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

def print_pending_passports(genus, database_path, max_passport_index=None):
    """
    Print indices of passports whose computation is currently ongoing (i.e., for which a state or unresolved file exists).
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
        passport_name = get_signature_to_underscore(G.signature()) + "_" + str(get_signature_pos(i,passport_list))
        storage_path = os.path.join(database_path,str(G.index()),str(G.genus()))
        unresolved_passport_file_path = storage_path + "/" + passport_name + "_unresolved_passport.sobj"
        if isfile(unresolved_passport_file_path):
            print(i)
            continue
        for letter in ascii_lowercase:
            state_file_path = storage_path + "/" + passport_name + "_" + letter + "_state.sobj"
            if isfile(state_file_path):
                print(i)
                break

def print_K_polynomials_latex(db_path, only_noncongruence_subgroups=False):
    """
    Prints the polynomials of all number fields in the database in a latex table.
    """
    for index in range(2,18):
        curr_path = db_path + "/" + str(index)
        curr_els = []
        for genus in [0,1]:
            if os.path.isdir(curr_path + f"/{genus}"):
                #Find all files ending with .sage
                for file in os.listdir(curr_path + f"/{genus}"):
                    if file.endswith(".sage"):
                        #Open the file
                        with open(curr_path + f"/{genus}/" + file, "r") as f:
                            #Read the file
                            data = f.readlines()
                            if only_noncongruence_subgroups and data[37][-3] == "u": #True -> congruence subgroup
                                continue
                            with open("dummy.sage", "w") as f2: #Very ugly but sage_eval(x = ...) does not work
                                f2.write("P.<v> = PolynomialRing(QQ)\n")
                                f2.write(data[41])
                            load("dummy.sage")
                            if K.degree() == 1:
                                K_latex = "\\QQ"
                            else:
                                K_latex = latex(K.polynomial())
                            #Get file name which is the name after the last / before .sage
                            name = file.split("/")[-1][:-5]
                            name = name.replace("_","\_")
                            curr_els.append((name, K_latex))
        curr_els.sort(key=lambda x: x[0])
        print("")
        print("\\section{$\\mu$ = " + str(index) + "}")
        print("\\begin{longtable}{|c|c|}") 
        print("\\hline")
        for (name, K_latex) in curr_els:
            print("\\emph{" + name + "} & $ " + K_latex + "$ \\\\")
        print("\\hline")
        print("\\end{longtable}")

def get_paths_of_all_labels(db_path, only_noncongruence_subgroups=False):
    """
    Get list of paths of all labels in the database.
    """
    paths = []
    for index in range(2,18):
        curr_path = db_path + "/" + str(index)
        for genus in [0,1]:
            if os.path.isdir(curr_path + f"/{genus}"):
                #Find all files ending with .sage
                for file in os.listdir(curr_path + f"/{genus}"):
                    if file.endswith(".sage"):
                        #Open the file
                        with open(curr_path + f"/{genus}/" + file, "r") as f:
                            #Read the file
                            if only_noncongruence_subgroups:
                                data = f.readlines()
                                if data[37][-3] == "u": #True -> congruence subgroup
                                    continue
                            path = curr_path + f"/{genus}/" + file
                            paths.append(path[:-5]) #Remove .sage
    return paths

def print_missing_passports(genus, database_path, max_orbit_size, max_passport_index=None):
    """
    Print indicies of all passports up to max_passport_index that have not been computed yet (and are also not pending).
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
        if len(passport_list[i]) > max_orbit_size:
            continue
        G = passport_list[i][0]
        entry_name = get_signature_to_underscore(G.signature()) + "_" + str(get_signature_pos(i,passport_list))
        storage_path = os.path.join(database_path,str(G.index()),str(G.genus()))
        pattern = re.compile(r'{}.*\.sobj$'.format(entry_name)) #Basically entry_name*.sobj
        found_file = False
        for filepath in os.listdir(storage_path):
            if pattern.match(filepath):
                found_file = True #Found either a finished or ongoing computation
                break
        if found_file == False:
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
