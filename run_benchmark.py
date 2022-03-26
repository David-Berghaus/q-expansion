"""
Routines for computing a passport and storing the result.
Note that this file is supposed to be called from the commandline like this:
sage create_database_entry.py <passport_index> <genus> <rigorous_trunc_order> <eisenstein_digit_prec> <max_weight>
"""

from sage.all_cmdline import *   # import sage library

import os
import sys
from pathlib import Path

from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from psage.modform.arithgroup.mysubgroup import MySubgroup

load("benchmarks.sage")

version = int(sys.argv[1])
#S = AutomorphicFormSpace(Gamma0(1),12)
S = AutomorphicFormSpace(MySubgroup(o2='(1)(2 4)(3 7)(5 10)(6 11)(8 14)(9 15)(12 13)(16 17)',o3='(1 7 4)(2 11 10)(3 15 14)(5)(6 12 13)(8 17 9)(16)'),4) #Noncongruence group of signature (17, 0, 3, 1, 2)
digit_prec = 200

print("version: ", version)
print("S: ", S)
print("digit_prec: ", digit_prec)

run_benchmark(S,digit_prec,version)
