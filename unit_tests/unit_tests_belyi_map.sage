from psage.modform.arithgroup.mysubgroup import MySubgroup

from classes.belyi_map import BelyiMap

def run_unit_tests_belyi_map():
    test_belyi_map()
    test_belyi_map2()
    
def test_belyi_map():
    U1 = MySubgroup(o2='(1 2)(3 4)(5)(6 7)',o3='(1)(2 3 5)(4 6 7)')
    B = BelyiMap(U1) #Test recognition of numberfield that is not QQ
    #Test q-expansions at principal cusp
    tmp_rig = B.get_hauptmodul_q_expansion(10)
    tmp_approx = B.get_hauptmodul_q_expansion_approx(25,50)
    CF = tmp_approx.get_cusp_expansion(Cusp(1,0))[0].parent()
    assert (tmp_rig.get_cusp_expansion(Cusp(1,0))[9]-tmp_approx.get_cusp_expansion(Cusp(1,0))[9]).abs() < CF(10)**(-42)
    assert (CF(-11.18048794365418621081119907300687097510914273497,74.25239079134411962653170377372871107862679259671)-tmp_approx.get_cusp_expansion(Cusp(1,0))[9]).abs() < CF(10)**(-42)
    #Test q-expansion at other cusp
    assert (tmp_rig.get_cusp_expansion(Cusp(2,1))[9]-tmp_approx.get_cusp_expansion(Cusp(2,1))[9]).abs() < CF(10)**(-44)
    assert (CF(-1.687158782764540195465087890625000000000000000000e9,-2.922237035627901077270507812500000000000000000000e9)-tmp_approx.get_cusp_expansion(Cusp(2,1))[24]).abs() < CF(10)**(-45)
    print("test_belyi_map ok")

def test_belyi_map2():
    B = BelyiMap(Gamma0(6)) #Test some higher dimensional S and M space computation
    C_rig = B.get_cuspforms(8,25)
    C_approx = B.get_cuspforms(8,25,50) #We also test that a higher precision is used for the series reversion
    C_sage = CuspForms(6,8).q_echelon_basis(25)
    assert (C_sage[0][24]-C_rig[0].get_cusp_expansion(Cusp(1,0))[24]).abs() == 0
    assert (C_sage[4][24]-C_rig[4].get_cusp_expansion(Cusp(1,0))[24]).abs() == 0
    assert (C_sage[0][24]-C_approx[0].get_cusp_expansion(Cusp(1,0))[24]).abs() == 0
    assert (C_sage[4][24]-C_approx[4].get_cusp_expansion(Cusp(1,0))[24]).abs() == 0
    #Now to another cusp
    assert (0.0001193725232434080170705685108977290047248894985521-C_rig[0].get_cusp_expansion(Cusp(0,1))[1]).abs() < 1e-44
    assert (0.004483167581161408321902149062642889803383630544025-C_rig[1].get_cusp_expansion(Cusp(0,1))[4]).abs() < 1e-44
    assert (0.0001193725232434080170705685108977290047248894985521-C_approx[0].get_cusp_expansion(Cusp(0,1))[1]).abs() < 1e-44
    assert (0.004483167581161408321902149062642889803383630544025-C_approx[1].get_cusp_expansion(Cusp(0,1))[4]).abs() < 1e-44
    M_rig = B.get_modforms(2,25)
    M_approx = B.get_modforms(2,25,50) #We also test that a higher precision is used for the series reversion
    M_sage = ModularForms(6,2).q_echelon_basis(25)
    assert (M_sage[0][24]-M_rig[0].get_cusp_expansion(Cusp(1,0))[24]).abs() == 0
    assert (M_sage[2][24]-M_rig[2].get_cusp_expansion(Cusp(1,0))[24]).abs() == 0
    assert (M_sage[0][24]-M_approx[0].get_cusp_expansion(Cusp(1,0))[24]).abs() == 0
    assert (M_sage[2][24]-M_approx[2].get_cusp_expansion(Cusp(1,0))[24]).abs() == 0
    print("test_belyi_map2 ok")
