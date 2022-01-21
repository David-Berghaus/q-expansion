from ast import copy_location
from psage.modform.arithgroup.mysubgroup import MySubgroup

from classes.belyi_map import BelyiMap

def run_unit_tests_belyi_map():
    test_belyi_map()
    test_belyi_map2()
    test_belyi_map3()
    test_belyi_map4()
    
def test_belyi_map():
    U1 = MySubgroup(o2='(1 2)(3 4)(5)(6 7)',o3='(1)(2 3 5)(4 6 7)')
    B = BelyiMap(U1) #Test recognition of numberfield that is not QQ
    tmp_rig = B.get_hauptmodul_q_expansion(10,only_principal_cusp_expansion=False)
    tmp_approx = B.get_hauptmodul_q_expansion_approx(25,50,only_principal_cusp_expansion=False)
    CF = ComplexField(tmp_approx.get_cusp_expansion(Cusp(1,0)).base_ring().precision()) #We work with CF here because otherwise the differences might be empty balls
    #Test q-expansions at principal cusp
    assert (tmp_rig.get_cusp_expansion(Cusp(1,0))[9]-CF(tmp_approx.get_cusp_expansion(Cusp(1,0))[9])).abs() < CF(10)**(-42)
    assert (CF(-11.18048794365418621081119907300687097510914273497,74.25239079134411962653170377372871107862679259671)-CF(tmp_approx.get_cusp_expansion(Cusp(1,0))[9])).abs() < CF(10)**(-42)
    #Test q-expansion at other cusp
    assert (tmp_rig.get_cusp_expansion(Cusp(2,1))[9]-CF(tmp_approx.get_cusp_expansion(Cusp(2,1))[9])).abs() < CF(10)**(-44)
    assert (CF(-204788.266511846997417732376120597452504942099776674366809059697369771255508569807611107026928336726,-354856.283075376626293839386472132772905098041315141845293562059967286426080264752766926391137344234)-CF(tmp_approx.get_cusp_expansion(Cusp(2,1))[9])).abs() < CF(10)**(-45)
    print("test_belyi_map ok")

def test_belyi_map2():
    B = BelyiMap(Gamma0(6)) #Test some higher dimensional S and M space computation
    C_rig = B.get_cuspforms(8,25,only_principal_cusp_expansion=False)
    C_approx = B.get_cuspforms(8,25,50,only_principal_cusp_expansion=False) #We also test that a higher precision is used for the series reversion
    C_sage = CuspForms(6,8).q_echelon_basis(25)
    CF = ComplexField(C_approx[0].get_cusp_expansion(Cusp(1,0)).base_ring().precision()) #We work with CF here because otherwise the differences might be empty balls
    assert (C_sage[0][24]-C_rig[0].get_cusp_expansion(Cusp(1,0))[24]).abs() == 0
    assert (C_sage[4][24]-C_rig[4].get_cusp_expansion(Cusp(1,0))[24]).abs() == 0
    assert (C_sage[0][24]-CF(C_approx[0].get_cusp_expansion(Cusp(1,0))[24])).abs() == 0
    assert (C_sage[4][24]-CF(C_approx[4].get_cusp_expansion(Cusp(1,0))[24])).abs() == 0
    #Now to another cusp
    assert (0.0001193725232434080170705685108977290047248894985521-C_rig[0].get_cusp_expansion(Cusp(0,1))[1]).abs() < 1e-44
    assert (0.004483167581161408321902149062642889803383630544025-C_rig[1].get_cusp_expansion(Cusp(0,1))[4]).abs() < 1e-44
    assert (0.0001193725232434080170705685108977290047248894985521-CF(C_approx[0].get_cusp_expansion(Cusp(0,1))[1])).abs() < 1e-44
    assert (0.004483167581161408321902149062642889803383630544025-CF(C_approx[1].get_cusp_expansion(Cusp(0,1))[4])).abs() < 1e-44
    M_rig = B.get_modforms(2,25,only_principal_cusp_expansion=False)
    M_approx = B.get_modforms(2,25,50,only_principal_cusp_expansion=False) #We also test that a higher precision is used for the series reversion
    M_sage = ModularForms(6,2).q_echelon_basis(25)
    assert (M_sage[0][24]-M_rig[0].get_cusp_expansion(Cusp(1,0))[24]).abs() == 0
    assert (M_sage[2][24]-M_rig[2].get_cusp_expansion(Cusp(1,0))[24]).abs() == 0
    assert (M_sage[0][24]-CF(M_approx[0].get_cusp_expansion(Cusp(1,0))[24])).abs() == 0
    assert (M_sage[2][24]-CF(M_approx[2].get_cusp_expansion(Cusp(1,0))[24])).abs() == 0
    print("test_belyi_map2 ok")

def test_belyi_map3():
    U1 = MySubgroup(o2='(1 2)(3 4)(5)(6 7)',o3='(1)(2 3 5)(4 6 7)')
    B = BelyiMap(U1) #Test degree two Kw example
    u_alg = QQbar(B._Kw.gen())
    v_alg = QQbar(B._Kv.gen())
    C_rig = B.get_cuspforms(4,25,only_principal_cusp_expansion=True) #This is constructed through Kw-arithmetic
    C_approx = B.get_cuspforms(4,25,50,only_principal_cusp_expansion=True)
    CF = ComplexField(167)
    correct_res = CF(-4.2933873728750574663589574362627879379968021749049681273915,-64.175078941013332742794027720411935007114227517392815852039)
    assert (CF(C_rig[0].get_cusp_expansion(Cusp(1,0))[14])-correct_res).abs() < CF(10)**(-39)
    assert (CF(C_approx[0].get_cusp_expansion(Cusp(1,0))[14])-correct_res).abs() < CF(10)**(-39)
    assert QQbar(C_rig[0].get_cusp_expansion(Cusp(1,0))[14])-C_rig[0].get_cusp_expansion(Cusp(1,0),factor_into_u_v=True)[14](u=u_alg,v=v_alg) == 0 #Test factorization into u and v
    C_rig = B.get_modforms(4,25,only_principal_cusp_expansion=True) #This is constructed through Kw-arithmetic
    C_approx = B.get_modforms(4,25,50,only_principal_cusp_expansion=True)
    correct_res = CF(-72.433016538957003294780826994583182363158868739244676993287,-224.75965449478140273855502676201462038401186709430507383978)
    assert (CF(C_rig[1].get_cusp_expansion(Cusp(1,0))[11])-correct_res).abs() < CF(10)**(-39)
    assert (CF(C_approx[1].get_cusp_expansion(Cusp(1,0))[11])-correct_res).abs() < CF(10)**(-39)
    assert QQbar(C_rig[1].get_cusp_expansion(Cusp(1,0))[11])-C_rig[1].get_cusp_expansion(Cusp(1,0),factor_into_u_v=True)[11](u=u_alg,v=v_alg) == 0 #Test factorization into u and v

def test_belyi_map4():
    """
    Unittest that tests the correct recognition and arithmetic of u when the cusp at infinity has width one.
    """
    U7 = MySubgroup(o2='(7 2)(3 4)(5)(6 1)',o3='(7)(2 3 5)(4 6 1)') #Group with cusp of width one at infinity
    B = BelyiMap(U7) #Degree two numberfield
    #Test q-expansions at principal cusp
    tmp_rig = B.get_hauptmodul_q_expansion(10,only_principal_cusp_expansion=True)
    tmp_approx = B.get_hauptmodul_q_expansion_approx(25,50,only_principal_cusp_expansion=True)
    CF = ComplexField(tmp_approx.get_cusp_expansion(Cusp(1,0)).base_ring().precision()) #We work with CF here because otherwise the differences might be empty balls
    assert (tmp_rig.get_cusp_expansion(Cusp(1,0))[9]-CF(tmp_approx.get_cusp_expansion(Cusp(1,0))[9])).abs() < CF(10)**(-42)
    assert (CF(97.97658170916759497398537044515919396276246573035,-7.354816147385987288161433019146684584608985406791)-CF(tmp_approx.get_cusp_expansion(Cusp(1,0))[9])).abs() < CF(10)**(-31)
    #Also test that u-v-factorization works correctly here
    u = B._u_interior_Kv
    assert (CF(97.97658170916759497398537044515919396276246573035,-7.354816147385987288161433019146684584608985406791)-CF(tmp_rig.get_cusp_expansion(Cusp(1,0),factor_into_u_v=True)[9](u=u))).abs() < CF(10)**(-31)
    #Test q-expansions at other cusp
    tmp_rig = B.get_hauptmodul_q_expansion(15,only_principal_cusp_expansion=False) #We recompute this here to make sure that both cases work
    tmp_approx = B.get_hauptmodul_q_expansion_approx(25,50,only_principal_cusp_expansion=False)
    CF = ComplexField(tmp_approx.get_cusp_expansion(Cusp(1,0)).base_ring().precision()) #We work with CF here because otherwise the differences might be empty balls
    assert (tmp_rig.get_cusp_expansion(Cusp(0,1))[9]-CF(tmp_approx.get_cusp_expansion(Cusp(0,1))[9])).abs() < CF(10)**(-38)
    assert (CF(-409614.8283118484107493257393238180260932476735580,125.7860863354206940660961873374166208245937784693)-CF(tmp_approx.get_cusp_expansion(Cusp(0,1))[9])).abs() < CF(10)**(-33)
    print("test_belyi_map4 ok")