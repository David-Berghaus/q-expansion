load("../passports.sage")

def run_unit_tests_passports():
    test_passport_genus_zero()
    test_passport_higher_genera()
    
def test_passport_genus_zero():
    passport = [MySubgroup(o2='(1 6)(2)(3 4)(5 7)',o3='(1 7 6)(2 3 5)(4)'),MySubgroup(o2='(1 4)(2)(3 5)(6 7)',o3='(1 5 4)(2 3 6)(7)')]
    res = compute_passport_data_genus_zero(passport,10,50,4)
    print("test_passport_genus_zero ok")
    #compute_passport_data_genus_zero already tests itself by comparing its result to numerical values so we don't need assert checks

def test_passport_higher_genera():
    passport = [MySubgroup(o2='(1)(2 5)(3 7)(4 8)(6 9)',o3='(1 2 6)(3 8 5)(4 9 7)')]
    res = compute_passport_data_higher_genera(passport,10,50,6)
    print("test_passport_higher_genera ok")
    #compute_passport_data_higher_genera already tests itself by comparing its result to numerical values so we don't need assert checks