load("../passports.sage")

def run_unit_tests_passports():
    test_passport()
    
def test_passport():
    passport = [MySubgroup(o2='(1 6)(2)(3 4)(5 7)',o3='(1 7 6)(2 3 5)(4)'),MySubgroup(o2='(1 4)(2)(3 5)(6 7)',o3='(1 5 4)(2 3 6)(7)')]
    res = compute_passport_data_genus_zero(passport,10,50,4)
    print("test_passport ok")
    #compute_passport_data_genus_zero already tests itself by comparing its result to numerical values so we don't need assert checks