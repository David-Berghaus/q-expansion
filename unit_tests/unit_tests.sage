load("unit_tests_point_matching_dp.sage")
load("unit_tests_point_matching_arb_wrap.sage")
load("unit_tests_gmres_arb_wrap.sage")
load("unit_tests_ir_arb_wrap.sage")
load("unit_tests_newton_genus_zero.sage")
load("unit_tests_fourier_expansion.sage")
load("unit_tests_eisenstein.sage")
load("unit_tests_belyi_map.sage")
load("unit_tests_passports.sage")

def run_unit_tests():
    run_unit_tests_point_matching_dp()
    run_unit_tests_point_matching_arb_wrap()
    run_unit_tests_gmres_arb_wrap()
    run_unit_tests_ir_arb_wrap()
    run_unit_tests_newton_genus_zero()
    run_unit_tests_fourier_expansion()
    run_unit_tests_eisenstein()
    run_unit_tests_belyi_map()
    run_unit_tests_passports()

run_unit_tests()