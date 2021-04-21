load("unit_tests_point_matching_dp.sage")
load("unit_tests_point_matching_arb_wrap.sage")
load("unit_tests_gmres_arb_wrap.sage")
load("unit_tests_ir_arb_wrap.sage")

def run_unit_tests():
    run_unit_tests_point_matching_dp()
    run_unit_tests_point_matching_arb_wrap()
    run_unit_tests_gmres_arb_wrap()
    run_unit_tests_ir_arb_wrap()