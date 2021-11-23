from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
from psage.groups.permutation_alg import MyPermutation
from psage.modform.arithgroup.mysubgroup import MySubgroup

from classes.fourier_expansion import get_cuspform_q_expansion_approx, get_modform_q_expansion_approx, get_cuspform_basis_approx, get_modform_basis_approx

def run_unit_tests_fourier_expansion():
    test_fourier_expansion_cuspform()
    test_fourier_expansion_cuspform_non_fft_non_horner()
    test_fourier_expansion_cuspform_prec_loss()
    test_fourier_expansion_modform()
    test_fourier_expansion_basis_cuspform()
    test_fourier_expansion_basis_modform()

def test_fourier_expansion_cuspform():
    S = AutomorphicFormSpace(Gamma0(11),weight=4) #Search for a multiplicity two new-form
    MF = get_cuspform_q_expansion_approx(S,50,label=1)
    assert abs(MF.get_cusp_expansion(Cusp(1,0),trunc_order=10)[3]-(-4)) < 1e-45
    print("test_fourier_expansion_cuspform ok")

def test_fourier_expansion_cuspform_non_fft_non_horner():
    """
    Test case where 'use_FFT=False' and 'use_Horner=False' and classical matrix-vector multiplication is used.
    """
    S = AutomorphicFormSpace(Gamma0(11),weight=4) #Search for a multiplicity two new-form
    MF = get_cuspform_q_expansion_approx(S,50,label=1,use_FFT=False,use_Horner=False)
    assert abs(MF.get_cusp_expansion(Cusp(1,0),trunc_order=10)[3]-(-4)) < 1e-45
    print("test_fourier_expansion_cuspform_non_fft_non_horner ok")

def test_fourier_expansion_cuspform_prec_loss(): #Test that prec_loss is working as intended
    S = AutomorphicFormSpace(Gamma0(1),weight=12)
    MF = get_cuspform_q_expansion_approx(S,100,label=0,prec_loss=30)
    assert abs(MF.get_cusp_expansion(Cusp(1,0))[45]-(-548895690)) < 1e-69 #this wouldn't work without setting prec_loss
    print("test_fourier_expansion_cuspform_prec_loss ok")

def test_fourier_expansion_modform():
    S = AutomorphicFormSpace(Gamma0(3),weight=4)
    MF = get_modform_q_expansion_approx(S,50,label=0)
    assert abs(MF.get_cusp_expansion(Cusp(1,0),trunc_order=10)[3]-240) < 1e-40
    print("test_fourier_expansion_modform ok")

def test_fourier_expansion_basis_cuspform():
    S = AutomorphicFormSpace(Gamma0(4),weight=10) #Search for multiplicity three new-forms
    b = get_cuspform_basis_approx(S,50,labels=[1,2])
    f0 = b[0].get_cusp_expansion(Cusp(1,0),trunc_order=10)
    f1 = b[1].get_cusp_expansion(Cusp(1,0),trunc_order=10)
    assert f0[0] == 0 and f0[1] == 0 and f0[3] == 0
    assert abs(f0[8]-256) < 1e-40
    assert f1[0] == 0 and f1[1] == 0 and f1[2] == 0
    assert abs(f1[9]-72) < 1e-40
    print("test_fourier_expansion_basis_cuspform ok")

def test_fourier_expansion_basis_modform():
    S = AutomorphicFormSpace(Gamma0(4),weight=10) #Search for multiplicity three new-forms
    b = get_modform_basis_approx(AutomorphicFormSpace(Gamma0(4),4),50,labels=[0,2])
    f0 = b[0].get_cusp_expansion(Cusp(1,0),trunc_order=10)
    f1 = b[1].get_cusp_expansion(Cusp(1,0),trunc_order=10)
    assert f0[0] == 1 and f0[1] == 0 and f0[2] == 0
    assert abs(f0[8]-2160) < 1e-40
    assert f1[0] == 0 and f1[1] == 0 and f1[2] == 1
    assert abs(f1[6]-28) < 1e-40
    print("test_fourier_expansion_basis_modform ok")