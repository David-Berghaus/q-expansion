# Code to compute Modular Forms and Eisenstein Series on Noncongruence Subgroups of PSL(2,&#8484;)
This repository hosts code that has been developed for my PhD-thesis. For further information, see [arXiv:2207.13365](https://arxiv.org/abs/2207.13365), [Noncongruence Database](https://github.com/David-Berghaus/noncong_database).

## Dependencies
The code has been developed using [SageMath](https://www.sagemath.org/) version 9.2 (other versions of Sage might or might not work) as well as [Psage](https://github.com/fredstro/psage).

## Installation
To install the code run
```bash
sage setup.py develop
```

## Example usage
### Belyi Maps
Let `G` be a genus zero subgroup
```python
sage: from psage.modform.arithgroup.mysubgroup import MySubgroup
sage: G = MySubgroup(o2='(1 2)(3 4)(5)(6 7)',o3='(1)(2 3 5)(4 6 7)')
```
To compute the associated `Belyi map` we run
```python
sage: from classes.belyi_map import BelyiMap
sage: B = BelyiMap(G)
sage: B
(x + 14*u)(x^2 + ((2*v - 10)*u)*x + (184/3*v - 164/3)*u^2)^3 / (x + (6*v - 16)*u)
```
The `Belyi map` object contains a variety of functions and parameters. For example, to determine the expressions of `v` and `u` from the above example we run
```python
sage: B._Kv
Number Field in v with defining polynomial x^2 - x + 1 with v = 0.50000000000000000? - 0.866025403784439?*I
sage: B.get_u_str()
'(125248356/96889010407*v - 199546416/96889010407)^(1/6)'
```
Because `G` is a genus zero subgroup, we can compute the Fourier expansion of its `Hauptmodul`
```python
sage: B.get_hauptmodul_q_expansion(10)
Hauptmodul of weight 0 with leading order expansion at infinity given by:
q_6^-1 + ((-80/3*v + 148/3)*u^2)*q_6 + ((-220*v + 64)*u^3)*q_6^2 + ((464*v + 992/3)*u^4)*q_6^3 + ((10928/3*v - 4576/3)*u^5)*q_6^4 + ((-875839/243*v - 486421/243)*u^6)*q_6^5 + ((5312*v + 48448/3)*u^7)*q_6^6 + ((-3552116/81*v - 40255396/243)*u^8)*q_6^7 + ((196309792/243*v - 143118584/243)*u^9)*q_6^8 + ((-786407200/243*v + 338569424/243)*u^10)*q_6^9 + O(q_6^10)
```
Moreover, the `Belyi map` can be used to compute Fourier expansions of spaces of `cusp forms` and `modular forms` rigorously. For this we run (the first parameter corresponds to the `weight` while the second parameter denotes the `truncation order`)
```python
sage: B.get_cuspforms(4,10)
[CuspForm of weight 4 with leading order expansion at infinity given by:
 q_6 + ((4*v - 6)*u)*q_6^2 + ((28/3*v - 140/3)*u^2)*q_6^3 + ((-84*v + 224)*u^3)*q_6^4 + ((-672*v - 112/3)*u^4)*q_6^5 + ((3776/3*v - 400/3)*u^5)*q_6^6 + ((-102701/243*v - 27791/243)*u^6)*q_6^7 + ((-1054412/243*v + 7219348/243)*u^7)*q_6^8 + ((538636/27*v - 1477280/81)*u^8)*q_6^9 + O(q_6^10)])

sage: B.get_modforms(4,10)
[ModForm of weight 4 with leading order expansion at infinity given by:
 1 + 240*q_6^6 + O(q_6^10),
 ModForm of weight 4 with leading order expansion at infinity given by:
 q_6 + ((-104/3*v + 16/3)*u^2)*q_6^3 + ((-380*v + 360)*u^3)*q_6^4 + ((-2432*v + 1424/3)*u^4)*q_6^5 + ((-15904/3*v - 5152/3)*u^5)*q_6^6 + ((-4550573/243*v - 6201935/243)*u^6)*q_6^7 + ((-47168/3*v - 355072/3)*u^7)*q_6^8 + ((2876464/9*v - 58380968/81)*u^8)*q_6^9 + O(q_6^10),
 ModForm of weight 4 with leading order expansion at infinity given by:
 q_6^2 + ((-2*v + 10)*u)*q_6^3 + ((-44*v + 52)*u^2)*q_6^4 + ((-304*v + 288)*u^3)*q_6^5 + ((-1632*v + 824)*u^4)*q_6^6 + ((-7552*v + 800)*u^5)*q_6^7 + ((-5732782/243*v - 2174842/243)*u^6)*q_6^8 + ((-325322/9*v - 2510672/27)*u^7)*q_6^9 + O(q_6^10)]
```
If we are interested in computing Fourier expansions using floating-point numbers with rigorous error bounds, we specify the additional parameter `digit_prec`
```python
sage: B.get_cuspforms(4,10,digit_prec=100)
[CuspForm of weight 4 with leading order expansion at infinity given by:
 1.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000*q_6 + ([0.40612586902423187183335550026487686820483221051224283744087953569362805591694546003920748787847540 +/- 4.83e-99] + [1.8008094153146945892280356287437496429473684304866141043864220502882601089687405949482072680314734 +/- 3.36e-98]*I)*q_6^2 + ([-0.57225269974370032162829055456220172105456532610557227329887099527821812923086466881592822071445 +/- 3.71e-96] + [-5.17403144371494882322435101719851179904582529787733639632235041741687734522376882410268813211592 +/- 5.20e-96]*I)*q_6^3 + ([5.4576579687657888278877119045753255017920315739170005707736919028138204813423995872513929184125 +/- 8.37e-95] + [-6.2828063634869774428704267575373877394923210587696171244007797432784081535418607599470865524615 +/- 3.11e-95]*I)*q_6^4 + ([-0.309068825274934491197609857874297334630224430653008054218849916216846626768493851407892260752 +/- 7.57e-94] + [-10.237445123663256040609619379942973631040702980202373806268081493946088592298259509418266439239 +/- 7.96e-94]*I)*q_6^5 + ([2.7571289971805312209164723044834801503807421684934001520786426794791312225744810264040245499 +/- 3.00e-92] + [-5.5407910114108948034342876705765995745842409826484591193637050689354290971474137146681017555 +/- 3.61e-92]*I)*q_6^6 + ([0.870011644807860670298940067488886433373849331486030480497071798315327797091347992758119491 +/- 1.23e-91] + [-0.152639655015122268606554783030309519846858068249675104350080381158850831875412785356928862 +/- 2.85e-91]*I)*q_6^7 + ([2.26566749722255701948458971144637187642563486125791422815896227471733802372053099487069485 +/- 1.57e-90] + [17.33470003266323465504954059357681226192948590462371810751914208942639532198253465082752925 +/- 3.45e-90]*I)*q_6^8 + ([-4.1900260618523146392424799564699111199182023563424398414215427932180154369708621521746875 +/- 3.67e-89] + [0.3209425108888763297411129942038273585424564998133946104228282854323256367001945505625534 +/- 2.70e-89]*I)*q_6^9 + O(q_6^10)]
```
