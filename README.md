# Code to compute Modular Forms and Eisenstein Series on Noncongruence Subgroups of PSL(2,&#8484;)
This repository hosts code that has been developed for my PhD-thesis. For further information, see [arXiv:2207.13365](https://arxiv.org/abs/2207.13365), [Noncongruence Database](https://github.com/David-Berghaus/noncong_database).

## Dependencies
The code has been developed using [SageMath](https://www.sagemath.org/) version 9.2 (other versions of Sage might or might not work) as well as [Psage](https://github.com/fredstro/psage). For the case that you cannot/do not want to install the later, we have also copied the required files into the folder `psage` to make this repository "self sufficient".

## Installation
To install the code run
```bash
sage setup.py develop `psage`
```
Once the compilation is finished, go to the folder `unit_tests` and run
```python
sage: load("unit_tests.sage")
sage: run_unit_tests()
```
to make sure that everything works correctly.

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
sage: j_G = B.get_hauptmodul_q_expansion(15)
sage: j_G
Hauptmodul of weight 0 with leading order expansion at infinity given by:
q_6^-1 + ((-80/3*v + 148/3)*u^2)*q_6 + ((-220*v + 64)*u^3)*q_6^2 + ((464*v + 992/3)*u^4)*q_6^3 + ((10928/3*v - 4576/3)*u^5)*q_6^4 + ((-875839/243*v - 486421/243)*u^6)*q_6^5 + ((5312*v + 48448/3)*u^7)*q_6^6 + ((-3552116/81*v - 40255396/243)*u^8)*q_6^7 + ((196309792/243*v - 143118584/243)*u^9)*q_6^8 + ((-786407200/243*v + 338569424/243)*u^10)*q_6^9 + O(q_6^10)
```
Fourier expansions are returned as instances of the `FourierExpansion` class. By default, the expansions up to 10-th order are printed. We can access the full expansion via
```python
sage: j_G.get_cusp_expansion(Cusp(1,0))
q_6^-1 + (-1937780208140/93936267*w^8 + 2647164/386569*w^2)*q_6 + (-5328895572385/31312089*w^9 - 110754064/386569*w^3)*q_6^2 + (11239125207212/31312089*w^10 + 1240786976/1159707*w^4)*q_6^3 + (264700776431924/93936267*w^11 + 1653813536/386569*w^5)*q_6^4 + (52706752/386569*w^6 + 3503356/386569)*q_6^5 + (15074131072/1159707*w^7 - 5163264/386569*w)*q_6^6 + (-13117827734704/93936267*w^8 + 42625392/386569*w^2)*q_6^7 + (-33447329283592/31312089*w^9 - 785239168/386569*w^3)*q_6^8 + (34579490081584/10437363*w^10 + 3145628800/386569*w^4)*q_6^9 + (1499731150802752/93936267*w^11 + 9987105280/386569*w^5)*q_6^10 + (-151696640/386569*w^6 + 18143414/386569)*q_6^11 + (22798246912/386569*w^7 - 
28366848/386569*w)*q_6^12 + (-67436944554616/93936267*w^8 + 171935448/386569*w^2)*q_6^13 + (-155884848364196/31312089*w^9 - 3547563584/386569*w^3)*q_6^14 + O(q_6^15)
```
Moreover, the `Belyi map` can be used to compute Fourier expansions of spaces of `cusp forms` and `modular forms` rigorously. For this we run (the first parameter corresponds to the `weight` while the second parameter denotes the `truncation order`)
```python
sage: B.get_cuspforms(4,15)
[CuspForm of weight 4 with leading order expansion at infinity given by:
 q_6 + ((4*v - 6)*u)*q_6^2 + ((28/3*v - 140/3)*u^2)*q_6^3 + ((-84*v + 224)*u^3)*q_6^4 + ((-672*v - 112/3)*u^4)*q_6^5 + ((3776/3*v - 400/3)*u^5)*q_6^6 + ((-102701/243*v - 27791/243)*u^6)*q_6^7 + ((-1054412/243*v + 7219348/243)*u^7)*q_6^8 + ((538636/27*v - 1477280/81)*u^8)*q_6^9 + O(q_6^10)])

sage: B.get_modforms(4,15)
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

### Hejhal's method
To compute Fourier expansions of modular forms numerically, using an optimized version of Hejhal's method, we can run
```python
sage: from psage.modform.maass.automorphic_forms import AutomorphicFormSpace
sage: from classes.fourier_expansion import get_hauptmodul_q_expansion_approx, get_cuspform_basis_approx, get_modform_basis_approx

sage: digit_prec = 50 #Compute coefficients to about 50 digits precision

sage: get_hauptmodul_q_expansion_approx(AutomorphicFormSpace(G,0),digit_prec)
Hauptmodul of weight 0 with leading order expansion at infinity given by:
1.000000000000000000000000000000000000000000000000*q_6^-1 + 0.0000000000000000000000000000000000000000000000000 + (-1.388983637482258528107138076635774326727017422309 + 5.016851501461998136341128266268848850335043646534*I)*q_6 + (6.998031044495970114687963449848038883098420993437 + 4.504136743636433523935354884887510289659450352209*I)*q_6^2 + (-3.512673610933919578738824203927641677521862598476 + 9.620911095204126577534732242676748912025407787146*I)*q_6^3 + (2.094356720009019678564632648311143927361526949870 - 16.23972128159629096243530616166816650495984318081*I)*q_6^4 + (8.870011644807860670298940067488886433373885562314 - 0.1526396550151222686065547830303095198469284518880*I)*q_6^5 + (5.896478618157735372967788703995928301000169284024 + 10.65387003159377190681326352327856535598346781938*I)*q_6^6 + (-7.102011599173434322441160409328590357785970189447 + 41.39193782819187367845381988814291456296782534338*I)*q_6^7 + (40.67104315040637664617150784384462137410245416033 + 37.61506266854387755161499630632108573866992332052*I)*q_6^8 + (-11.18048794365418621081119907300687097511058015156 + 74.25239079134411962653170377372871107862602489451*I)*q_6^9 + O(q_6^10)

sage: get_cuspform_basis_approx(AutomorphicFormSpace(G,4),digit_prec)
[CuspForm of weight 4 with leading order expansion at infinity given by:
 0.0000000000000000000000000000000000000000000000000 + 1.000000000000000000000000000000000000000000000000*q_6 + (0.4061258690242318718333555002648768682048322105118 + 1.800809415314694589228035628743749642947368430484*I)*q_6^2 + (-0.5722526997437003216282905545622017210545653261049 - 5.174031443714948823224351017198511799045825297871*I)*q_6^3 + (5.457657968765788827887711904575325501792031573910 - 6.282806363486977442870426757537387739492321058763*I)*q_6^4 + (-0.3090688252749344911976098578742973346302244306591 - 10.23744512366325604060961937994297363104070298019*I)*q_6^5 + (2.757128997180531220916472304483480150380742168491 - 5.540791011410894803434287670576599574584240982661*I)*q_6^6 + (0.8700116448078606702989400674888864333738493314920 - 0.1526396550151222686065547830303095198468580682568*I)*q_6^7 + (2.265667497222557019484589711446371876425634861251 + 17.33470003266323465504954059357681226192948590460*I)*q_6^8 + (-4.190026061852314639242479956469911119918202356384 + 0.3209425108888763297411129942038273585424564997926*I)*q_6^9 + O(q_6^10)]
 
 sage: get_modform_basis_approx(AutomorphicFormSpace(G,4),digit_prec)
[ModForm of weight 4 with leading order expansion at infinity given by:
 1.000000000000000000000000000000000000000000000000 + 0.0000000000000000000000000000000000000000000000000*q_6 + 0.0000000000000000000000000000000000000000000000000*q_6^2 + (5.173867968909958021375225882611907408025101824397e-48 - 7.970323278073582313339687657733545079142527762251e-49*I)*q_6^3 + (1.521600339715448373107679223977448577684990119567e-48 - 1.690006571972357747090536448860321771465182930549e-47*I)*q_6^4 + (-2.481823997698783086166639991159310772652454396656e-47 + 8.824315575114062602244189376941089893248108439657e-48*I)*q_6^5 + (239.9999999999999999999999999999999999999999999997 + 2.710574002751411907670920075732540763231046905499e-49*I)*q_6^6 + (-2.496990767074187374580873494609326250935329826270e-48 - 9.332762177652844612402267175999807446596701979134e-47*I)*q_6^7 + (-6.223461591778605168274369519418530624655598066330e-47 + 1.816222685565602054010243890294696748061962944949e-47*I)*q_6^8 + (5.721009063778472937971074082829849712664081280362e-47 + 1.110029227581787666878157992006010749856859285721e-45*I)*q_6^9 + O(q_6^10),
 ModForm of weight 4 with leading order expansion at infinity given by:
 0.0000000000000000000000000000000000000000000000000 + 1.000000000000000000000000000000000000000000000000*q_6 + 0.0000000000000000000000000000000000000000000000000*q_6^2 + (-3.922472674451917699470857262395952095563173685918 - 0.3143598845059013737664455018593258974215676047610*I)*q_6^3 + (15.56961126657719867821959419302920548111307097088 
- 2.223337024813179898668839840812346812291088352924*I)*q_6^4 + (-8.570691348242511613465697697226770028194863625042 - 31.94540342790802704797863241436137033115271776987*I)*q_6^5 + (-23.38257653178727285203545774821251364908427233504 + 22.92846755091636610947226438242966273592271078635*I)*q_6^6 + (67.45591183459758628268347857974629132905021352851 + 16.13619210159863982412150563463272066952499578644*I)*q_6^7 + (-29.10058892204735348602258852954910002606495774264 - 74.35665548155587899612687971664107082944201905954*I)*q_6^8 + (-101.1055988790414513948803130837117586125454121268 + 92.82761528176435592570644430612204498481043277230*I)*q_6^9 + O(q_6^10),
 ModForm of weight 4 with leading order expansion at infinity given by:
 0.0000000000000000000000000000000000000000000000000 + 0.0000000000000000000000000000000000000000000000000*q_6 + 1.000000000000000000000000000000000000000000000000*q_6^2 + (-2.168735504573075062036992828741925592016852219463 - 2.349498803263025433159017870009486907820436436223*I)*q_6^3 + (-3.350219974708217377842566707833750374508608359810 
+ 4.859671559209047449457905515339185901624257693116*I)*q_6^4 + (12.45568901326175894257567535442336438489045677671 - 1.778669619850543918935071872649877449832870682352*I)*q_6^5 + (-11.92883054653896394660571697221726303840162214564 - 17.20576977087227245013909048171313460360696771426*I)*q_6^6 + (-16.54277398308318732549883382690088090228445301094 + 33.24474606846536882060572602345959744750544589586*I)*q_6^7 + (52.19043085789084480106078249708879803483242526495 - 5.647667235559523938442526972121452234333748525241*I)*q_6^8 + (-37.33346907663913873807894495454271160265929787356 - 62.23738027826932901253870681918951924693253664783*I)*q_6^9 + O(q_6^10)]
```
 
 ### Eisenstein series
 To compute Eisenstein series numerically we run
```python
sage: from eisenstein.eisenstein_computation import compute_eisenstein_series
sage: eisenstein_series = compute_eisenstein_series(cuspforms,modforms)
sage: eisenstein_series
[ModForm of weight 4 with leading order expansion at infinity given by:
 1.000000000000000000000000000000000000000000000000 + 0.0000000000000000000000000000000000000000000000000*q_6 + 0.0000000000000000000000000000000000000000000000000*q_6^2 + (5.173867968909958021375225882611907408025101824397e-48 - 7.970323278073582313339687657733545079142527762251e-49*I)*q_6^3 + (1.521600339715448373107679223977448577684990119567e-48 - 1.690006571972357747090536448860321771465182930549e-47*I)*q_6^4 + (-2.481823997698783086166639991159310772652454396656e-47 + 8.824315575114062602244189376941089893248108439657e-48*I)*q_6^5 + (239.9999999999999999999999999999999999999999999997 + 2.710574002751411907670920075732540763231046905499e-49*I)*q_6^6 + (-2.496990767074187374580873494609326250935329826270e-48 - 9.332762177652844612402267175999807446596701979134e-47*I)*q_6^7 + (-6.223461591778605168274369519418530624655598066330e-47 + 1.816222685565602054010243890294696748061962944949e-47*I)*q_6^8 + (5.721009063778472937971074082829849712664081280362e-47 + 1.110029227581787666878157992006010749856859285721e-45*I)*q_6^9 + O(q_6^10), 
 ModForm of weight 4 with leading order expansion at infinity given by:
 0.0000000000000000000000000000000000000000000000000 + 1.000000000000000000000000000000000000000000000000*q_6 + (-3.959109788096710882221472451439743149088870430971 - 7.795413438144329803802188079360076192337431485404*I)*q_6^2 + (-13.65152525436440654764120628976826327472286261639 + 25.89375372074392402064795167715072324903597684807*I)*q_6^3 + (66.71664893834763385362532291318665737264181200048 + 4.653039549767070150131102051030926738554270863798*I)*q_6^4 + (-71.74959669473780057680001369295817005615242718986 - 122.0007006414948858372588629873477037325319571613*I)*q_6^5 + (-110.2811156399130224563304862496987597331715929018 + 184.0381650063516921508005160942826614490270642449*I)*q_6^6 + (392.1071104831046577268891085343816117261130389296 + 13.47435515105882973830847314756741893350040983067*I)*q_6^7 + (-279.7541356387557678479425397874621868205459842135 - 458.8429069014838697794072043238070499869638942437*I)*q_6^8 + (-438.4644066102331779945844497978370662214171744273 + 630.2620632595361783359588630807033599555620303355*I)*q_6^9 + O(q_6^10)]
```
These results are normalized to a reduced row echelon basis. To obtain a more canonical normalization for which the leading order coefficients are 1 for one cusp and zero at the others we can run
```python
sage: from eisenstein.eisenstein_computation import echelon_basis_to_eisenstein_basis
sage: echelon_basis_to_eisenstein_basis(eisenstein_series)
[ModForm of weight 4 with leading order expansion at infinity given by:
 1.000000000000000000000000000000000000000000000000 + (0.09600261606478410651446132739708791352637943378200 + 0.1600435981787040832542649659356889835823386710947*I)*q_6 + (0.8675211189862639012472634923855411991140381793971 - 1.382010259439953007974070885363760895325463350190*I)*q_6^2 + (-5.454711653514514669605066818403045202708050477885 + 0.3010288745927257385684503677435384241093328711083*I)*q_6^3 + (5.660283641144634922635684610893993488920887039381 + 11.12427652394911701792105552660416447074543118740*I)*q_6^4 + (12.63728212669975718840986391144647481930810247703 - 23.19545004621678474830962404877918535941417230443*I)*q_6^5 + (199.9585942662027830631160760894444598037093026086 + 0.01835873819872895534684551182952741287722945056718*I)*q_6^6 + (35.48682410246816884848043828010572232945921714869 + 64.04780617745837623362792493106383419497741385541*I)*q_6^7 + (46.57774094302660467443812829326969555762071980609 - 88.82297789831234121028781356031453273507262880430*I)*q_6^8 + (-142.9631384854657842949046435700282511408229961404 - 9.666614427888156543010183076143852706846395848912*I)*q_6^9 + O(q_6^10),
 ModForm of weight 4 with leading order expansion at infinity given by:
 0.0000000000000000000000000000000000000000000000000 + (-3.456094178332227834520607786295164886949659616155 - 5.761569534433346997153538773684803408964192159400*I)*q_6 + (-31.23076028350550044490148572587948316810537445821 + 49.75236933983830828706655187309539223171668060684*I)*q_6^2 + (196.3696195265225281057824054625096272974898172039 - 10.83703948533812658846421323876738326793598336013*I)*q_6^3 + (-203.7702110812068572148846459921837656011519334179 - 400.4739548621682126451579989577499209468355227462*I)*q_6^4 + (-454.9421565611912587827551008120730934950916891726 + 835.0362016638042509391464657560506729389102029594*I)*q_6^5 + (1441.490606416699809727821260779999447066465106076 - 0.6609145751542423924864384258629868635802602219901*I)*q_6^6 + (-1277.525667688854078545295778083806003860531817354 - 2305.721022388501544410605297518298031019186898795*I)*q_6^7 + (-1676.798673948957768279772618557709040074345913016 + 3197.627204339244283570361288171323178462614636954*I)*q_6^8 + (5146.672985476768234616567168521017041069627861052 + 347.9981194039736355483665907411786974464702505948*I)*q_6^9 + O(q_6^10)]
```
