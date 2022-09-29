# Sage script to reproduce the database data

# Define some helper functions that are needed to transform u-v-factored expansions to L
def convert_from_K_to_L(expression_in_K, v_L):
    """
    Given an expression in K, convert the expression efficiently to L by plugging in v(w).
    """
    if expression_in_K == 0:
        return 0
    coeffs = list(expression_in_K.polynomial())
    res = coeffs[-1]
    for i in range(len(coeffs)-2,-1,-1): #Horner's method
        res = res*v_L+coeffs[i]
    return res

def transform_u_v_factored_q_expansion_to_L(q_expansion, L, v_L, u_interior_K, principal_cusp_width):
    """
    Given a q_expansion that has coefficients of the form (expression_in_K)*u**u_pow, convert the coefficients to L.
    """
    if principal_cusp_width == 1:
        u = L(u_interior_K)
    else:
        u = L.gen()
    leading_order_exponent = q_expansion.valuation()
    coeffs = list(q_expansion)
    coeffs_L = []
    for coeff in coeffs:
        if coeff == 0 or coeff == 1:
            u_pow = 0
        else:
            u_pow = coeff.degree()
        coeffs_L.append(convert_from_K_to_L(coeff[u_pow],v_L)*u**u_pow)
    P = LaurentSeriesRing(L,q_expansion.variable())
    return P(coeffs_L).shift(leading_order_exponent).O(q_expansion.prec())

res = {} #The dictionary in which we store the results
res["G"] = ArithmeticSubgroup_Permutation(S2="(2,4)(3,5)(6,7)",S3="(1,2,3)(4,5,6)")
principal_cusp_width = 3
res["monodromy_group"] = "S7"
P.<T> = PolynomialRing(QQ)
K.<v> = NumberField(P([-1, 1]),embedding=1.00000000000000+0.000000000000000*1j)
res["K"] = K
res["v"] = K.gen()
u_interior_K = 2/823543
L.<w> = NumberField(P([-2, 0, 0, 823543]),embedding=0.0134415052246602+0.000000000000000*1j) #Base ring of the q-expansions in raw format
v_L = 1

# The pretty form of q-expansions involves powers of u, where u is given by
u_str = "(2/823543)^(1/3)"
res["u_str"] = u_str
# With embedding
res["u"] = QQbar(L.gen())
Pu.<u> = PolynomialRing(K)
Pq.<q_3> = LaurentSeriesRing(Pu)

# Store the q-expansions in pretty form into res
res["q_expansions"] = {}
weights = [0, 2, 4]
for weight in weights:
    res["q_expansions"][weight] = {}
res["q_expansions"][0]["hauptmodul_pretty"] = Pq([0, 148932*u^2, 35666932*u^3, 7392301056*u^4, 1360018740528*u^5, 231535133811191*u^6, 36898595322937344*u^7, 5585391081472928292*u^8, 810336872638231259496*u^9, 113443808393800865820672*u^10, 15418082849911364113493184*u^11, 16327731265829819340180438819/8*u^12, 264171668904610362297445220352*u^13, 67010545446339775585882819153941/2*u^14, 4174104253199792518383073583171333*u^15, 511591978382014069189403568768769536*u^16, 61781937959578124679247522427279856408*u^17]).O(17)
res["q_expansions"][0]["hauptmodul_pretty"] += Pq.gen()**(-1) #Don't forget the first coeff, which we cannot construct from a list
res["q_expansions"][2]["modforms_pretty"] = []
res["q_expansions"][2]["modforms_pretty"].append(Pq([1, -462*u, -84420*u^2, 807828*u^3, -891458736*u^4, -82305718992*u^5, 5155138704870*u^6, -807981764899218*u^7, -57396539567144736*u^8, -829520378016134700*u^9, -368800915551641445600*u^10, -38181751584385599534960*u^11, 3403490699635925493047523/2*u^12, -626937472524463105637199687/4*u^13, -12963545628243868733327853615*u^14, 53103234035865914749482329871/2*u^15, -135694801849991608603081593600864*u^16, -8835738290946071781371659355069250*u^17]).O(18))
res["q_expansions"][4]["cuspforms_pretty"] = []
res["q_expansions"][4]["cuspforms_pretty"].append(Pq([0, 1, 18*u, -8640*u^2, -1823860*u^3, -124044624*u^4, 17247354624*u^5, 155396886869*u^6, 243506186396892*u^7, -4171668253982784*u^8, 1198246183061776128*u^9, -162373510774205223120*u^10, -9417599666207795520000*u^11, -1178361267593053377469187/8*u^12, -54207930629208500080095411/2*u^13, 2267486458269250660022120760*u^14, -511267509110035007571720758662*u^15, 71053413143588231341797908343702*u^16, 3726442093739814209803125723117120*u^17]).O(19))
res["q_expansions"][4]["modforms_pretty"] = []
res["q_expansions"][4]["modforms_pretty"].append(Pq([1, 0, 0, 98825160*u^3, 0, 0, 366240459338460*u^6, 0, 0, 469178525829958565880*u^9, 0, 0, 503685257516490257269597095*u^12, 0, 0, 357983664191053218690168156516135*u^15]).O(18))
res["q_expansions"][4]["modforms_pretty"].append(Pq([0, 1, 0, -648*u^2, -2679868*u^3, -243999216*u^4, 50397461280*u^5, -1673081058283*u^6, 52997921199936*u^7, 34089621896320200*u^8, -890087190481230336*u^9, -383910426241434871920*u^10, 34928013196996800473664*u^11, -16287935659753267777718531/8*u^12, -203468035110912882423694080*u^13, 32150711324824180991718931377*u^14, -2054677545970584161767797137698*u^15, -59687780250618171823181598653310*u^16]).O(18))
res["q_expansions"][4]["modforms_pretty"].append(Pq([0, 0, 1, -444*u, 47556*u^2, 6664144*u^3, -1841672592*u^4, 101582108064*u^5, 10583792510942*u^6, -2125627230572388*u^7, 116018520752389248*u^8, 12307606414846091600*u^9, -2463645159066921999648*u^10, 104927599945557044446176*u^11, 39192015510290807196365861/4*u^12, -3320358318506103370188534513/2*u^13, 85745002047808286344226465502*u^14, 7263399633011466842498861499834*u^15]).O(18))

# Now consider the q-expansions in raw format defined over L
Pq.<q_3> = LaurentSeriesRing(L)
for weight in weights:
    if "hauptmodul_pretty" in res["q_expansions"][weight]:
        res["q_expansions"][weight]["hauptmodul_raw"] = transform_u_v_factored_q_expansion_to_L(res["q_expansions"][weight]["hauptmodul_pretty"],L,v_L,u_interior_K,principal_cusp_width)
    if "cuspforms_pretty" in res["q_expansions"][weight]:
        res["q_expansions"][weight]["cuspforms_raw"] = [transform_u_v_factored_q_expansion_to_L(cuspform,L,v_L,u_interior_K,principal_cusp_width) for cuspform in res["q_expansions"][weight]["cuspforms_pretty"]]
    if "modforms_pretty" in res["q_expansions"][weight]:
        res["q_expansions"][weight]["modforms_raw"] = [transform_u_v_factored_q_expansion_to_L(modform,L,v_L,u_interior_K,principal_cusp_width) for modform in res["q_expansions"][weight]["modforms_pretty"]]

# Add the Eisenstein scaling constants as well
res["q_expansions"][2]["eisenstein_basis_factors"] = [[1]]
res["q_expansions"][2]["eisenstein_canonical_normalizations"] = {0: [1.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000]}
res["q_expansions"][4]["eisenstein_basis_factors"] = {0: [1, 0, 1.75980923022701600090333071689577426066151771111808524114463933521704778949884796506392761835170974e-99 - 4.39952307556754000225832679223943565165379427779521310286159833804261947374711991265981904587927434e-100*I], 1: [0, 1, 8.22682055381244703090631807007847837576719352691577684112465115062155030586764576785329335559316087 - 5.74318959104143157570434967624471116052511198105345286842408320981886879940691350659008669793930784e-100*I]}
res["q_expansions"][4]["eisenstein_canonical_normalizations"] = {0: [1.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000, 3.16706173692649208844927033697277120152412813660255716494080821612177184792287077002281233466379612 + 2.82891885401576307851733123960097444625576227488696921718094793956844855977959070253175798229442381e-100*I], 1: [0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000, -1.78147222702115179975271456454718380085732207683893840527920462156849666445661480813783193824838531 - 1.28133937281628828797685665977079487890947990554737278169000963757306373496924737485688350184186350e-100*I]}

# Now add the curve
F.<x> = FunctionField(L)
res["curve"] = (F([-6729651447039507456/678223072849*w, 6772700364362496/96889010407, 5164155191808/117649*w^2, -51543431136/117649*w, -66250656/117649, 833868*w^2, 1848*w, 1]))/(F([37949472/343*w, 2299968/2401, 1280664*w^2, 1848*w, 1]))

# Now add the embeddings
res["embeddings"] = {('(1)(2 4)(3 5)(6 7)', '(1 2 5)(3 6 4)(7)', '(1 2 3)(4 5 6 7)'): 1}
