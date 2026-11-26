"""X_(iii) (static electric, Fermi radial background gauge, 1.5 gauge hE = hF = 0) evaluated on
MP data with an extremal body, for several L."""
import pickle, sys, sympy as sp
M, q = sp.symbols('M q', positive=True)
cond = pickle.load(open('ck2iii_struct_radial.pkl', 'rb'))['cond']
A2, F2 = sp.symbols('A2 F2')
rel = [c for c in cond if c.has(A2)][0]
a2sol = sp.solve(sp.numer(sp.together(rel)), A2)[0]
hE, hF = sp.symbols('hEk1 hFk1')
X = sp.factor(sp.expand(6 * (M * a2sol - q * F2)).subs({hE: 0, hF: 0}))
print('X = 6(M a2 - q e2) =', X)
for Lv in sys.argv[1:]:
    mp = pickle.load(open(f'ck2iii_mp_L{Lv}.pkl', 'rb'))
    sub = {sp.Symbol('E'): mp['E'], sp.Symbol('Ao2s0_2'): mp['A0'][1], sp.Symbol('S4o2s0_0'): mp['hij'][0][sp.Symbol('s4R')]}
    val = sp.simplify(X.subs(sub).subs(q, M))
    print(f'L = {Lv}:  X(MP, q = M) =', val)
