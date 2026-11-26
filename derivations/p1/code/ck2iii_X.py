"""Class (iii), static axisymmetric electric: read X in the gauge of checkpoint 1.5
(l = 1 parts E, F of the slot-1 (first-order r^2) h_ij vanish).
Solve the order-3 conditions + the two gauge conditions for a2 (= -ho3s2_1/2 per unit z)."""
import pickle, sympy as sp
import __main__
import ckpt2_axi as C
__main__.Field = C.Field
import trunc as T
from trunc import x, y, z, R, nn
d = pickle.load(open('ckpt2_axi.pkl', 'rb'))
cond3 = [sp.numer(sp.together(c)) for c in d['log'][3]['cond']]
M, q, E = C.M, C.q, sp.Symbol('E')
h, A, syms1 = C.ansatz(3, 1)
# l = 1 projection of the slot-1 h_ij: coefficient of R^2 after x -> R n
coef = [[T.radial_part(h[i + 1][j + 1], 2) for j in range(3)] for i in range(3)]
Dv, Ev, Fv = T.l1_tensor(coef)
gauge = [sp.expand(Ev[2]), sp.expand(Fv[2])]
print('gauge conditions (l=1 E_z, F_z of slot-1 h_ij):', gauge)
sol3 = d['log'][3]['sol']
w = lambda a: sp.Rational(3, a + 2) if a % 2 else 0          # P1 weight of cos^a
lam1 = lambda pre: sum(w(a) * sp.sympify(sol3.get(sp.Symbol(f'{pre}o3s2_{a}'), sp.Symbol(f'{pre}o3s2_{a}')))
                       for a in range(0, 9))
H1, A1 = sp.expand(lam1('h')), sp.expand(lam1('A'))          # l = 1 parts: h00 -> H1 z, A0 -> A1 z
a2v, e2v = sp.symbols('a2v e2v')
extra = [sp.expand(a2v + H1 / 2), sp.expand(e2v - A1)]
allsyms = set().union(*[c.free_symbols for c in cond3 + gauge + extra])
s1 = sorted([s_ for s_ in allsyms if 'o3s1' in str(s_) or 'o3s2' in str(s_)], key=str)
solr, condr, freer = C.solve_ordered(cond3 + gauge + extra, s1, [e2v, a2v])
print('conditions on (e2, a2):')
for c in condr:
    print('   ', sp.factor(c))
print('free slot parameters left:', len(freer))
pickle.dump({'cond': condr}, open('ck2iii_X.pkl', 'wb'))
