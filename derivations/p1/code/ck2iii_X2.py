"""Class (iii) readout: pure background (first-order data f = 0), gauge of checkpoint 1.5
(l = 1 parts E_z = F_z = 0 of slot-1 h_ij), a2, e2 = l = 1 projections of slot 2.
Remaining free parameters (higher-l gauge, which do not shift a2, and smooth first-order data,
which class (iii) excludes) are set to zero; we report whether any of them still enters a2."""
import pickle, sympy as sp
import __main__
import ckpt2_axi as C
__main__.Field = C.Field
import trunc as T
from trunc import x, y, z, R
d = pickle.load(open('ckpt2_axi.pkl', 'rb'))
M, q, E = C.M, C.q, sp.Symbol('E')
nof = {sp.Symbol('Ao2s1_1'): 0}
cond3 = [sp.expand(sp.numer(sp.together(c)).subs(nof)) for c in d['log'][3]['cond']]
sol3 = {k: sp.sympify(v).subs(nof) for k, v in d['log'][3]['sol'].items()}
h, A, _ = C.ansatz(3, 1)
coef = [[T.radial_part(h[i + 1][j + 1], 2) for j in range(3)] for i in range(3)]
Dv, Ev, Fv = T.l1_tensor(coef)
gauge = [sp.expand(Ev[2]), sp.expand(Fv[2])]
w = lambda a: sp.Rational(3, a + 2) if a % 2 else 0
lam1 = lambda pre: sp.expand(sum(w(a) * sol3.get(sp.Symbol(f'{pre}o3s2_{a}'), sp.Symbol(f'{pre}o3s2_{a}'))
                                 for a in range(0, 9)))
a2v, e2v = sp.symbols('a2v e2v')
eqs = cond3 + gauge + [sp.expand(a2v + lam1('h') / 2), sp.expand(e2v - lam1('A'))]
syms = set().union(*[e.free_symbols for e in eqs])
params = sorted([s_ for s_ in syms if 'o3s' in str(s_)], key=str)
print(len(params), 'slot parameters;', len(eqs), 'equations')
tail = [sp.Symbol('ho3s2_1'), sp.Symbol('Ao3s2_1')]
head = [p_ for p_ in params if p_ not in tail]
cols = head + tail + [a2v, e2v]
A_, b_ = sp.linear_eq_to_matrix(eqs, cols)
gens = sorted(A_.free_symbols | b_.free_symbols, key=str)
dom = sp.QQ.frac_field(*gens)
from sympy.polys.matrices import DomainMatrix
Aug = A_.row_join(b_)
dm = DomainMatrix([[dom.from_sympy(Aug[i, j]) for j in range(Aug.cols)] for i in range(Aug.rows)], Aug.shape, dom)
rr, piv = dm.rref(); rr = rr.to_Matrix()
free = [j for j in range(len(cols)) if j not in piv]
print('free columns:', [cols[j] for j in free])
ia = cols.index(a2v)
if ia in piv:
    row = rr.row(list(piv).index(ia))
    # a2 = -(sum_{free j} row[j] col_j) + b
    dep = {cols[j]: sp.factor(-row[j]) for j in free if row[j] != 0}
    print('a2 depends on:', dep)
    a2expr = sum(v * k for k, v in dep.items()) + sp.factor(row[-1])
    X = sp.factor(sp.expand((6 * M * a2expr - 6 * q * e2v)))
    Xfree0 = sp.factor(X.subs({k: 0 for k in dep if k != e2v}))
    print('6 M a2 - 6 q e2 =', X)
    print('with the free (non-e2) parameters set to zero:', Xfree0)
    pickle.dump({'X': X, 'X0': Xfree0, 'dep': dep}, open('ck2iii_X2.pkl', 'wb'))
else:
    print('a2 still not a pivot')
