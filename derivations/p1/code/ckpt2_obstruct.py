"""Find the obstruction to the order-3 equations of ckpt2_axi: a left null
vector w of the coefficient matrix with w.b != 0, reported by (component, N,
monomial).  Uses the saved order-2 solution (run ckpt2_axi.py --radial --log
--keepall first)."""
import pickle, random, sympy as sp
import __main__
import ckpt2_axi as C
__main__.Field = C.Field
import trunc as T
from trunc import x, y, z, R
from sympy.polys.matrices import DomainMatrix
dd = pickle.load(open('ckpt2_axi_o2.pkl', 'rb'))
F = dd['F']; dcols = dd['log'][2]['free']
m = 3; unk = []
for k in range(0, C.EPS + 1):
    h, A, syms = C.ansatz(m, k); F.add((k, m, 0), h, A); unk += syms
    if k >= 2:
        hL, AL, sL = C.ansatz(m, k, withlog=True); F.add((k, m, 0), hL, AL); unk += sL
T.set_limits(C.EPS, m, 0)
g00, gss, A0 = C.rn_background(C.EPS)
G = T.mat(lambda a_, b: T.add(g00 if a_ == b == 0 else (gss if a_ == b else {}), F.h[a_][b]))
AA = [T.add(A0 if a_ == 0 else {}, F.A[a_]) for a_ in range(4)]
Ein, Mx = T.einstein_maxwell(G, AA)
Ls = sp.Symbol('LOGR')
rows, labels = [], []
names = [f'E{a_}{b}' for a_ in range(4) for b in range(a_, 4)] + [f'Mx{a_}' for a_ in range(4)]
for nm, cpt in zip(names, [Ein[a_][b] for a_ in range(4) for b in range(a_, 4)] + Mx):
    for N in range(0, 4):
        ex = cpt.get((N, m, 0), 0)
        if ex == 0: continue
        e = sp.expand(sp.sympify(ex).subs(y, 0).subs(sp.log(R), Ls))
        lo = max([-t_.as_powers_dict().get(R, 0) for t_ in sp.Add.make_args(e)] + [0])
        num = sp.expand(sp.rem(sp.Poly(sp.expand(e * R**lo), x), sp.Poly(x**2 + z**2 - R**2, x)).as_expr())
        for mon, cf in sp.Poly(num, x, z, R, Ls).terms():
            rows.append(cf); labels.append((nm, N, mon, lo))
cols = unk + list(dcols)
Am, bm = sp.linear_eq_to_matrix(rows, cols)
gens = sorted(Am.free_symbols | bm.free_symbols, key=str)
num = {g: sp.Rational(random.randint(3, 17), random.randint(3, 17)) for g in gens}
An = Am.subs(num); bn = bm.subs(num)
# left null space of An (rows): solve An^T w = 0
ns = An.T.nullspace()
print('left null space dim', len(ns))
for w in ns:
    val = (w.T * bn)[0]
    if val != 0:
        print('obstruction w.b =', val)
        for i, wi in enumerate(w):
            if wi != 0:
                nm, N, mon, lo = labels[i]
                print(f'   {wi}  x {nm} N={N} monomial x^{mon[0]} z^{mon[1]} R^{mon[2]-lo} log^{mon[3]} ; b_i = {bn[i]}')
        break
