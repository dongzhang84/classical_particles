"""Class (iii) (pure background, static, electric, axisymmetric E along z) with a STRUCTURED
l = 1 ansatz at third order, as in ck2nz.py.
Lower orders: the near-zone order-1 and order-2 fields from ckpt2_axi.py --light (pickled after
order 2), with the first-order datum f = 0 (class (iii)); free parameters there are zero except
the background slot-0 data (dE, curvature), which stay symbolic.
Order 3: general l = 1 ansatz (12 vectors along z; slots 1, 2, 3; powers 2, 1, 0), with a2, e2 at
slot 2; equations projected on l = 1; solved with slot-1 coefficients and (e2, a2) as data."""
import pickle, sys, time, sympy as sp
import __main__
import ckpt2_axi as C
__main__.Field = C.Field
import trunc as T
from trunc import t, x, y, z, R
from ckpt15_l1 import l1_form, NAMES, half, solve_ordered
from ck2nz import projected_equations
dd = pickle.load(open('ckpt2_axi_o2.pkl' if C.EPS == 3 else f'ckpt2_axi_o2_eps{C.EPS}.pkl', 'rb'))
F = dd['F'].subs({sp.Symbol('Ao2s1_1'): 0})
print('order-2 symbols kept:', sorted(set().union(*[sp.sympify(v).free_symbols for row in F.h for c in row for v in c.values()]) - {C.M, C.q, R, x, y, z}, key=str))
zv = lambda name: [0, 0, sp.Symbol(name)]
V = {}
a2, e2 = sp.symbols('A2 F2')
for k, p in [(kk, 3 - kk) for kk in range(0, C.EPS + 1)]:
    U = {nm: zv(f'{nm}k{k}') for nm in NAMES}
    if k == 2:
        U['h00'] = [0, 0, -2 * a2]; U['a0'] = [0, 0, e2]
    h, A = l1_form(U, p)
    F.add((k, 3, 0), h, A)
    V[k] = U
T.set_limits(C.EPS, 3, 0)
g00, gss, A0 = C.rn_background(C.EPS)
G = T.mat(lambda a_, b: T.add(g00 if a_ == b == 0 else (gss if a_ == b else {}), F.h[a_][b]))
AA = [T.add(A0 if a_ == 0 else {}, F.A[a_]) for a_ in range(4)]
t0 = time.time()
Ein, Mx = T.einstein_maxwell(G, AA)
print(f'Einstein-Maxwell to lambda^3: {time.time()-t0:.0f}s', flush=True)
eqs, tags = [], []
for N in range(C.EPS + 1):
    e_, t_ = projected_equations(Ein, Mx, (N, 3, 0), range(-8, 3))
    eqs += e_; tags += t_
print(len(eqs), 'projected equations', flush=True)
k1 = [V[1][nm][2] for nm in NAMES]
data = [s for s in k1 if isinstance(s, sp.Symbol)] + [e2, a2]
resp = [V[k][nm][2] for k in list(range(C.EPS, 1, -1)) + [0] for nm in NAMES if isinstance(V[k][nm][2], sp.Symbol) and V[k][nm][2] not in (a2, e2)]
sol, cond, free = solve_ordered(eqs, resp, data)
print('conditions (slot-1 coefficients, e2, a2):')
for c in cond:
    print('   ', c)
pickle.dump({'cond': cond}, open('ck2iii_struct.pkl', 'wb'))
