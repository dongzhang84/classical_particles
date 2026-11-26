"""Compare the near-zone l = 1 relations of ckpt2_axi.py with the test-body
condition  M d_z sqrt(-g00) = q d_z A_0  at x = 0 (static), built from the
smooth parts of the far-zone fields:
    g00 = -1 - 2 a z + lam (c1 + l1 z) + lam^2 (c2 + l2 z),
    A0  =  E z + lam e1 z + lam^2 e2 z ,
where (near-zone order m, slot k):  c_k = h00 constant at (k, k),
l_k = h00 z-coefficient at (k+1, k),  e_k = A0 z-coefficient at (k+1, k)."""
import pickle, sympy as sp
import __main__
from ckpt2_axi import Field
__main__.Field = Field
M, q = sp.symbols('M q', positive=True)
d = pickle.load(open('ckpt2_axi.pkl', 'rb'))
log = d['log']
S = lambda n: sp.Symbol(n)
sol = {}
for m in sorted(log):
    sol.update({k: v for k, v in log[m]['sol'].items()})
def val(name):
    s = S(name)
    v = sol.get(s, s)
    for _ in range(4):
        v = sp.sympify(v).subs(sol)
    return sp.factor(v)
a, E, lam = S('a'), S('E'), S('lam')
c1, l1, e1 = val('ho1s1_0'), val('ho2s1_1'), val('Ao2s1_1')
print('a =', val('a'))
print('c1 =', c1, '  l1 =', l1, '  e1 =', e1)
maxo = max(log)
if maxo >= 3:
    c2, l2, e2 = val('ho2s2_0'), val('ho3s2_1'), val('Ao3s2_1')
    print('c2 =', c2, '\nl2 =', l2, '\ne2 =', e2)
else:
    c2, l2, e2 = S('c2'), S('l2'), S('e2')
C1, L1, E1, C2, L2, E2 = sp.symbols('C1 L1 E1 C2 L2 E2')
g00 = lambda zz: -1 - 2 * val('a') * zz + lam * (C1 + L1 * zz) + lam**2 * (C2 + L2 * zz)
zz = sp.Symbol('zz')
expr = M * sp.diff(sp.sqrt(-g00(zz)), zz) - q * sp.diff(E * zz + lam * E1 * zz + lam**2 * E2 * zz, zz)
expr = sp.series(expr.subs(zz, 0), lam, 0, 3).removeO()
for k in (1, 2):
    tb = sp.solve(sp.expand(expr).coeff(lam, k), [L1, L2][k - 1])[0]
    tb = tb.subs({C1: c1, E1: e1, L1: l1})
    if k == 2:
        tb = tb.subs({C2: c2, E2: e2})
    got = [l1, l2][k - 1]
    print(f'order {k}: near-zone l{k} - test-body l{k} =', sp.factor(sp.expand(got - tb)))
