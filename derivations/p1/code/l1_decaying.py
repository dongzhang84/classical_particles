"""In the gauge a = 0, K = 0 the remaining polar fields H0, H2, h1 are the
gauge invariants. Reduce the equations to one ODE and look at its second
(decaying) solution near the horizon."""
import pickle
import sympy as sp
from l1_response import r, th, M, q, H0, H2, K, h1, a, b, c
eqs = pickle.load(open('l1_eqs.pkl', 'rb'))
pol = [eqs['Ep'][i, j] for i in range(4) for j in range(i, 4) if eqs['Ep'][i, j] != 0] + [m for m in eqs['Mp'] if m != 0]
red = []
for e in pol:
    e = e.subs({a: 0, K: 0, b: 0, c: 0}).doit()
    e = sp.simplify(e.subs(th, sp.pi / 3))
    if e != 0:
        red.append(sp.factor(sp.numer(sp.together(e))))
for e in red:
    print(sp.collect(sp.expand(e), [H0.diff(r, 2), H0.diff(r), H0, H2.diff(r), H2, h1.diff(r), h1]))
    print('---')

D = r**2 - 2 * M * r + q**2
u, up, upp, uppp, v = sp.symbols('u up upp uppp v')   # H0, H0', H0'', H0''', H2
h1e = (r**2 * (r - M) * v - r * (r**2 - 3 * M * r + 2 * q**2) * u + r**2 * D * up) / (2 * D)
# check the algebraic equation reproduces red[2]
H2p = up + 4 * h1e / r**2                          # from the last equation
e2 = red[1]
sub = {H0.diff(r, 2): upp, H0.diff(r): up, H0: u, H2.diff(r): H2p, H2: v, h1.diff(r): sp.Symbol('h1p'), h1: h1e}
# h1' : differentiate h1e totally, using v' = H2p
h1p = sp.diff(h1e, r) + sp.diff(h1e, u) * up + sp.diff(h1e, up) * upp + sp.diff(h1e, v) * H2p
sub[sp.Symbol('h1p')] = h1p
e2s = sp.factor(e2.subs({H0.diff(r, 2): upp, H0.diff(r): up, H2.diff(r): H2p, h1.diff(r): h1p}).subs({H0: u, H2: v, h1: h1e}))
print('eq2 reduced:', e2s)
vsol = sp.solve(e2s, v)
print('H2 from eq2:', vsol)
