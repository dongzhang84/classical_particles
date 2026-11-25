"""Level 4: the 1/r part of the first-order near-zone field is fixed only up to
a gauge transformation xi ~ eps * r^0. Check that such pieces do not change the
combination W.Q = (2/3) Q^D + Q^E that gives the motion (nor the identity
-Q^A + Q^D + Q^F)."""
import sympy as sp
import trunc as T
from trunc import t, x, y, z, R
from level4_l1 import build, M, q, xs, n, dot, cross, half, dl
import level4_Q as LQ

V = sp.symbols('V_x V_y V_z')
S = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f's{min(i,j)}{max(i,j)}'))
gauges = {
    'xi_i = M (V.n) n_i': [0] + [M * dot(V, n) * n[i] for i in range(3)],
    'xi_i = M V_i': [0] + [M * V[i] for i in range(3)],
    'xi_0 = M V.n': [-M * dot(V, n), 0, 0, 0],
    'xi_i = M n_i (l=0)': [0] + [M * n[i] for i in range(3)],
    'xi_i = M S_ij n_j (l=2)': [0] + [M * sum(S[i, j] * n[j] for j in range(3)) for i in range(3)],
    'xi_0 = M (l=0)': [-M, 0, 0, 0],
}
T.set_limits(1, 1, 1)
h2, A2 = LQ.field2('EB')
for name, xi in gauges.items():
    # h = L_xi eta (flat part only: eps * r^-1), as a mu1 field at order eps
    h = [[sp.expand(sum(T.eta[k, j] * T.d_expr(xi[k], i) + T.eta[i, k] * T.d_expr(xi[k], j) for k in range(4)))
          for j in range(4)] for i in range(4)]
    h1 = T.mat(lambda i, j: T.S({(1, 1, 0): h[i][j]}))
    A1 = [{} for _ in range(4)]
    g, A = build(T.mat(lambda i, j: T.add(h1[i][j], h2[i][j])), [T.add(A1[i], A2[i]) for i in range(4)])
    Ein, Mx = T.einstein_maxwell(g, A)
    out = LQ.project(Ein, Mx, (1, 1, 1), -2)
    WQ = [sp.factor(sp.Rational(2, 3) * out['D'][k] + out['E'][k]) for k in range(3)]
    ID = [sp.factor(-out['A'][k] + out['D'][k] + out['F'][k]) for k in range(3)]
    print(f'{name:26s} W.Q = {WQ}   identity = {ID}', flush=True)
