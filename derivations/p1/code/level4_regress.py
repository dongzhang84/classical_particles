"""Regression: the linear operator columns of level 2 recomputed with trunc.py
(sign conventions of the rows: R_ab - 8 pi T_ab and d_b(sqrt(-g) F^ab))."""
import sympy as sp
import trunc as T
from trunc import t, x, y, z, R
from level4_l1 import build, M, q, xs, n, dot, half, dl
import level4_Q as LQ

V = sp.symbols('V_x V_y V_z')
T.set_limits(1, 1, 0)
def run(name, h, a):
    hh = T.mat(lambda i, j: T.S({(0, 1, 0): h[i][j]}))
    AA = [T.S({(0, 1, 0): a[i]}) for i in range(4)]
    g, A = build(hh, AA)
    Ein, Mx = T.einstein_maxwell(g, A)
    for key, rp, tp in [((1, 1, 0), -2, 0), ((0, 1, 0), -2, 0)]:
        get = lambda e: T.radial_part(e.get(key, 0), rp)
        S00 = get(Ein[0][0]); Tij = [[get(Ein[i+1][j+1]) for j in range(3)] for i in range(3)]
        Dv, Ev, Fv = T.l1_tensor(Tij)
        X0 = T.l1_scalar(get(Mx[0]))
        res = [sp.factor(T.l1_scalar(S00)[0] / V[0]), sp.factor(Dv[0] / V[0]), sp.factor(Ev[0] / V[0]),
               sp.factor(Fv[0] / V[0]), sp.factor(X0[0] / V[0])]
        print(name, key, 'A D E F X =', res)

z4 = [[0] * 4 for _ in range(4)]
h = [row[:] for row in z4]; h[0][0] = sp.sqrt(1) * R * dot(V, n)
run('DA', h, [0, 0, 0, 0])
h = [row[:] for row in z4]
for i in range(3):
    for j in range(3):
        h[i+1][j+1] = R * dl(i, j) * dot(V, n)
run('DF', h, [0, 0, 0, 0])
run('E0', z4, [dot(V, xs), 0, 0, 0])
run('HA', [[dot(V, n) if (i, j) == (0, 0) else 0 for j in range(4)] for i in range(4)], [0, 0, 0, 0])
run('HP', z4, [dot(V, n), 0, 0, 0])
