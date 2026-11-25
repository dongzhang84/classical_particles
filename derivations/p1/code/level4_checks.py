"""Level 4 checks: first-order residuals vanish on shell, and products of two
field-type first-order fields (E x E, E x B, B x B) give nothing in block I."""
import sys
import sympy as sp
import trunc as T
from trunc import x, y, z, R
from level4_l1 import build
import level4_Q as LQ
from series import DIRECTIONS


def onshell_zero(e):
    for nv in DIRECTIONS[:4]:
        v = sp.expand(sp.sympify(e).subs({x: R * nv[0], y: R * nv[1], z: R * nv[2]}, simultaneous=True))
        if sp.simplify(v) != 0:
            return False
    return True


T.set_limits(1, 1, 1)
h1, A1 = LQ.field1()
h2, A2 = LQ.field2('EBG')
g, A = build(T.mat(lambda i, j: T.add(h1[i][j], h2[i][j])), [T.add(A1[i], A2[i]) for i in range(4)])
Ein, Mx = T.einstein_maxwell(g, A)
for key in [(0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1)]:
    comps = [Ein[i][j].get(key, 0) for i in range(4) for j in range(4)] + [Mx[i].get(key, 0) for i in range(4)]
    print(key, 'first-order / flat residual vanishes on shell:', all(onshell_zero(c) for c in comps), flush=True)

# products of two field-type pieces: put the E/B field in both slots
for w1, w2 in [('E', 'E'), ('E', 'B'), ('B', 'B')]:
    ha, Aa = LQ.field2(w1)
    hb, Ab = LQ.field2(w2)
    swap = lambda dct: {(k[0], k[2], k[1]): v for k, v in dct.items()}   # mu2 -> mu1
    ha = T.mat(lambda i, j: swap(ha[i][j])); Aa = [swap(Aa[i]) for i in range(4)]
    g, A = build(T.mat(lambda i, j: T.add(ha[i][j], hb[i][j])), [T.add(Aa[i], Ab[i]) for i in range(4)])
    Ein, Mx = T.einstein_maxwell(g, A)
    for rp in (-2,):
        out = LQ.project(Ein, Mx, (1, 1, 1), rp)
        blockI = {k: [sp.factor(e) for e in out[k]] for k in ('A', 'D', 'E', 'F', 'X')}
        other = {k: [sp.factor(e) for e in out[k]] for k in ('B0i_static', 'C0i_static', 'mag_0i')}
        print(w1, 'x', w2, 'block I:', blockI)
        print('        block II (static 0i):', other, flush=True)
