"""Checkpoint 1 (q = 0, vacuum): the static l = 2 near-zone second-order
solution of Gralla 2012, eqs. (11)-(13), with undetermined coefficients."""
import sympy as sp
import trunc as T
from trunc import t, x, y, z, R

M = sp.Symbol('M')
xs = (x, y, z); n = [xi / R for xi in xs]
Eij = sp.Matrix([[1, 2, 0], [2, -3, 1], [0, 1, 2]])          # numeric STF (traceless symmetric)
Bij = sp.Matrix([[2, 0, 1], [0, -1, 3], [1, 3, -1]])
c = sp.symbols('c1:4'); d = sp.symbols('d1:4'); e = sp.symbols('e1:4'); f_ = sp.symbols('f2:4')
Enn = sum(Eij[i, j] * n[i] * n[j] for i in range(3) for j in range(3))
Bn = [sum(Bij[k, l] * n[l] for l in range(3)) for k in range(3)]
lev = sp.LeviCivita
h = [[0] * 4 for _ in range(4)]
h[0][0] = Enn * (c[0] * R**2 + c[1] * M * R + c[2] * M**2)
for i in range(3):
    v = sum(lev(i, j, k) * n[j] * Bn[k] for j in range(3) for k in range(3)) * (d[0] * R**2 + d[1] * M * R + d[2] * M**2)
    h[0][i+1] = h[i+1][0] = v
    for j in range(3):
        h[i+1][j+1] = Enn * ((1 if i == j else 0) * (e[0] * R**2 + e[1] * M * R + e[2] * M**2)
                             + n[i] * n[j] * (f_[0] * M * R + f_[1] * M**2))
T.set_limits(2, 1, 0)                       # keep M^0..M^2 in the background
psi = 1 + M / (2 * R)
g00 = sp.series(-((1 - M / (2 * R)) / psi)**2, M, 0, 3).removeO()
gss = sp.series(psi**4, M, 0, 3).removeO()
def split(expr):
    expr = sp.expand(expr)
    return {(k, 0, 0): expr.coeff(M, k) * M**k for k in range(3) if expr.coeff(M, k) != 0}
g = T.mat(lambda a, b: T.add(split(g00) if a == b == 0 else (split(gss) if a == b else {}),
                             {(0, 1, 0): h[a][b]} if h[a][b] != 0 else {}))
# trunc keys: first slot counts powers of M; move the M^k factor into the key
def rekey(dct):
    out = {}
    for k, v in dct.items():
        v = sp.expand(v)
        for p in range(3):
            cp = v.coeff(M, p) if p else v.subs(M, 0)
            if cp != 0:
                kk = (p, k[1], k[2])
                out[kk] = out.get(kk, 0) + cp
    return out
g = T.mat(lambda a, b: rekey(g[a][b]))
Ein, Mx = T.einstein_maxwell(g, [{}, {}, {}, {}])
eqs = []
from series import DIRECTIONS
for a in range(4):
    for b in range(a, 4):
        for p in range(3):
            ex = Ein[a][b].get((p, 1, 0), 0)
            if ex == 0:
                continue
            for nv in DIRECTIONS[:5]:
                v = sp.expand(sp.sympify(ex).subs({x: R * nv[0], y: R * nv[1], z: R * nv[2]}, simultaneous=True))
                eqs += [v.coeff(R, kk) for kk in range(-6, 3) if v.coeff(R, kk) != 0]
sol = sp.solve(eqs + [c[0] + 1, d[0] - sp.Rational(2, 3)], list(c) + list(d) + list(e) + list(f_), dict=True)
print(sol)
