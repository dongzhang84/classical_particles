"""Checkpoint 1.5, near-zone second order, l = 2 sector on RN (to eps^2):
the static tidal response of a charged black hole.  Data: the vacuum tidal
field at rbar^2 (eps^0) in the RW-type form of checkpoint 1; unknowns: the most
general static l = 2 metric and potential terms at rbar^1 (eps^1) and rbar^0
(eps^2).  Reports the solution and the free (gauge) parameters."""
import sympy as sp
import trunc as T
from trunc import t, x, y, z, R
from series import DIRECTIONS
from ckpt15_l1 import rn_background, solve_ordered, put, empty
M, q = sp.symbols('M q')
xs = (x, y, z); n = [xi / R for xi in xs]
Eij = sp.Matrix([[1, 2, 0], [2, -3, 1], [0, 1, 2]]); Bij = sp.Matrix([[2, 0, 1], [0, -1, 3], [1, 3, -1]])
lev = sp.LeviCivita; kd = lambda i, j: 1 if i == j else 0
Enn = sum(Eij[i, j] * n[i] * n[j] for i in range(3) for j in range(3))
En = [sum(Eij[i, j] * n[j] for j in range(3)) for i in range(3)]
Bn = [sum(Bij[k, l] * n[l] for l in range(3)) for k in range(3)]
nBn = [sum(lev(i, j, k) * n[j] * Bn[k] for j in range(3) for k in range(3)) for i in range(3)]
Bnn = sum(Bij[i, j] * n[i] * n[j] for i in range(3) for j in range(3))
# magnetic h_ij structure n_(i eps_j)kl n^k B^l_m n^m
mag_ij = lambda i, j: (n[i] * nBn[j] + n[j] * nBn[i]) / 2
def form(c, p):
    """General static l = 2 field: c = dict of 10 coefficients, times R^p."""
    h = [[0] * 4 for _ in range(4)]
    h[0][0] = c['00'] * Enn * R**p
    for i in range(3):
        h[0][i+1] = h[i+1][0] = (c['0iB'] * nBn[i] + c['0iE'] * En[i] + c['0iEn'] * n[i] * Enn) * R**p
        for j in range(3):
            h[i+1][j+1] = (c['d'] * kd(i, j) * Enn + c['nn'] * n[i] * n[j] * Enn + c['E'] * Eij[i, j]
                           + c['nE'] * (n[i] * En[j] + n[j] * En[i]) / 2 + c['mag'] * mag_ij(i, j)) * R**p
    A = [c['A0'] * Enn * R**p] + [(c['AB'] * nBn[i] + c['AE'] * En[i] + c['AEn'] * n[i] * Enn) * R**p for i in range(3)]
    return h, A
keys = ['00', '0iB', '0iE', '0iEn', 'd', 'nn', 'E', 'nE', 'mag', 'A0', 'AB', 'AE', 'AEn']
T.set_limits(2, 1, 0)
P = empty()
data = {k: 0 for k in keys}; data.update({'00': -1, '0iB': sp.Rational(2, 3), 'd': -1})
h, A = form(data, 2); put(P, 0, h, A)
C1 = {k: sp.Symbol(f'{k}_1') for k in keys}; C2 = {k: sp.Symbol(f'{k}_2') for k in keys}
for C, p, k in [(C1, 1, 1), (C2, 0, 2)]:
    h, A = form(C, p); put(P, k, h, A)
g00, gss, A0 = rn_background(2)
G = T.mat(lambda a, b: T.add(g00 if a == b == 0 else (gss if a == b else {}), P[0][a][b]))
AA = [T.add(A0 if a == 0 else {}, P[1][a]) for a in range(4)]
Ein, Mx = T.einstein_maxwell(G, AA)
eqs = []
for cpt in [Ein[a][b] for a in range(4) for b in range(a, 4)] + Mx:
    for N in range(3):
        ex = cpt.get((N, 1, 0), 0)
        if ex == 0: continue
        for nv in DIRECTIONS[:6]:
            e = sp.expand(sp.expand(sp.sympify(ex).subs({x: R*nv[0], y: R*nv[1], z: R*nv[2]}, simultaneous=True)) * R**6)
            eqs += [cf for cf in sp.Poly(e, R).coeffs() if cf != 0]
unk = list(C2.values()) + list(C1.values())
sol, cond, free = solve_ordered(eqs, unk, [])
print('conditions:', cond)
print('free:', free)
for k in keys:
    for C, tag in [(C1, 'rbar^1'), (C2, 'rbar^0')]:
        v = sol.get(C[k], C[k])
        if v != 0: print(tag, k, sp.factor(v))
print('check 00_1 in sol:', sol.get(C1['00'], 'absent'), ' 0iB_1:', sol.get(C1['0iB'], 'absent'), ' A0_1:', sol.get(C1['A0'], 'absent'))
print('n eqs', len(eqs), 'example', eqs[:2])
