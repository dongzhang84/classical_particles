"""Adjudicate level 1 against GW (106): pure-gauge perturbations must give zero
linearized Ricci, so their (D, H, K) content must be a null vector of the
correct matrix.  For each gauge vector xi_a (time-even, l=1), compute
h = L_xi g with g = eta + (2M/r) delta to O(M), read off its D/H/K pieces,
and apply both matrices.
"""
import sympy as sp
from series import t, x, y, z, r, X, DIRECTIONS, radial_coeff
M = sp.Symbol('M', positive=True)
c = sp.symbols('c_x c_y c_z')
xs = (x, y, z)
cx = sum(ci * xi for ci, xi in zip(c, xs))
psi = 2 * M / r
g = sp.diag(-1 + psi, 1 + psi, 1 + psi, 1 + psi)
ginv = sp.diag(-(1 + psi), 1 - psi, 1 - psi, 1 - psi)  # O(M)

def lie(xi_low):
    xi_up = [sp.expand(sum(ginv[a, b] * xi_low[b] for b in range(4))) for a in range(4)]
    h = sp.zeros(4, 4)
    for a in range(4):
        for b in range(4):
            e = sum(xi_up[cc] * sp.diff(g[a, b], X[cc]) for cc in range(4))
            e += sum(g[cc, b] * sp.diff(xi_up[cc], X[a]) + g[a, cc] * sp.diff(xi_up[cc], X[b]) for cc in range(4))
            h[a, b] = sp.series(sp.expand(e), M, 0, 2).removeO()
    return h

def V(n):  # symbols
    return sp.symbols(f'{n}_x {n}_y {n}_z')
names = ['DA', 'DC', 'DE', 'DF', 'HA', 'HD', 'HE', 'HF', 'KB', 'KC']
U = {k: V(k) for k in names}
dot = lambda a, b: sum(p * q for p, q in zip(a, b))
d = lambda i, j: 1 if i == j else 0

def decompose(h):
    """Fit h's l=1 even pieces at the orders used in GW (98)-(101)."""
    eqs = []
    for n in DIRECTIONS:
        eqs.append(radial_coeff(h[0, 0], n, 0, 1) - dot(U['DA'], n))
        eqs.append(radial_coeff(h[0, 0], n, 0, 0) - dot(U['HA'], n))
        for i in range(3):
            eqs.append(radial_coeff(h[0, i+1], n, 1, 0) - U['DC'][i])
            eqs.append(radial_coeff(h[0, i+1], n, 1, -1) - (dot(U['KB'], n) * n[i] + U['KC'][i]))
            for j in range(i, 3):
                eqs.append(radial_coeff(h[i+1, j+1], n, 0, 1)
                           - (sp.Rational(1, 2) * (n[i] * U['DE'][j] + n[j] * U['DE'][i]) + d(i, j) * dot(U['DF'], n)))
                eqs.append(radial_coeff(h[i+1, j+1], n, 0, 0)
                           - (dot(U['HD'], n) * n[i] * n[j] + sp.Rational(1, 2) * (U['HE'][i] * n[j] + U['HE'][j] * n[i])
                              + dot(U['HF'], n) * d(i, j)))
    unk = [s for k in names for s in U[k]]
    sol = sp.solve(eqs, unk, dict=True)
    assert sol, 'gauge perturbation not of the ansatz form at these orders'
    return [sp.factor(sol[0].get(U[k][0], 0)) for k in names]

h_ = sp.Rational(1, 2)
ours = sp.Matrix([
 [3*M/2, -M, M/2, M/2, 1, 0, 0, 0, 1, -1],
 [0, 3*M, 0, 0, 0, 0, 0, 0, 3*h_, 3*h_],
 [0, -M, 0, 0, 0, 0, 0, 0, -h_, -h_],
 [3*M, 0, 0, 3*M, 3*h_, -3*h_, 0, -3*h_, 3, 0],
 [0, -2*M, 0, -2*M, -1, 1, 0, 1, -1, 1],
 [-3*M/2, M, M/2, -5*M/2, -h_, 3*h_, 0, 3*h_, -1, 0]])
gw = sp.Matrix([
 [-3*M, 0, -M, -M, -2, 0, 0, 0, -2, 2],
 [0, -6*M, 0, 0, 0, 0, 0, 0, -3, -3],
 [0, 2*M, 0, 0, 0, 0, 0, 0, 1, 1],
 [-6*M, 0, 0, -6*M, -3, 3, 0, 3, -6, 0],
 [0, 4*M, 0, 4*M, 2, -2, 0, -2, 2, -2],
 [0, 3*M, -3*M, 9*M, 1, -3, 0, -3, 2, 0]]) * sp.Rational(-1, 2)

tests = {
 'xi_0 = t (c.x)':        [t * cx, 0, 0, 0],
 'xi_i = (c.x) x_i':      [0] + [cx * xi for xi in xs],
 'xi_i = r^2 c_i':        [0] + [r**2 * ci for ci in c],
}
for name, xi in tests.items():
    v = sp.Matrix(decompose(lie(xi)))   # x-components, with c_y = c_z = 0 below
    v = v.subs({c[1]: 0, c[2]: 0})
    print(f'{name}:  (D,H,K) =', list(v.T))
    print('   ours * v =', list(sp.simplify(ours * v).T))
    print('   GW   * v =', list(sp.simplify(gw * v).T))
