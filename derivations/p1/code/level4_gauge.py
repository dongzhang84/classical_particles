"""Level 4 gauge check: the derived first-order equation of motion

  M zeta''_i = -M R_0i0j zeta^j + q zeta^k d_k F_i0 + q F_ij zeta'^j
               + q f_i0 - (1/2) q h_00 E_i - M d_t h_i0 + (1/2) M d_i h_00

must be invariant under a far-zone gauge transformation
  h -> h + L_xi g,  f -> f + L_xi F,  zeta^i -> zeta^i - xi^i(t, 0)
(background in accelerated Fermi coordinates, M a = q E along the world line)."""
import itertools
import sympy as sp

t = sp.Symbol('t')
xs = sp.symbols('x y z')
X = (t,) + xs
M, q = sp.symbols('M q')
E = [sp.Function(f'E{i}')(t) for i in range(3)]
a = [q * e / M for e in E]
# background Riemann (only R_0i0j, R_0ijk, R_ijkl at x = 0 needed), field gradients
R0i0j = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'S{min(i,j)}{max(i,j)}'))
Fij = sp.Matrix(3, 3, lambda i, j: 0 if i == j else (sp.Symbol(f'B{3-i-j}') * (1 if (i, j) in [(0, 1), (1, 2), (2, 0)] else -1)))
dE = sp.Matrix(3, 3, lambda i, k: sp.Symbol(f'dE{i}{k}'))     # d_k F_i0
dF = [[[sp.Symbol(f'dF{min(i,j)}{max(i,j)}{k}') * (1 if i < j else -1) if i != j else 0
        for k in range(3)] for j in range(3)] for i in range(3)]   # d_k F_ij
W = [[[sp.Symbol(f'W{i}{j}{k}') for k in range(3)] for j in range(3)] for i in range(3)]  # R_0ijk part
Rs = [[[[sp.Symbol(f'Q{i}{j}{k}{l}') for l in range(3)] for k in range(3)] for j in range(3)] for i in range(3)]

ax = sum(a[i] * xs[i] for i in range(3))
g = sp.zeros(4, 4)
g[0, 0] = -(1 + ax)**2 - sum(R0i0j[i, j] * xs[i] * xs[j] for i in range(3) for j in range(3))
for i in range(3):
    g[0, i+1] = g[i+1, 0] = -sp.Rational(2, 3) * sum(W[j][i][k] * xs[j] * xs[k] for j in range(3) for k in range(3))
    for j in range(3):
        g[i+1, j+1] = (1 if i == j else 0) - sp.Rational(1, 3) * sum(Rs[i][k][j][l] * xs[k] * xs[l] for k in range(3) for l in range(3))
F = sp.zeros(4, 4)
for i in range(3):
    F[i+1, 0] = E[i] + sum(dE[i, k] * xs[k] for k in range(3))
    F[0, i+1] = -F[i+1, 0]
    for j in range(3):
        F[i+1, j+1] = Fij[i, j] + sum(dF[i][j][k] * xs[k] for k in range(3))

# gauge vector: general Taylor expansion to second order in x, arbitrary in t
xi = []
for mu in range(4):
    e = sp.Function(f'Z{mu}')(t)
    for k in range(3):
        e += sp.Function(f'Z{mu}_{k}')(t) * xs[k]
        for l in range(k, 3):
            e += sp.Function(f'Z{mu}_{k}{l}')(t) * xs[k] * xs[l]
    xi.append(e)

D = lambda f, m: sp.diff(f, X[m])
Lg = sp.Matrix(4, 4, lambda m, nu: sum(xi[r] * D(g[m, nu], r) + g[r, nu] * D(xi[r], m) + g[m, r] * D(xi[r], nu) for r in range(4)))
LF = sp.Matrix(4, 4, lambda m, nu: sum(xi[r] * D(F[m, nu], r) + F[r, nu] * D(xi[r], m) + F[m, r] * D(xi[r], nu) for r in range(4)))
at0 = lambda e: sp.expand(e).subs({xx: 0 for xx in xs})

dzeta = [-xi[i+1].subs({xx: 0 for xx in xs}) for i in range(3)]
for i in range(3):
    lhs = M * sp.diff(dzeta[i], t, 2)
    rhs = (-M * sum(R0i0j[i, j] * dzeta[j] for j in range(3))
           + q * sum(dE[i, k] * dzeta[k] for k in range(3))
           + q * sum(Fij[i, j] * sp.diff(dzeta[j], t) for j in range(3))
           + q * at0(LF[i+1, 0])
           - sp.Rational(1, 2) * q * at0(Lg[0, 0]) * E[i]
           - M * at0(D(Lg[i+1, 0], 0))
           + sp.Rational(1, 2) * M * at0(D(Lg[0, 0], i+1)))
    print(i, sp.simplify(sp.expand(lhs - rhs)))
