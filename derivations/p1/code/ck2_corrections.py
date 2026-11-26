"""Checkpoint 2: difference between the two definitions of the second-order regular field,
   ours  (R = total - transported singular field: j^R contains j^H - L_xi h^H)  and
   Gralla's (j^R_G = j - (j^P - L_xi h^P)), for a charged body in electrovac.
The difference delta j = -L_xi h^H, delta A = -L_xi A^H(1) enters the second-order perturbed
Lorentz equation linearly:
   dZ_i = 1/2 d_i dj_00 - d_0 dj_0i - 1/2 a0_i dj_00 + (q/M) dF_i0     (at x = 0)
(dj_ik terms cancel using a0 = q E0/M).  Fields (P gauge, first order, near gamma):
   h^H_00 = -2 a1(t).x - dE_kl x^k x^l,  h^H_0i = O(x^2),  h^H_ij = O(x^2),
   A^H_0 = f(t).x + (1/2) G_kl x^k x^l,  A^H_i = (1/2)(b x x)_i + O(x^2),  M a1 = q f.
xi^mu = Z^mu(t) + w^mu_nu(t) x^nu + O(x^2)  (arbitrary smooth gauge data)."""
import sympy as sp
t = sp.Symbol('t'); X = sp.symbols('x y z'); co = (t,) + X
q, M = sp.symbols('q M', positive=True)
f = [sp.Function(f'f{i}')(t) for i in range(3)]
a1 = [q * fi / M for fi in f]
a0 = sp.symbols('a0_1:4'); E0 = [M * ai / q for ai in a0]
dE = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'dE{min(i,j)}{max(i,j)}'))
G = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'G{min(i,j)}{max(i,j)}'))
b = sp.symbols('b1:4')
Z = [sp.Function(f'Z{m}')(t) for m in range(4)]
w = [[sp.Function(f'w{m}{n}')(t) for n in range(4)] for m in range(4)]
xi = [Z[m] + sum(w[m][n] * co[n] for n in range(4)) for m in range(4)]
h = sp.zeros(4, 4)
h[0, 0] = -2 * sum(a1[i] * X[i] for i in range(3)) - sum(dE[i, j] * X[i] * X[j] for i in range(3) for j in range(3))
bxx = [b[1] * X[2] - b[2] * X[1], b[2] * X[0] - b[0] * X[2], b[0] * X[1] - b[1] * X[0]]
A = [sum(f[i] * X[i] for i in range(3)) + sp.Rational(1, 2) * sum(G[i, j] * X[i] * X[j] for i in range(3) for j in range(3))] \
    + [sp.Rational(1, 2) * bxx[i] for i in range(3)]
d = lambda e, m: sp.diff(e, co[m])
Lh = sp.Matrix(4, 4, lambda m, n: sum(xi[r] * d(h[m, n], r) + h[r, n] * d(xi[r], m) + h[m, r] * d(xi[r], n) for r in range(4)))
LA = [sum(xi[r] * d(A[n], r) + A[r] * d(xi[r], n) for r in range(4)) for n in range(4)]
dj, dA = -Lh, [-v for v in LA]
F = lambda Av, m, n: d(Av[n], m) - d(Av[m], n)
at0 = {xx: 0 for xx in X}
dZ = []
for i in range(3):
    e = (sp.Rational(1, 2) * d(dj[0, 0], i + 1) - d(dj[0, i + 1], 0) - sp.Rational(1, 2) * a0[i] * dj[0, 0]
         + q / M * F(dA, i + 1, 0))
    dZ.append(sp.simplify(sp.expand(e.subs(at0))))
# expected closed form, in terms of the displacement Zh = xi|_gamma = Z + w_{.0} t
Zh = [Z[m] + w[m][0] * t for m in range(4)]
Zd = [sp.diff(Zh[m], t) for m in range(4)]
Fik = lambda i, k: sp.diff(A[k + 1], X[i]) - sp.diff(A[i + 1], X[k])
expect = []
for i in range(3):
    ex = (sum(dE[i, k] * Zh[k + 1] for k in range(3)) + a1[i] * Zd[0]
          - a0[i] * sum(a1[k] * Zh[k + 1] for k in range(3))
          - q / M * (sum(Zh[k + 1] * G[i, k] for k in range(3)) + sum(Fik(i, k).subs(at0) * Zd[k + 1] for k in range(3))))
    expect.append(ex)
Zs = sp.symbols('Zh0:4'); Zds = sp.symbols('Zhd0:4')
rep = {}
for m in range(4):
    rep[sp.diff(Z[m], t)] = Zds[m] - sp.diff(w[m][0] * t, t)
    rep[Z[m]] = Zs[m] - w[m][0] * t
dZr = [sp.simplify(sp.expand(e.subs(rep, simultaneous=True))) for e in dZ]
print('after writing Z = Zh - w_{.0} t: depends on w?', any(e.has(*[w[m][n] for m in range(4) for n in range(4)]) for e in dZr))
print('matches closed form:', all(sp.simplify(dZ[i] - expect[i]) == 0 for i in range(3)))
for i in range(3):
    print(f'dZ_{i+1} =', dZ[i])
