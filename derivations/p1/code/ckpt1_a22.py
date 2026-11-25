"""Checkpoint 1 (q = 0, vacuum): the r^0 term a22 of the second-order
far-zone field in Gralla's P gauge (Gralla 2012, eqs. (46)-(48)).

In the near zone, a20 rbar^2 + a21 rbar + a22 is the static l = 2 part of the
second-order near-zone perturbation gbar^(2) on the Schwarzschild background
(isotropic coordinates).  It must solve the linearized vacuum equations order
by order in M.  We

  1. fix the M^0 and M^1 coefficients to the P-gauge values (eqs. (31)-(33),
     reproduced in checkpoint 1) and solve for the most general static l = 2
     M^2 coefficients;
  2. test whether the printed eqs. (46)-(48) lie in that solution space;
  3. compute a22 from the RW solution and the gauge vector v of eq. (29)
     (with the factor M): gbar_P^(2) = gbar_RW^(2) - L_zeta gbar^(0),
     zeta_i = M rbar^2 (2 E_ij n^j - n_i E_jk n^j n^k), and check it.
"""
import sympy as sp
import trunc as T
from trunc import t, x, y, z, R
from series import DIRECTIONS

M = sp.Symbol('M')
xs = (x, y, z); n = [xi / R for xi in xs]
Eij = sp.Matrix([[1, 2, 0], [2, -3, 1], [0, 1, 2]])          # numeric STF
Bij = sp.Matrix([[2, 0, 1], [0, -1, 3], [1, 3, -1]])
lev = sp.LeviCivita
kd = lambda i, j: 1 if i == j else 0
Enn = sum(Eij[i, j] * n[i] * n[j] for i in range(3) for j in range(3))
En = [sum(Eij[i, j] * n[j] for j in range(3)) for i in range(3)]
Bn = [sum(Bij[k, l] * n[l] for l in range(3)) for k in range(3)]
nBn = [sum(lev(i, j, k) * n[j] * Bn[k] for j in range(3) for k in range(3)) for i in range(3)]

# M^2 unknowns: 00: c Enn ; i0: d (eps n B n) ; ij: e delta Enn + f n n Enn + g E_ij + k n_(i E_j)n
c, d, e, f, g, k = sp.symbols('c d e f g k')


def tensor(p00, pi0, pd, pnn, pE, pnE):
    """Static l = 2 ansatz from radial profiles (functions of R, M)."""
    h = [[0] * 4 for _ in range(4)]
    h[0][0] = Enn * p00
    for i in range(3):
        h[0][i+1] = h[i+1][0] = nBn[i] * pi0
        for j in range(3):
            h[i+1][j+1] = (kd(i, j) * Enn * pd + n[i] * n[j] * Enn * pnn + Eij[i, j] * pE
                           + (n[i] * En[j] + n[j] * En[i]) / 2 * pnE)
    return h


T.set_limits(2, 1, 0)
psi = 1 + M / (2 * R)
g00 = sp.series(-((1 - M / (2 * R)) / psi)**2, M, 0, 3).removeO()
gss = sp.series(psi**4, M, 0, 3).removeO()


def rekey(expr, slot):
    expr = sp.expand(expr)
    out = {}
    for p in range(3):
        cp = expr.coeff(M, p) if p else expr.subs(M, 0)
        if cp != 0:
            out[(p, slot, 0)] = cp
    return out


def residuals(h):
    """Linearized vacuum equations for h on Schwarzschild, M^0..M^2, as a list."""
    G = T.mat(lambda a, b: T.add(rekey(g00 if a == b == 0 else (gss if a == b else 0), 0),
                                 rekey(h[a][b], 1)))
    Ein, _ = T.einstein_maxwell(G, [{}, {}, {}, {}])
    eqs = []
    for a in range(4):
        for b in range(a, 4):
            for p in range(3):
                ex = Ein[a][b].get((p, 1, 0), 0)
                if ex == 0:
                    continue
                for nv in DIRECTIONS[:5]:
                    v = sp.expand(sp.sympify(ex).subs({x: R * nv[0], y: R * nv[1], z: R * nv[2]},
                                                      simultaneous=True))
                    eqs += [v.coeff(R, kk) for kk in range(-6, 3) if v.coeff(R, kk) != 0]
    return eqs


r2 = sp.Rational
# ---- 1. P gauge through M^1 (eqs. (31)-(33)), general M^2
hP = tensor(-R**2 + 2 * M * R + c * M**2,
            r2(2, 3) * R**2 + r2(2, 3) * M * R + d * M**2,
            -R**2 - 2 * M * R + e * M**2,
            f * M**2,
            -4 * M * R + g * M**2,
            k * M**2)
eqs = residuals(hP)
sol = sp.solve(eqs, [c, d, e, f, g, k], dict=True)
print('1. P gauge, general M^2 solution:', sol)

# ---- 2. printed eqs. (46)-(48)
printed = {c: -3, d: r2(10, 3), e: -6, f: -1, g: 8, k: 0}
res = [sp.simplify(q.subs(printed)) for q in eqs]
print('2. printed (46)-(48) satisfy the equations:', all(q == 0 for q in res))
bad = sorted({q for q in res if q != 0}, key=str)[:6]
print('   sample nonzero residuals:', bad)
# which single components are off? test each printed value against the solution
if sol:
    for var, val in printed.items():
        s = sol[0].get(var, var)
        print(f'   {var}: printed {val}, required {s}')

# ---- 3. RW solution + gauge transformation by zeta
hRW = tensor(-R**2 + 2 * M * R - r2(3, 2) * M**2,
             r2(2, 3) * R**2 + r2(2, 3) * M * R - r2(1, 6) * M**2,
             -R**2 - 4 * M * R - 5 * M**2,
             2 * M * R + 4 * M**2, 0, 0)
assert all(sp.simplify(q) == 0 for q in residuals(hRW)), 'RW solution'
hRW_printed = tensor(-R**2 + 2 * M * R + r2(3, 2) * M**2,                 # eqs. (11)-(13) as printed
                     r2(2, 3) * R**2 + r2(2, 3) * M * R - r2(1, 6) * M**2,
                     -R**2 - 4 * M * R - 5 * M**2, 2 * M * R - 4 * M**2, 0, 0)
print('   printed (11)-(13) satisfy the equations:',
      all(sp.simplify(q) == 0 for q in residuals(hRW_printed)))
# zeta_i = M R^2 (2 E n - n Enn): -L_zeta gbar0, gbar0 = Schwarzschild to O(M)
zeta = [0] + [M * R**2 * (2 * En[i] - n[i] * Enn) for i in range(3)]       # lower index, flat
gb = [[(g00 if a == b == 0 else (gss if a == b else 0)) for b in range(4)] for a in range(4)]
gb = [[sp.expand(sp.series(v, M, 0, 2).removeO()) if v != 0 else 0 for v in row] for row in gb]
# zeta^a = eta^{ab} zeta_b  (corrections are O(M^2) in L_zeta g, i.e. beyond a22's M^2? no:
# zeta^a is O(M); g^{ab} - eta^{ab} = O(M) would make it O(M^2) and L_zeta g ~ O(M^2) at r^0 -- keep it)
ginv0 = [[(1 / g00 if a == b == 0 else (1 / gss if a == b else 0)) for b in range(4)] for a in range(4)]
zup = [sp.expand(sp.series(sum(ginv0[a][b] * zeta[b] for b in range(4)), M, 0, 3).removeO())
       for a in range(4)]
D = lambda f_, a: T.d_expr(f_, a)
Lz = [[sp.expand(sum(zup[cc] * D(gb[a][b], cc) + gb[cc][b] * D(zup[cc], a) + gb[a][cc] * D(zup[cc], b)
                     for cc in range(4))) for b in range(4)] for a in range(4)]
hP3 = [[sp.expand(sp.series(hRW[a][b] - Lz[a][b], M, 0, 3).removeO()) for b in range(4)]
       for a in range(4)]
assert all(sp.simplify(q) == 0 for q in residuals(hP3)), 'gauge-transformed RW solution'
print('3. RW - L_zeta gbar0 solves the equations; its M^0, M^1, M^2 coefficients:')
nv = DIRECTIONS[0]
probe = {x: R * nv[0], y: R * nv[1], z: R * nv[2]}
# read off profiles by matching against the ansatz on several directions
cc_ = sp.symbols('p0:6')
for pw in (0, 1, 2):
    amb = tensor(*[cc_[i] * M**pw * R**(2 - pw) for i in range(6)])
    eqf = []
    for nv in DIRECTIONS[:5]:
        sub = {x: R * nv[0], y: R * nv[1], z: R * nv[2]}
        for a in range(4):
            for b in range(a, 4):
                lhs = sp.expand(sp.sympify(hP3[a][b]).subs(sub, simultaneous=True)).coeff(M, pw)
                rhs = sp.expand(sp.sympify(amb[a][b]).subs(sub, simultaneous=True)).coeff(M, pw)
                eqf.append(sp.expand(lhs - rhs))
    print(f'   M^{pw}:', sp.solve(eqf, cc_, dict=True))
