"""Level 3 (q = 0, background field F != 0): the near-zone first-order Maxwell
field and its stress tensor (calc.tex Sec. 5).

Near-zone background: Schwarzschild in isotropic coordinates,
    ds^2 = -alpha^2 dt^2 + psi^4 dx^i dx^i,  psi = 1 + M/2r,  alpha = (1 - M/2r)/psi.
Near-zone first-order potential (uniform field E_i = F_i0, B_k = eps_kij F_ij / 2 at large r):
    A_0 = (1 - M/2r)^2 E.x,     A_i = (1/2) psi^4 (B x x)_i.
Checks:
  1. Maxwell's equations hold exactly.
  2. The O(M) part equals the r^0 term of the far-zone singular field,
     b^alpha = -M F_{alpha c} n^c  (Zimmerman & Poisson 2014, from U^alpha_{|beta gamma mu}).
  3. 8 pi T_ab = 2 F_ac F_b^c - (1/2) g_ab F^2 (exact):
       T_00, T_ij contain only EE and BB terms, T_0i only EB terms;
       T_00, T_ij are even functions of x  ->  no electric l = 1 part;
       T_0i is an even function of x       ->  its l = 1 part is listed, order by order in M.
"""
import itertools
import sympy as sp

t, x, y, z = sp.symbols('t x y z', real=True)
X = (t, x, y, z)
xs = (x, y, z)
M, R = sp.symbols('M R', positive=True)      # R stands for r = |x|
E = sp.symbols('E_x E_y E_z')
B = sp.symbols('B_x B_y B_z')
lev = sp.LeviCivita


def d(f, a):
    """Partial derivative in (t, x, y, z), treating R as |x|."""
    if a == 0:
        return sp.diff(f, t)
    return sp.diff(f, X[a]) + sp.diff(f, R) * X[a] / R


r = R
psi = 1 + M / (2 * r)
alpha = (1 - M / (2 * r)) / psi
g = sp.diag(-alpha**2, psi**4, psi**4, psi**4)
ginv = sp.diag(-1 / alpha**2, psi**-4, psi**-4, psi**-4)
sqrtg = alpha * psi**6

Ex = sum(E[i] * xs[i] for i in range(3))
BxX = [sum(lev(i, j, k) * B[j] * xs[k] for j in range(3) for k in range(3)) for i in range(3)]
A = [(1 - M / (2 * r))**2 * Ex] + [sp.Rational(1, 2) * psi**4 * BxX[i] for i in range(3)]
F = sp.Matrix(4, 4, lambda a, b: d(A[b], a) - d(A[a], b))
Fup = ginv * F * ginv
def onshell(e):
    """e reduced modulo x^2 + y^2 + z^2 = R^2 (zero iff e vanishes for R = |x|)."""
    num = sp.expand(sp.numer(sp.together(sp.expand(e))))
    return sp.expand(sp.rem(num, x**2 + y**2 + z**2 - R**2, x))

# 1. Maxwell
for a in range(4):
    div = sum(d(sqrtg * Fup[a, b], b) for b in range(4))
    assert onshell(div) == 0, ('Maxwell', a)
print('1. Maxwell equations hold exactly')

# 2. O(M) part of A versus the singular-field r^0 term
n = [xi / R for xi in xs]
Fbg = sp.zeros(4, 4)
for i in range(3):
    Fbg[i+1, 0], Fbg[0, i+1] = E[i], -E[i]
    for j in range(3):
        Fbg[i+1, j+1] = sum(lev(i, j, k) * B[k] for k in range(3))
A_M = [sp.diff(a, M).subs(M, 0) for a in A]
b_up = [-sum(Fbg[al, c+1] * n[c] for c in range(3)) for al in range(4)]   # b^alpha
b_dn = [-b_up[0]] + b_up[1:]
assert all(onshell(A_M[al] - b_dn[al]) == 0 for al in range(4))
print('2. O(M) part of A:', [sp.simplify(v) for v in A_M], '= singular-field r^0 term')

# 3. stress tensor
F2 = sum(F[a, b] * Fup[a, b] for a in range(4) for b in range(4))
T8 = sp.Matrix(4, 4, lambda a, b: 2 * sum(F[a, c] * (F * ginv)[b, c] for c in range(4)) - g[a, b] * F2 / 2)

subsE = {Bi: 0 for Bi in B}
subsB = {Ei: 0 for Ei in E}
scale = sp.Symbol('s')
for (a, b) in [(0, 0), (1, 1), (1, 2), (0, 1), (0, 2)]:
    e = sp.factor(T8[a, b])
    # bilinear EB part: scale E -> s E and keep s^1
    eb = sp.expand(sp.numer(sp.together(e)).subs({Ei: scale * Ei for Ei in E})).coeff(scale, 1)
    if a == 0 and b == 0 or (a > 0 and b > 0):
        assert onshell(eb) == 0, ('EB in', a, b)
    else:
        assert onshell(e.subs(subsE)) == 0 and onshell(e.subs(subsB)) == 0, ('EE/BB in', a, b)
    par = onshell(e.subs({x: -x, y: -y, z: -z}, simultaneous=True) - e)
    assert par == 0, ('parity', a, b)
print('3. T_00, T_ij: only EE, BB; T_0i: only EB; all components even in x')

# l = 1 part of T_0i, order by order in M: V_i(n) -> (QB.n) n_i + QC_i, from
#   <V_i> = QB/3 + QC,   <(V.n) n_i> = (QB + QC)/3   (sphere averages)
nn = sp.symbols('n_x n_y n_z')


def sphere_avg(poly):
    poly = sp.Poly(sp.expand(poly), *nn)
    tot = 0
    for (p, q, s), c in poly.terms():
        if p % 2 or q % 2 or s % 2:
            continue
        df = lambda k: sp.factorial2(k - 1) if k > 0 else 1
        tot += c * df(p) * df(q) * df(s) / sp.factorial2(p + q + s + 1)
    return sp.simplify(tot)


ExB = [sum(lev(i, j, k) * E[j] * B[k] for j in range(3) for k in range(3)) for i in range(3)]
T0 = [sp.factor(T8[0, i+1]) for i in range(3)]
rho = R
T0n = [sp.series(Ti.subs({x: rho * nn[0], y: rho * nn[1], z: rho * nn[2]}), M, 0, 3).removeO() for Ti in T0]
for p in range(3):
    V = [sp.simplify(sp.expand(Ti).coeff(M, p) * rho**p) for Ti in T0n]
    avg = [sphere_avg(v) for v in V]
    vn = sum(V[i] * nn[i] for i in range(3))
    avgn = [sphere_avg(vn * nn[i]) for i in range(3)]
    QB = [sp.simplify(sp.Rational(3, 2) * (3 * avgn[i] - avg[i])) for i in range(3)]
    QC = [sp.simplify(avg[i] - QB[i] / 3) for i in range(3)]
    cB = sp.simplify(QB[0] / ExB[0]) if QB[0] != 0 else 0
    cC = sp.simplify(QC[0] / ExB[0]) if QC[0] != 0 else 0
    assert all(sp.simplify(QB[i] - cB * ExB[i]) == 0 and sp.simplify(QC[i] - cC * ExB[i]) == 0 for i in range(3))
    print(f'   8 pi T_0i, M^{p}/r^{p} term, l = 1 part: ({cB}) (ExB).n n_i + ({cC}) (ExB)_i')

# 4. the component formulas written in calc.tex, eq. (T2comp)
cE = [F[i+1, 0] for i in range(3)]
cF = [[F[i+1, j+1] for j in range(3)] for i in range(3)]
E2 = sum(v**2 for v in cE)
FF = sum(cF[i][j]**2 for i in range(3) for j in range(3))
Finv2 = -2 * cE[0]**0 * E2 / (alpha**2 * psi**4) + FF / psi**8
assert onshell(T8[0, 0] - (E2 / psi**4 + alpha**2 * FF / (2 * psi**8))) == 0
for i in range(3):
    assert onshell(T8[0, i+1] + 2 * sum(cE[j] * cF[i][j] for j in range(3)) / psi**4) == 0
    for j in range(3):
        rhs = (-2 * cE[i] * cE[j] / alpha**2 + 2 * sum(cF[i][k] * cF[j][k] for k in range(3)) / psi**4
               - (psi**4 / 2) * (1 if i == j else 0) * Finv2)
        assert onshell(T8[i+1, j+1] - rhs) == 0
print('4. component formulas for 8 pi T^(2) confirmed')
