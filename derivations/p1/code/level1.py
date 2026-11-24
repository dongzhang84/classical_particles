"""Level 1: reproduce GW 2008 Eq. (106) (vacuum, q = 0, F = 0, no spin).

Near-zone background: linearized Schwarzschild in isotropic coordinates,
    gbar = eta + (2M/r) delta_ab            (GW Eq. (86), O(M) part; the M^2/r^2
                                             part only affects order 1/r^3).
Second-order near-zone perturbation, l = 1, electric parity, even under time
reversal (GW Eqs. (98)-(101)):
    h_00 = r DA.n + HA.n
    h_i0 = t DC_i + (t/r) (KB.n n_i + KC_i)
    h_ij = r ( n_(i DE_j) + delta_ij DF.n ) + HD.n n_i n_j + HE_(i n_j) + HF.n delta_ij
Compute the linearized Ricci tensor and extract (GW Eqs. (94)-(96))
    R_00 |_{1/r^2}   = RA.n
    R_ij |_{1/r^2}   = RD.n n_i n_j + RE_(i n_j) + RF.n delta_ij
    R_i0 |_{t/r^3}   = RB.n n_i + RC_i
Then compare with the matrix in GW Eq. (106) (spin column dropped).
"""
import sys, time
import sympy as sp
from series import (t, x, y, z, r, eta, linearized_ricci, DIRECTIONS,
                    radial_coeff, fit_l1)

M = sp.Symbol('M', positive=True)
xs = (x, y, z)
nvec = [xi / r for xi in xs]

def vec(name):
    return sp.symbols(f'{name}_x {name}_y {name}_z')

DA, DC, DE, DF = vec('DA'), vec('DC'), vec('DE'), vec('DF')
HA, HD, HE, HF = vec('HA'), vec('HD'), vec('HE'), vec('HF')
KB, KC = vec('KB'), vec('KC')
inputs = [DA, DC, DE, DF, HA, HD, HE, HF, KB, KC]
labels = ['DA', 'DC', 'DE', 'DF', 'HA', 'HD', 'HE', 'HF', 'KB', 'KC']

dot = lambda a, b: sum(ai * bi for ai, bi in zip(a, b))
delta = lambda i, j: 1 if i == j else 0

h = sp.zeros(4, 4)
h[0, 0] = r * dot(DA, nvec) + dot(HA, nvec)
for i in range(3):
    h0i = t * DC[i] + (t / r) * (dot(KB, nvec) * nvec[i] + KC[i])
    h[0, i + 1] = h[i + 1, 0] = h0i
    for j in range(3):
        hij = (r * (sp.Rational(1, 2) * (nvec[i] * DE[j] + nvec[j] * DE[i])
                    + delta(i, j) * dot(DF, nvec))
               + dot(HD, nvec) * nvec[i] * nvec[j]
               + sp.Rational(1, 2) * (HE[i] * nvec[j] + HE[j] * nvec[i])
               + dot(HF, nvec) * delta(i, j))
        h[i + 1, j + 1] = hij

b1 = sp.diag(2 / r, 2 / r, 2 / r, 2 / r)

t0 = time.time()
R0, R1 = linearized_ricci(b1, h, M)
R = R0 + M * R1
print(f'linearized Ricci done in {time.time()-t0:.0f}s', flush=True)

dirs = DIRECTIONS
# R_00 at t^0 rho^-2
v00 = {n: radial_coeff(R[0, 0], n, 0, -2) for n in dirs}
RA = vec('RA')
solA, resA = fit_l1(v00, lambda n: dot(RA, n), list(RA))
# R_ij at t^0 rho^-2
RD, RE, RF = vec('RD'), vec('RE'), vec('RF')
vij = {}
for n in dirs:
    for i in range(3):
        for j in range(i, 3):
            vij[(n, i, j)] = radial_coeff(R[i + 1, j + 1], n, 0, -2)
def ans_ij(key):
    n, i, j = key
    return (dot(RD, n) * n[i] * n[j] + sp.Rational(1, 2) * (RE[i] * n[j] + RE[j] * n[i])
            + dot(RF, n) * delta(i, j))
eqs = [sp.expand(ans_ij(k) - v) for k, v in vij.items()]
solDEF = sp.solve(eqs, list(RD) + list(RE) + list(RF), dict=True)
assert solDEF, 'R_ij not pure l=1'
solDEF = solDEF[0]
resDEF = [sp.simplify(e.subs(solDEF)) for e in eqs]
# R_i0 at t^1 rho^-3
RB, RC = vec('RB'), vec('RC')
vi0 = {(n, i): radial_coeff(R[0, i + 1], n, 1, -3) for n in dirs for i in range(3)}
eqs = [sp.expand(dot(RB, k[0]) * k[0][k[1]] + RC[k[1]] - v) for k, v in vi0.items()]
solBC = sp.solve(eqs, list(RB) + list(RC), dict=True)
assert solBC, 'R_i0 not pure l=1'
solBC = solBC[0]
resBC = [sp.simplify(e.subs(solBC)) for e in eqs]

assert all(e == 0 for e in resA + resDEF + resBC), 'residuals nonzero'
sol = {**solA, **solDEF, **solBC}

# matrix: row = output vector (x-component), column = input vector
rows = [('RA', RA), ('RB', RB), ('RC', RC), ('RD', RD), ('RE', RE), ('RF', RF)]
print('\nRows: R^A..R^F ; columns:', ' '.join(labels))
mat = []
for name, Rv in rows:
    e = sp.expand(sol[Rv[0]])
    row = [sp.factor(e.coeff(v[0])) for v in inputs]
    # isotropy check: no cross-component mixing
    for v in inputs:
        assert e.coeff(v[1]) == 0 and e.coeff(v[2]) == 0
    mat.append(row)
    print(name, row)

# GW Eq. (106), spin column dropped, overall factor -1/2
gw = sp.Matrix([
    [-3*M, 0, -M, -M, -2, 0, 0, 0, -2, 2],
    [0, -6*M, 0, 0, 0, 0, 0, 0, -3, -3],
    [0, 2*M, 0, 0, 0, 0, 0, 0, 1, 1],
    [-6*M, 0, 0, -6*M, -3, 3, 0, 3, -6, 0],
    [0, 4*M, 0, 4*M, 2, -2, 0, -2, 2, -2],
    [0, 3*M, -3*M, 9*M, 1, -3, 0, -3, 2, 0],
]) * sp.Rational(-1, 2)
ours = sp.Matrix(mat)
print('\nGW (106) x (-1/2):')
sp.pprint(gw)
print('\nDifference (ours - GW):')
sp.pprint(sp.simplify(ours - gw))
