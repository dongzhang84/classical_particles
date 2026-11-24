"""Level 2: charged body, uncharged vacuum background (q != 0, F = 0).

Near zone through O(lambda^2), l = 1, time-reversal-even sector, as in GW
Sec. VI, now with the electromagnetic potential.

Background (near zone, lambda^0), kept to the orders that matter at 1/r^2:
    g = eta + (2M/r) delta            (linearized RN; q^2/r^2 terms enter at 1/r^3)
    A = -(q/r) dt                     (Coulomb; the qM/r^2 term enters at 1/r^3)
Second-order near-zone perturbation (l = 1, time-even):
  metric: as in level 1 (DA, DC, DE, DF, HA, HD, HE, HF, KB, KC);
  potential:
    a_0 = E0 . x + HP . n
    a_i = t EC_i + (t/r) (KBp . n n_i + KCp_i)
  where E0_i = d_i a^tail_0 and EC_i = d_0 a^tail_i on the worldline (the far-zone
  tail field, so the tail electric field is f_0i = EC_i - E0_i), and HP, KBp, KCp
  are unknown homogeneous pieces (the EM analogues of H, K).
Field equations at this order:
  Einstein:  R_ab - 8 pi T_ab = 0, with T_ab linearized in eps at O(q):
             8 pi T_ab -> 2 (F0_ac f_b^c + f_ac F0_b^c - (1/2) eta_ab F0_cd f^cd)
  Maxwell:   nabla_b F^ab = 0, linearized in eps, kept at O(1), O(M), O(q).
Rows extracted (time-even l = 1):
  Einstein R_00 |1/r^2 -> A ; R_i0 |t/r^3 -> B, C ; R_ij |1/r^2 -> D, E, F
  Maxwell  M^0  |1/r^2 -> X ; M^i  |t/r^3 -> Y n n_i + Z_i
Then: left null space of all unknown columns (H, K, HP, KBp, KCp) gives the
universal relations among (DA, DC, DE, DF, E0, EC).
"""
import itertools, time, sys
import sympy as sp
from series import (t, x, y, z, r, X, eta, linearized_ricci, DIRECTIONS,
                    radial_coeff)

M, q = sp.symbols('M q', positive=True)
xs = (x, y, z)
nvec = [xi / r for xi in xs]
half = sp.Rational(1, 2)

def vec(name):
    return sp.symbols(f'{name}_x {name}_y {name}_z')

knames = ['DA', 'DC', 'DE', 'DF', 'E0', 'EC']          # known (far-zone data)
unames = ['HA', 'HD', 'HE', 'HF', 'KB', 'KC', 'HP', 'KBp', 'KCp']  # unknown
V = {k: vec(k) for k in knames + unames}
dot = lambda a, b: sum(ai * bi for ai, bi in zip(a, b))
delta = lambda i, j: 1 if i == j else 0


def metric_pert(V):
    h = sp.zeros(4, 4)
    h[0, 0] = r * dot(V['DA'], nvec) + dot(V['HA'], nvec)
    for i in range(3):
        h[0, i+1] = h[i+1, 0] = t * V['DC'][i] + (t / r) * (dot(V['KB'], nvec) * nvec[i] + V['KC'][i])
        for j in range(3):
            h[i+1, j+1] = (r * (half * (nvec[i] * V['DE'][j] + nvec[j] * V['DE'][i]) + delta(i, j) * dot(V['DF'], nvec))
                           + dot(V['HD'], nvec) * nvec[i] * nvec[j]
                           + half * (V['HE'][i] * nvec[j] + V['HE'][j] * nvec[i])
                           + dot(V['HF'], nvec) * delta(i, j))
    return h


def pot_pert(V):
    a = [0] * 4
    a[0] = dot(V['E0'], xs) + dot(V['HP'], nvec)
    for i in range(3):
        a[i+1] = t * V['EC'][i] + (t / r) * (dot(V['KBp'], nvec) * nvec[i] + V['KCp'][i])
    return a


def field(a):
    return sp.Matrix(4, 4, lambda m, n: sp.diff(a[n], X[m]) - sp.diff(a[m], X[n]))


def einstein_rows(h, a):
    b1 = sp.diag(2 / r, 2 / r, 2 / r, 2 / r)
    R0, R1 = linearized_ricci(b1, h, M)
    Ric = R0 + M * R1
    F0 = field([-q / r, 0, 0, 0])
    f = field(a)
    ei = eta  # raise with eta: enough at O(eps q)
    F0u = ei * F0 * ei   # F0^{ab}
    F0m = F0 * ei        # F0_a^b
    fm = f * ei
    F0f = sum(F0[c, d] * (ei * f * ei)[c, d] for c in range(4) for d in range(4))
    T8 = sp.Matrix(4, 4, lambda m, n: 2 * (sum(F0[m, c] * fm[n, c] for c in range(4))
                                           + sum(f[m, c] * F0m[n, c] for c in range(4))
                                           - half * eta[m, n] * F0f))
    return Ric - T8


def maxwell_rows(h, a):
    """eps-linear part of nabla_b F^{ab} = (-g)^{-1/2} d_b[(-g)^{1/2} F^{ab}],
    kept at O(1), O(M), O(q) (O(qM) and higher enter at 1/r^3).

      O(1): d_b f^{ab}, indices raised with eta
      O(M): from f with g = eta + M b1:  d_b[ sqrt(-g) g^ac g^bd f_cd ] at O(M);
            the (-g)^{-1/2} prefactor multiplies d_b f^{ab}, which is already
            O(1/r^2) for the H', K' pieces, so its O(M) part is O(1/r^3).
      O(q): linearized Maxwell operator (C1 Eq. (3)) with background Coulomb F0:
            d_b Y^{ab},  Y^{ab} = (1/2) h F0^{ab} - h^{ac} F0_c^b - h^{bd} F0^a_d
    """
    b1 = sp.diag(2 / r, 2 / r, 2 / r, 2 / r)
    f = field(a)
    fu = eta * f * eta
    fM = -(eta * b1 * eta * f * eta + eta * f * eta * b1 * eta) \
         + sp.Rational(1, 2) * sum(eta[c, c] * b1[c, c] for c in range(4)) * fu
    F0 = field([-q / r, 0, 0, 0])
    F0u = eta * F0 * eta
    hu = eta * h * eta
    trh = sum(eta[c, c] * h[c, c] for c in range(4))
    Y = sp.Matrix(4, 4, lambda m, n: half * trh * F0u[m, n]
                  - sum(hu[m, c] * (F0 * eta)[c, n] for c in range(4))
                  - sum(hu[n, d] * (eta * F0)[m, d] for d in range(4)))
    out = []
    for m in range(4):
        expr = sum(sp.diff(fu[m, n] + M * fM[m, n] + Y[m, n], X[n]) for n in range(4))
        out.append(expr)
    return out


def fit(values_by_dir, ansatz_fn, unknowns):
    eqs = [sp.expand(ansatz_fn(k) - v) for k, v in values_by_dir.items()]
    sol = sp.solve(eqs, unknowns, dict=True)
    assert sol, 'not pure l=1'
    sol = sol[0]
    assert all(sp.simplify(e.subs(sol)) == 0 for e in eqs)
    return sol


def extract(V):
    h = metric_pert(V)
    a = pot_pert(V)
    t0 = time.time()
    E = einstein_rows(h, a)
    print(f'Einstein done {time.time()-t0:.0f}s', flush=True)
    Mx = maxwell_rows(h, a)
    print(f'Maxwell done {time.time()-t0:.0f}s', flush=True)
    dirs = DIRECTIONS
    R = {k: vec('R' + k) for k in 'ABCDEFXYZ'}
    sol = {}
    sol.update(fit({n: radial_coeff(E[0, 0], n, 0, -2) for n in dirs},
                   lambda n: dot(R['A'], n), list(R['A'])))
    vij = {(n, i, j): radial_coeff(E[i+1, j+1], n, 0, -2) for n in dirs for i in range(3) for j in range(i, 3)}
    sol.update(fit(vij, lambda k: dot(R['D'], k[0]) * k[0][k[1]] * k[0][k[2]]
                   + half * (R['E'][k[1]] * k[0][k[2]] + R['E'][k[2]] * k[0][k[1]])
                   + dot(R['F'], k[0]) * delta(k[1], k[2]),
                   list(R['D']) + list(R['E']) + list(R['F'])))
    vi0 = {(n, i): radial_coeff(E[0, i+1], n, 1, -3) for n in dirs for i in range(3)}
    sol.update(fit(vi0, lambda k: dot(R['B'], k[0]) * k[0][k[1]] + R['C'][k[1]],
                   list(R['B']) + list(R['C'])))
    sol.update(fit({n: radial_coeff(Mx[0], n, 0, -2) for n in dirs},
                   lambda n: dot(R['X'], n), list(R['X'])))
    vmi = {(n, i): radial_coeff(Mx[i+1], n, 1, -3) for n in dirs for i in range(3)}
    sol.update(fit(vmi, lambda k: dot(R['Y'], k[0]) * k[0][k[1]] + R['Z'][k[1]],
                   list(R['Y']) + list(R['Z'])))
    return R, sol


if __name__ == '__main__':
    R, sol = extract(V)
    cols = knames + unames
    rows = list('ABCDEFXYZ')
    mat = []
    for rn in rows:
        e = sp.expand(sol[R[rn][0]])
        for k in cols:
            assert e.coeff(V[k][1]) == 0 and e.coeff(V[k][2]) == 0, (rn, k)
        mat.append([sp.factor(e.coeff(V[k][0])) for k in cols])
    Mt = sp.Matrix(mat)
    print('columns:', cols)
    for rn, row in zip(rows, mat):
        print(rn, row)
    nk = len(knames)
    U = Mt[:, nk:]
    print('\nuniversal relations (coefficients on', knames, '):')
    for c in (U.T).nullspace():
        rel = [sp.factor(v) for v in (c.T * Mt[:, :nk])]
        if any(v != 0 for v in rel):
            print('  rows', [sp.nsimplify(v) for v in c.T], '->', rel)
    sp.srepr(Mt)
    with open('level2_matrix.txt', 'w') as fh:
        fh.write(str(cols) + '\n')
        for rn, row in zip(rows, mat):
            fh.write(rn + ' ' + str(row) + '\n')
