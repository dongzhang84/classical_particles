"""Level 1, column by column, split into the pieces used in the hand derivation
(calc.tex Sec. 3). Independent of series.linearized_ricci: the O(M) operator is
implemented from the explicit formulas written in the text.

Background  gbar = eta + M b,  b_ab = 2 Phi delta_ab,  Phi = 1/r.
  R[h] = R0[h] + M (R1a[h] + R1b[h]) + O(M^2), with
  R0_ab  = 1/2 (d^c d_a h_bc + d^c d_b h_ac - box h_ab - d_a d_b h)
  X^c_ab = -Phi delta^{cd} P_dab - h^{cd} G_dab          (eps*M Christoffel)
  R1a_ab = d_c X^c_ab - d_b X^c_ca
  R1b_ab = G^c_cd C^d_ab + C^c_cd G^d_ab - G^c_bd C^d_ca - C^c_bd G^d_ca
  P_dab = d_a h_bd + d_b h_ad - d_d h_ab,  C^c_ab = 1/2 eta^cd P_dab,
  G_dab = delta_bd Phi_a + delta_ad Phi_b - delta_ab Phi_d,  G^c_ab = eta^cd G_dab.
Every output is linear in one constant vector V and isotropic, so it has the
form  R_00 = a V.n,  R_i0 = b V.n n_i + c V_i,  R_ij = d V.n n_i n_j + e V_(i n_j) + f V.n delta_ij.
"""
import itertools
import sympy as sp
from series import t, x, y, z, r, X, eta

M = sp.Symbol('M')
xs = (x, y, z)
n = [xi / r for xi in xs]
V = sp.symbols('V_x V_y V_z')
Vn = sum(a * b for a, b in zip(V, n))
half = sp.Rational(1, 2)
dl = lambda i, j: 1 if i == j else 0
Phi = 1 / r


def col(name):
    h = sp.zeros(4, 4)
    for i in range(3):
        for j in range(3):
            if name == 'DE':
                h[i+1, j+1] = r * half * (n[i] * V[j] + n[j] * V[i])
            if name == 'DF':
                h[i+1, j+1] = r * dl(i, j) * Vn
            if name == 'HD':
                h[i+1, j+1] = Vn * n[i] * n[j]
            if name == 'HE':
                h[i+1, j+1] = half * (V[i] * n[j] + V[j] * n[i])
            if name == 'HF':
                h[i+1, j+1] = Vn * dl(i, j)
    if name == 'DA':
        h[0, 0] = r * Vn
    if name == 'HA':
        h[0, 0] = Vn
    for i in range(3):
        v = {'DC': t * V[i], 'KB': (t / r) * Vn * n[i], 'KC': (t / r) * V[i]}.get(name, 0)
        h[0, i+1] = h[i+1, 0] = v
    return h


def d(e, a):
    return sp.diff(e, X[a])


def pieces(h):
    hu = eta * h * eta
    tr = sum(eta[a, a] * h[a, a] for a in range(4))
    R0 = sp.zeros(4, 4)
    for a, b in itertools.product(range(4), repeat=2):
        R0[a, b] = half * (sum(eta[c, c] * (d(d(h[b, c], a), c) + d(d(h[a, c], b), c)
                                            - d(d(h[a, b], c), c)) for c in range(4))
                           - d(d(tr, a), b))
    P = [[[d(h[bb, dd], a) + d(h[a, dd], bb) - d(h[a, bb], dd) for bb in range(4)]
          for a in range(4)] for dd in range(4)]          # P[d][a][b]
    Ph = [0] + [d(Phi, i) for i in (1, 2, 3)]
    G = [[[dl(bb, dd) * Ph[a] + dl(a, dd) * Ph[bb] - dl(a, bb) * Ph[dd] for bb in range(4)]
          for a in range(4)] for dd in range(4)]          # G[d][a][b]
    Cu = [[[half * eta[c, c] * P[c][a][bb] for bb in range(4)] for a in range(4)] for c in range(4)]
    Gu = [[[eta[c, c] * G[c][a][bb] for bb in range(4)] for a in range(4)] for c in range(4)]
    Xu = [[[-Phi * P[c][a][bb] - sum(hu[c, dd] * G[dd][a][bb] for dd in range(4))
            for bb in range(4)] for a in range(4)] for c in range(4)]
    R1a = sp.zeros(4, 4)
    R1b = sp.zeros(4, 4)
    for a, b in itertools.product(range(4), repeat=2):
        R1a[a, b] = (sum(d(Xu[c][a][b], c) for c in range(4))
                     - d(sum(Xu[c][c][a] for c in range(4)), b))
        R1b[a, b] = sum(Gu[c][c][dd] * Cu[dd][a][b] + Cu[c][c][dd] * Gu[dd][a][b]
                        - Gu[c][b][dd] * Cu[dd][c][a] - Cu[c][b][dd] * Gu[dd][c][a]
                        for c in range(4) for dd in range(4))
    return R0, R1a, R1b, Cu, Xu


from series import radial_coeff, rho

DIRS = [(sp.Rational(2, 7), sp.Rational(3, 7), sp.Rational(6, 7)),
        (sp.Rational(1, 3), sp.Rational(-2, 3), sp.Rational(2, 3))]
a_, b_, c_, d_, e_, f_ = sp.symbols('a b c d e f')
SEL = {'00': [(0, 0)], 'i0': [(0, 1), (0, 2), (0, 3)],
       'ij': [(i, j) for i in (1, 2, 3) for j in (1, 2, 3)]}


def powers(R, comp, tp):
    """rho powers present in the given components at t^tp, along DIRS[0]."""
    out = set()
    for k in SEL[comp]:
        e = R[k].subs({x: rho * DIRS[0][0], y: rho * DIRS[0][1], z: rho * DIRS[0][2]})
        e = sp.expand(sp.expand(sp.simplify(e)).coeff(t, tp))
        for term in sp.Add.make_args(e):
            if term != 0:
                out.add(term.as_coeff_exponent(rho)[1])
    return sorted(out)


def form(R, comp, tp, rp):
    eqs, unk = [], {'00': [a_], 'i0': [b_, c_], 'ij': [d_, e_, f_]}[comp]
    for nn in DIRS:
        vn = sum(a * b for a, b in zip(V, nn))
        for (i, j) in SEL[comp]:
            val = radial_coeff(R[i, j], nn, tp, rp)
            if comp == '00':
                ans = a_ * vn
            elif comp == 'i0':
                ans = b_ * vn * nn[j-1] + c_ * V[j-1]
            else:
                ans = (d_ * vn * nn[i-1] * nn[j-1] + e_ * half * (V[i-1] * nn[j-1] + V[j-1] * nn[i-1])
                       + f_ * vn * dl(i, j))
            q = sp.expand(ans - val)
            eqs += [q.coeff(v) for v in V]
    s = sp.solve([q for q in eqs if q != 0], unk, dict=True)
    assert s, ('no isotropic fit', comp)
    return {str(k): v for k, v in s[0].items()}


if __name__ == '__main__':
    cols = ['DA', 'DC', 'DE', 'DF', 'HA', 'HD', 'HE', 'HF', 'KB', 'KC']
    for name in cols:
        h = col(name)
        R0, R1a, R1b, Cu, Xu = pieces(h)
        print(f'=== {name} ===')
        if name.startswith('D'):
            for c, a, b in itertools.product(range(4), repeat=3):
                if a <= b and sp.simplify(Cu[c][a][b]) != 0:
                    print(f'  C^{c}_{a}{b} =', sp.simplify(Cu[c][a][b]))
            for c, a, b in itertools.product(range(4), repeat=3):
                if a <= b and sp.simplify(Xu[c][a][b]) != 0:
                    print(f'  X^{c}_{a}{b} =', sp.factor(Xu[c][a][b]))
        for lab, R in (('R0', R0), ('R1a', R1a), ('R1b', R1b)):
            for comp in ('00', 'i0', 'ij'):
                for tp in (0, 1):
                    for rp in powers(R, comp, tp):
                        print(f'  {lab}_{comp} t^{tp} r^{rp}:', form(R, comp, tp, rp))
