"""Level 2, column by column: the electromagnetic pieces (calc.tex Sec. 4).

Independent of level2.py: every piece is obtained by series-expanding the exact
expressions, not from the linearized formulas used there.

  Maxwell  M^a = d_b[ sqrt(-g) g^{ac} g^{bd} F_cd ]   (the (-g)^{-1/2} prefactor
           only affects 1/r^3), with g = eta + M b + eps h, F = F0 + eps f,
           F0 from A0 = -q/r dt.  Pieces kept:
             flat : eps^1 M^0 q^0   (from f)
             O(M) : eps^1 M^1 q^0   (from f in the background b)
             O(q) : eps^1 M^0 q^1   (from h acting on F0)
  Einstein source  8 pi T_ab = 2 F_ac F_b^c - (1/2) g_ab F^2 at eps^1 q^1 M^0
           (indices with g = eta + eps h; the h F0 F0 part is O(q^2)).
Rows: M^0 |1/r^2 -> X ;  M^i |t/r^3 -> Y n n_i + Z_i ;
      -8 pi T_00 |1/r^2 -> A ; -8 pi T_ij |1/r^2 -> D, E, F  (sign: rows are R - 8 pi T).
"""
import itertools
import sympy as sp
from series import t, x, y, z, r, X, eta, radial_coeff

M, q, eps = sp.symbols('M q epsilon')
xs = (x, y, z)
n = [xi / r for xi in xs]
V = sp.symbols('V_x V_y V_z')
dot = lambda a, b: sum(ai * bi for ai, bi in zip(a, b))
Vn = dot(V, n)
half = sp.Rational(1, 2)
dl = lambda i, j: 1 if i == j else 0
DIRS = [(sp.Rational(2, 7), sp.Rational(3, 7), sp.Rational(6, 7)),
        (sp.Rational(1, 3), sp.Rational(-2, 3), sp.Rational(2, 3))]


def metric_col(name):
    h = sp.zeros(4, 4)
    if name == 'DA':
        h[0, 0] = r * Vn
    for i in range(3):
        if name == 'DC':
            h[0, i+1] = h[i+1, 0] = t * V[i]
        for j in range(3):
            if name == 'DE':
                h[i+1, j+1] = r * half * (n[i] * V[j] + n[j] * V[i])
            if name == 'DF':
                h[i+1, j+1] = r * dl(i, j) * Vn
    return h


def pot_col(name):
    a = [0] * 4
    if name == 'E0':
        a[0] = dot(V, xs)
    if name == 'HP':
        a[0] = Vn
    for i in range(3):
        a[i+1] = {'EC': t * V[i], 'KBp': (t / r) * Vn * n[i], 'KCp': (t / r) * V[i]}.get(name, 0)
    return a


def field(a):
    return sp.Matrix(4, 4, lambda m, k: sp.diff(a[k], X[m]) - sp.diff(a[m], X[k]))


F0 = field([-q / r, 0, 0, 0])
b = sp.diag(2 / r, 2 / r, 2 / r, 2 / r)


def maxwell(h, a):
    g = eta + M * b + eps * h
    # inverse and sqrt(-g) to first order in M and eps, including M*eps
    ginv = (eta - M * eta * b * eta - eps * eta * h * eta
            + M * eps * (eta * b * eta * h * eta + eta * h * eta * b * eta))
    trb = sum(eta[c, c] * b[c, c] for c in range(4))
    trh = sum(eta[c, c] * h[c, c] for c in range(4))
    trbh = sum((eta * b * eta)[c, d] * h[c, d] for c in range(4) for d in range(4))
    sq = 1 + half * (M * trb + eps * trh) + M * eps * (half * half * trb * trh - half * trbh)
    Ftot = F0 + eps * field(a)
    Fu = ginv * Ftot * ginv.T
    out = []
    for m in range(4):
        e = sp.expand(sum(sp.diff(sq * Fu[m, k], X[k]) for k in range(4)))
        e1 = e.coeff(eps, 1)
        out.append({'flat': e1.coeff(M, 0).coeff(q, 0),
                    'M': e1.coeff(M, 1).coeff(q, 0),
                    'q': e1.coeff(M, 0).coeff(q, 1)})
    return out


def einstein_source(h, a):
    g = eta + eps * h
    ginv = eta - eps * eta * h * eta
    Ftot = F0 + eps * field(a)
    Fmix = Ftot * ginv                      # F_a^b
    F2 = sum(Ftot[c, d] * (ginv * Ftot * ginv.T)[c, d] for c in range(4) for d in range(4))
    T8 = sp.Matrix(4, 4, lambda m, k: 2 * sum(Ftot[m, c] * Fmix[k, c] for c in range(4)) - half * g[m, k] * F2)
    return T8.applyfunc(lambda e: sp.expand(sp.expand(e).coeff(eps, 1)).coeff(q, 1))


a_, b_, c_, d_, e_, f_ = sp.symbols('a b c d e f')


def fit_scalar(expr, tp, rp):
    eqs = []
    for nn in DIRS:
        val = radial_coeff(expr, nn, tp, rp)
        eqs += [sp.expand(a_ * dot(V, nn) - val).coeff(v) for v in V]
    s = sp.solve(eqs, [a_], dict=True)
    assert s, 'scalar fit'
    return s[0][a_]


def fit_vector(comps, tp, rp):
    eqs = []
    for nn in DIRS:
        for i in range(3):
            val = radial_coeff(comps[i], nn, tp, rp)
            eqs += [sp.expand(b_ * dot(V, nn) * nn[i] + c_ * V[i] - val).coeff(v) for v in V]
    s = sp.solve(eqs, [b_, c_], dict=True)
    assert s, 'vector fit'
    return s[0][b_], s[0][c_]


def fit_tensor(T, tp, rp):
    eqs = []
    for nn in DIRS:
        for i in range(3):
            for j in range(3):
                val = radial_coeff(T[i+1, j+1], nn, tp, rp)
                ans = (d_ * dot(V, nn) * nn[i] * nn[j] + e_ * half * (V[i] * nn[j] + V[j] * nn[i])
                       + f_ * dot(V, nn) * dl(i, j))
                eqs += [sp.expand(ans - val).coeff(v) for v in V]
    s = sp.solve(eqs, [d_, e_, f_], dict=True)
    assert s, 'tensor fit'
    return s[0][d_], s[0][e_], s[0][f_]


if __name__ == '__main__':
    zero_a = [0, 0, 0, 0]
    for name in ['DA', 'DC', 'DE', 'DF', 'E0', 'EC', 'HP', 'KBp', 'KCp']:
        h = metric_col(name) if name.startswith('D') else sp.zeros(4, 4)
        a = pot_col(name) if not name.startswith('D') else zero_a
        mx = maxwell(h, a)
        print(f'=== {name} ===')
        for part in ('flat', 'M', 'q'):
            X0 = fit_scalar(mx[0][part], 0, -2)
            Y, Z = fit_vector([mx[i][part] for i in (1, 2, 3)], 1, -3)
            if any(v != 0 for v in (X0, Y, Z)):
                print(f'  Maxwell {part:4s}: X={X0}  Y={Y}  Z={Z}')
        if not name.startswith('D'):
            T8 = einstein_source(h, a)
            A = -fit_scalar(T8[0, 0], 0, -2)
            Dd, Ee, Ff = [-v for v in fit_tensor(T8, 0, -2)]
            B, C = [-v for v in fit_vector([T8[0, i] for i in (1, 2, 3)], 1, -3)]
            print(f'  -8piT: A={A} B={B} C={C} D={Dd} E={Ee} F={Ff}')
