"""Level 4, near-zone second order: the bilinear source from products of two
first-order near-zone fields, projected on the l = 1 rows.

mu1: the constant-metric field c_ab (from the far-zone regular field and the
     shift), extended into the near zone as the gauge transform of RN by
     xi^0 = -c_0i x^i - c_00 t/2, xi^i = c_ij x^j / 2  (static, body at rest).
mu2: the uniform-field / acceleration field solved in level4_l1.py, with its
     free (gauge) parameters kept symbolic.
"""
import sys, time
import sympy as sp
import trunc as T
from trunc import t, x, y, z, R, nn
from level4_l1 import M, q, E, B, xs, n, dot, cross, half, dl, build

g_ = sp.symbols('hC_1:4'); hE = sp.symbols('hE_1:4'); hF = sp.symbols('hF_1:4')
hMij = sp.symbols('hMij_1:4'); aC = sp.symbols('aC_1:4')
c = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'c{min(i,j)}{max(i,j)}'))
a = [q * Ei / M for Ei in E]


def field2(which):
    """First-order field solved in level4_l1 (a = qE/M imposed)."""
    h = [[0] * 4 for _ in range(4)]
    eh = [[0] * 4 for _ in range(4)]
    A0 = [0] * 4
    eA = [0] * 4
    if 'E' in which:
        h[0][0] = -2 * dot(a, xs)
        eh[0][0] = 5 * q * dot(E, n)
        for i in range(3):
            for j in range(3):
                eh[i+1][j+1] = ((q * dot(E, n) - dot(hF, n)) * n[i] * n[j]
                                + half * (hE[i] * n[j] + hE[j] * n[i]) + dot(hF, n) * dl(i, j))
        A0[0] = dot(E, xs)
        eA[0] = -(2 * M**2 + q**2) / (2 * M) * dot(E, n)
    if 'B' in which:
        for i in range(3):
            eh[0][i+1] = eh[i+1][0] = 2 * q * cross(B, n)[i]
            A0[i+1] = half * cross(B, xs)[i]
            eA[i+1] = M * cross(B, n)[i]
    if 'G' in which:   # free gauge pieces found in level4_l1
        for i in range(3):
            eh[0][i+1] = eh[i+1][0] = eh[0][i+1] + g_[i]
            for j in range(3):
                w = cross(hMij, n)
                eh[i+1][j+1] += half * (n[i] * w[j] + n[j] * w[i])
            eA[i+1] += aC[i]
    hh = T.mat(lambda i, j: T.S({(0, 0, 1): h[i][j], (1, 0, 1): eh[i][j]}))
    AA = [T.S({(0, 0, 1): A0[i], (1, 0, 1): eA[i]}) for i in range(4)]
    return hh, AA


def field1():
    """Gauge transform of the O(eps) background by the linear xi."""
    xi = [-sum(c[0, i+1] * xs[i] for i in range(3)) - half * c[0, 0] * t] + \
         [half * sum(c[i+1, j+1] * xs[j] for j in range(3)) for i in range(3)]
    gb = [[(T.eta[i, j] if True else 0) for j in range(4)] for i in range(4)]
    bb = [[(2 * M / R if i == j else 0) for j in range(4)] for i in range(4)]
    Ab = [-q / R, 0, 0, 0]
    de = T.d_expr

    def lie_g(g):
        return [[sp.expand(sum(xi[k] * de(g[i][j], k) + g[k][j] * de(xi[k], i) + g[i][k] * de(xi[k], j)
                               for k in range(4))) for j in range(4)] for i in range(4)]
    h0 = lie_g(gb)
    h1 = lie_g(bb)
    a1 = [sp.expand(sum(xi[k] * de(Ab[i], k) + Ab[k] * de(xi[k], i) for k in range(4))) for i in range(4)]
    for i in range(4):
        for j in range(4):
            assert sp.expand(h0[i][j] - c[i, j]) == 0
    hh = T.mat(lambda i, j: T.S({(0, 1, 0): h0[i][j], (1, 1, 0): h1[i][j]}))
    AA = [T.S({(1, 1, 0): a1[i]}) for i in range(4)]
    return hh, AA


def project(Ein, Mx, key, rpow=-2):
    get = lambda e: T.radial_part(e.get(key, 0), rpow)
    S00 = get(Ein[0][0])
    V0i = [get(Ein[0][i+1]) for i in range(3)]
    Tij = [[get(Ein[i+1][j+1]) for j in range(3)] for i in range(3)]
    M0 = get(Mx[0])
    Mi = [get(Mx[i+1]) for i in range(3)]
    out = {'A': T.l1_scalar(S00)}
    out['D'], out['E'], out['F'] = T.l1_tensor(Tij)
    out['B0i_static'], out['C0i_static'] = T.l1_vector(V0i)
    out['X'] = T.l1_scalar(M0)
    out['Mi_P'], out['Mi_C'] = T.l1_vector(Mi)
    out['mag_0i'] = T.l1_magnetic_vector(V0i)
    return out


if __name__ == '__main__':
    which = sys.argv[1] if len(sys.argv) > 1 else 'EBG'
    T.set_limits(1, 1, 1)
    h1, A1 = field1()
    h2, A2 = field2(which)
    hsum = T.mat(lambda i, j: T.add(h1[i][j], h2[i][j]))
    Asum = [T.add(A1[i], A2[i]) for i in range(4)]
    g, A = build(hsum, Asum)
    t0 = time.time()
    Ein, Mx = T.einstein_maxwell(g, A)
    print('computed', round(time.time() - t0, 1), 's', flush=True)
    # first-order residuals must vanish (both fields solve the O(eps) equations)
    for key in [(0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1)]:
        bad = [sp.simplify(Ein[i][j].get(key, 0)) for i in range(4) for j in range(4)] + \
              [sp.simplify(Mx[i].get(key, 0)) for i in range(4)]
        print(key, 'residual zero:', all(b == 0 for b in bad), flush=True)
    for rp in (-2, -1, -3):
        out = project(Ein, Mx, (1, 1, 1), rp)
        print(f'--- (eps mu1 mu2), R^{rp} ---')
        for k, v in out.items():
            vv = [sp.factor(e) for e in v]
            if any(e != 0 for e in vv):
                print(f'  {k}:', vv)
