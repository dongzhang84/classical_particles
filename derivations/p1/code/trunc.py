"""Truncated multivariate perturbation series for the level 4 calculation.

A quantity is a dict {(i, j, k): expr}, the coefficient of eps^i mu1^j mu2^k,
with i <= LIM[0], j <= LIM[1], k <= LIM[2]. eps counts powers of the body's
(M, q) in the near-zone background; mu1, mu2 mark two first-order near-zone
fields whose bilinear product we want.

Spatial dependence: polynomials in x, y, z and powers of R, where R stands for
|x|; the derivative d() applies the chain rule dR/dx^i = x^i/R.
"""
import itertools
import sympy as sp

t, x, y, z = sp.symbols('t x y z', real=True)
X = (t, x, y, z)
R = sp.Symbol('R', positive=True)
eta = sp.diag(-1, 1, 1, 1)
LIM = [1, 1, 1]


def set_limits(e, m1, m2):
    LIM[:] = [e, m1, m2]


def ok(k):
    return all(k[i] <= LIM[i] for i in range(3))


def S(d=None):
    return {k: v for k, v in (d or {}).items() if v != 0 and ok(k)}


def add(*qs):
    out = {}
    for qq in qs:
        for k, v in qq.items():
            out[k] = out.get(k, 0) + v
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}


def scal(c, qq):
    return {k: sp.expand(c * v) for k, v in qq.items() if sp.expand(c * v) != 0}


def mul(p, qq):
    out = {}
    for k1, v1 in p.items():
        for k2, v2 in qq.items():
            k = (k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2])
            if ok(k):
                out[k] = out.get(k, 0) + v1 * v2
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}


def d_expr(f, a):
    if a == 0:
        return sp.diff(f, t)
    return sp.diff(f, X[a]) + sp.diff(f, R) * X[a] / R


def d(qq, a):
    out = {k: sp.expand(d_expr(v, a)) for k, v in qq.items()}
    return {k: v for k, v in out.items() if v != 0}


def const(c):
    return {(0, 0, 0): c} if c != 0 else {}


def mat(f):
    return [[f(a, b) for b in range(4)] for a in range(4)]


def matmul(A, B):
    return mat(lambda a, b: add(*[mul(A[a][c], B[c][b]) for c in range(4)]))


def inverse(g):
    """Inverse of g = eta + delta as a Neumann series, truncated."""
    delta = mat(lambda a, b: add(g[a][b], const(-eta[a, b])))
    et = mat(lambda a, b: const(eta[a, b]))
    ed = matmul(et, delta)
    term = et
    out = et
    for n in range(1, sum(LIM) + 1):
        term = matmul(ed, term)
        term = mat(lambda a, b: scal(-1, term[a][b]))
        out = mat(lambda a, b: add(out[a][b], term[a][b]))
    return out


def sqrt_minus_det(g, ginv_unused=None):
    """sqrt(-det g) for g = eta + delta, to second order in delta."""
    delta = mat(lambda a, b: add(g[a][b], const(-eta[a, b])))
    et = mat(lambda a, b: const(eta[a, b]))
    ed = matmul(et, delta)
    tr1 = add(*[ed[a][a] for a in range(4)])
    ed2 = matmul(ed, ed)
    tr2 = add(*[ed2[a][a] for a in range(4)])
    ed3 = matmul(ed2, ed)
    tr3 = add(*[ed3[a][a] for a in range(4)])
    # log det = tr1 - tr2/2 + tr3/3 ; sqrt = exp(L/2)
    L = add(tr1, scal(-sp.Rational(1, 2), tr2), scal(sp.Rational(1, 3), tr3))
    h = scal(sp.Rational(1, 2), L)
    h2 = mul(h, h)
    h3 = mul(h2, h)
    return add(const(1), h, scal(sp.Rational(1, 2), h2), scal(sp.Rational(1, 6), h3))


def christoffel(g, ginv):
    dg = [[[d(g[a][b], c) for c in range(4)] for b in range(4)] for a in range(4)]
    G = [[[None] * 4 for _ in range(4)] for _ in range(4)]
    low = [[[add(dg[dd][b][c], dg[dd][c][b], scal(-1, dg[b][c][dd])) for c in range(4)]
            for b in range(4)] for dd in range(4)]
    for a in range(4):
        for b in range(4):
            for c in range(b, 4):
                v = scal(sp.Rational(1, 2), add(*[mul(ginv[a][dd], low[dd][b][c]) for dd in range(4)]))
                G[a][b][c] = v
                G[a][c][b] = v
    return G


def ricci(G):
    Rm = mat(lambda a, b: {})
    for b in range(4):
        for c in range(b, 4):
            terms = []
            for a in range(4):
                terms.append(d(G[a][b][c], a))
                terms.append(scal(-1, d(G[a][b][a], c)))
                for dd in range(4):
                    terms.append(mul(G[a][a][dd], G[dd][b][c]))
                    terms.append(scal(-1, mul(G[a][c][dd], G[dd][b][a])))
            v = add(*terms)
            Rm[b][c] = v
            Rm[c][b] = v
    return Rm


def field(A):
    return mat(lambda a, b: add(d(A[b], a), scal(-1, d(A[a], b))))


def einstein_maxwell(g, A):
    """Residuals: E_ab = R_ab - 2 F_ac F_b^c + 1/2 g_ab F_cd F^cd ;
    Mx^a = d_b( sqrt(-g) g^ac g^bd F_cd )."""
    ginv = inverse(g)
    G = christoffel(g, ginv)
    Ric = ricci(G)
    F = field(A)
    Fmix = matmul(F, ginv)                 # F_a^b
    Fup = matmul(ginv, Fmix)               # F^ab = g^ac F_c^b ... (F_cd g^db)
    F2 = add(*[mul(F[a][b], Fup[a][b]) for a in range(4) for b in range(4)])
    E = mat(lambda a, b: {})
    for a in range(4):
        for b in range(a, 4):
            v = add(Ric[a][b],
                    scal(-2, add(*[mul(F[a][c], Fmix[b][c]) for c in range(4)])),
                    scal(sp.Rational(1, 2), mul(g[a][b], F2)))
            E[a][b] = v
            E[b][a] = v
    sg = sqrt_minus_det(g)
    dens = mat(lambda a, b: mul(sg, Fup[a][b]))
    Mx = [add(*[d(dens[a][b], b) for b in range(4)]) for a in range(4)]
    return E, Mx


# ---------------------------------------------------------------- l = 1 parts
nn = sp.symbols('n_x n_y n_z')


def _df(k):
    return sp.factorial2(k - 1) if k > 0 else 1


def sphere_avg(expr):
    poly = sp.Poly(sp.expand(expr), *nn)
    tot = 0
    for (p, q, s), c in poly.terms():
        if p % 2 or q % 2 or s % 2:
            continue
        tot += c * _df(p) * _df(q) * _df(s) / sp.factorial2(p + q + s + 1)
    return sp.expand(tot)


def radial_part(expr, rpow, tpow=0):
    """Coefficient of t^tpow R^rpow after x = R n (n not normalized: polynomial)."""
    e = sp.expand(expr).coeff(t, tpow) if tpow else sp.expand(expr).subs(t, 0)
    e = sp.expand(e.subs({x: R * nn[0], y: R * nn[1], z: R * nn[2]}, simultaneous=True))
    out = 0
    for term in sp.Add.make_args(e):
        c, p = term.as_coeff_exponent(R)
        if p == rpow:
            out += c
    return sp.expand(out)


def l1_scalar(s):
    return [3 * sphere_avg(s * nn[k]) for k in range(3)]


def l1_vector(v):
    """V_i = (P.n) n_i + C_i  ->  (P, C)."""
    avg = [sphere_avg(v[i]) for i in range(3)]
    vn = sum(v[i] * nn[i] for i in range(3))
    avgn = [sphere_avg(vn * nn[i]) for i in range(3)]
    # avg = P/3 + C ; avgn = (P + C)/3
    P = [sp.expand(sp.Rational(3, 2) * (3 * avgn[i] - avg[i])) for i in range(3)]
    C = [sp.expand(avg[i] - P[i] / 3) for i in range(3)]
    return P, C


def l1_tensor(T):
    """T_ij = D.n n_i n_j + E_(i n_j) + F.n delta_ij  ->  (D, E, F)."""
    tr = sum(T[i][i] for i in range(3))
    tnn = sum(T[i][j] * nn[i] * nn[j] for i in range(3) for j in range(3))
    p1 = [sphere_avg(tr * nn[k]) for k in range(3)]
    p2 = [sphere_avg(tnn * nn[k]) for k in range(3)]
    p3 = [sphere_avg(sum(T[k][j] * nn[j] for j in range(3))) for k in range(3)]
    Fv = [sp.expand(sp.Rational(3, 2) * (p1[k] - p2[k])) for k in range(3)]
    Ev = [sp.expand(3 * (p3[k] - p2[k])) for k in range(3)]
    Dv = [sp.expand(3 * p2[k] - Ev[k] - Fv[k]) for k in range(3)]
    return Dv, Ev, Fv


def l1_magnetic_vector(v):
    """V_i = (W x n)_i  ->  W  (projection on the magnetic l=1 subspace)."""
    # <eps_ijk n_j V_k> ... use <(n x V)_i> = <n x (W x n)> = W - <n (W.n)> = 2W/3
    nxv = [nn[(i + 1) % 3] * v[(i + 2) % 3] - nn[(i + 2) % 3] * v[(i + 1) % 3] for i in range(3)]
    return [sp.expand(sp.Rational(3, 2) * sphere_avg(c)) for c in nxv]
