"""C1 用到的三个恒等式的数值检验（任意度规、任意 A，在一个点上精确有理数计算）。

约定：signature (-+++)，G=c=1，Gauss 单位。
  F = dA，  4π T_ab = F_ac F_b^c - (1/4) g_ab F_cd F^cd，
  M^a ≡ ∇_b F^{ab}（源方程 M^a = 4π j^a）。

检验：
 (a) 4π ∇^a T_ab = - F_b^c M_c          ⇒  ∇^a(G_ab - 8π T_ab) = 2 F_b^c M_c
 (b) δ(∇_b F^{ab}) = ∇_b Y^{ab} - (1/2) h ∇_b F^{ab}，
     Y^{ab} = f^{ab} + (1/2) h F^{ab} - 2 h^{c[a} F_c^{b]}，h = g^{cd} h_cd
     （A → A + ε a，g → g + ε h，f = da；右边的 ∇、升降指标都用背景 g；背景不必满足任何方程。
      背景满足无源 Maxwell 方程时最后一项为零。）
 (c) δ(4π T_ab) = 2 F_(a^c f_b)c - h^{cd} F_ac F_bd - (1/4) h_ab F^2
                  - (1/2) g_ab ( f_cd F^cd - h^{ce} F_c^d F_ed )
"""
import random
import sympy as sp

random.seed(1)
X = sp.symbols('x0:4')
eps = sp.Symbol('eps')
P = {X[0]: sp.Rational(1, 3), X[1]: sp.Rational(-1, 2), X[2]: sp.Rational(1, 5), X[3]: sp.Rational(2, 7)}
eta = sp.diag(-1, 1, 1, 1)
N = range(4)


def rpoly(scale=sp.Rational(1, 10)):
    """随机二次多项式。"""
    c = lambda: sp.Rational(random.randint(-9, 9), 10)
    e = c()
    for i in N:
        e += c() * X[i]
        for j in N:
            if j >= i:
                e += c() * X[i] * X[j]
    return scale * e


def sym_field(scale):
    m = sp.zeros(4, 4)
    for i in N:
        for j in N:
            if j >= i:
                m[i, j] = m[j, i] = rpoly(scale)
    return m


g0 = eta + sym_field(sp.Rational(1, 10))
h = sym_field(1)
A = [rpoly(1) for _ in N]
a = [rpoly(1) for _ in N]


def jets(gexpr, Aexpr):
    """返回在点 P 上的 g, ∂g, F, ∂F（F 来自 A 的一、二阶导数）。"""
    g = gexpr.subs(P)
    dg = [[[sp.diff(gexpr[i, j], X[k]).subs(P) for k in N] for j in N] for i in N]
    F = [[(sp.diff(Aexpr[j], X[i]) - sp.diff(Aexpr[i], X[j])).subs(P) for j in N] for i in N]
    dF = [[[sp.diff(sp.diff(Aexpr[j], X[i]) - sp.diff(Aexpr[i], X[j]), X[k]).subs(P)
            for k in N] for j in N] for i in N]
    return g, dg, F, dF


def geom(g, dg):
    gi = g.inv()
    dgi = [[[-sum(gi[i, p] * dg[p][q][k] * gi[q, j] for p in N for q in N) for k in N] for j in N] for i in N]
    Gam = [[[sp.Rational(1, 2) * sum(gi[i, l] * (dg[l][j][k] + dg[l][k][j] - dg[j][k][l]) for l in N)
             for k in N] for j in N] for i in N]
    return gi, dgi, Gam


def divF_up(g, dg, F, dF):
    """M^a = ∇_b F^{ab} = (1/√-g) ∂_b(√-g F^{ab})。"""
    gi, dgi, Gam = geom(g, dg)
    Fu = [[sum(gi[a_, c] * gi[b, d] * F[c][d] for c in N for d in N) for b in N] for a_ in N]
    dFu = [[[sum(dgi[a_][c][k] * gi[b, d] * F[c][d] + gi[a_, c] * dgi[b][d][k] * F[c][d]
                 + gi[a_, c] * gi[b, d] * dF[c][d][k] for c in N for d in N)
             for k in N] for b in N] for a_ in N]
    trG = [sum(Gam[b][b][k] for b in N) for k in N]  # ∂_k ln√-g
    return [sum(dFu[a_][b][b] + trG[b] * Fu[a_][b] for b in N) for a_ in N], Fu, gi, dgi, Gam


# ---------- (a) ----------
g, dg, F, dF = jets(g0, A)
Mu, Fu, gi, dgi, Gam = divF_up(g, dg, F, dF)
Fsq = sum(F[c][d] * Fu[c][d] for c in N for d in N)
# 4πT_ab 及其偏导（用 jets 链式法则）
def T4pi(g_, gi_, F_):
    Fu_ = [[sum(gi_[a_, c] * gi_[b, d] * F_[c][d] for c in N for d in N) for b in N] for a_ in N]
    F2 = sum(F_[c][d] * Fu_[c][d] for c in N for d in N)
    return [[sum(F_[a_][c] * F_[b][d] * gi_[c, d] for c in N for d in N) - sp.Rational(1, 4) * g_[a_, b] * F2
             for b in N] for a_ in N]
t = sp.Symbol('t')
dT = [[[None] * 4 for _ in N] for _ in N]
for k in N:
    gk = g + t * sp.Matrix(4, 4, lambda i, j: dg[i][j][k])
    gik = gi + t * sp.Matrix(4, 4, lambda i, j: dgi[i][j][k])
    Fk = [[F[i][j] + t * dF[i][j][k] for j in N] for i in N]
    Tk = T4pi(gk, gik, Fk)
    for i in N:
        for j in N:
            dT[i][j][k] = sp.diff(Tk[i][j], t).subs(t, 0)
T0 = T4pi(g, gi, F)
divT = [sum(gi[a_, c] * (dT[a_][b][c] - sum(Gam[d][c][a_] * T0[d][b] + Gam[d][c][b] * T0[a_][d] for d in N))
            for a_ in N for c in N) for b in N]
Mdown = [sum(g[c, d] * Mu[d] for d in N) for c in N]
rhs = [-sum(F[b][c] * gi[c, d] * Mdown[d] for c in N for d in N) for b in N]
print("(a) 4π∇^aT_ab + F_b^c M_c =", [sp.simplify(divT[b] - rhs[b]) for b in N])

# ---------- (b) ----------
gE = g0 + eps * h
AE = [A[i] + eps * a[i] for i in N]
gE_, dgE, FE, dFE = jets(gE, AE)
MuE = divF_up(gE_, dgE, FE, dFE)[0]
lhs = [sp.diff(MuE[i], eps).subs(eps, 0) for i in N]

hj = h.subs(P); dh = [[[sp.diff(h[i, j], X[k]).subs(P) for k in N] for j in N] for i in N]
_, _, f, df = jets(g0, a)
# Y^{ab}，以及它的偏导
def Yup(g_, gi_, h_, F_, f_):
    Fu_ = [[sum(gi_[p, c] * gi_[q, d] * F_[c][d] for c in N for d in N) for q in N] for p in N]
    fu_ = [[sum(gi_[p, c] * gi_[q, d] * f_[c][d] for c in N for d in N) for q in N] for p in N]
    hu_ = gi_ * h_ * gi_
    trh = sum(gi_[p, q] * h_[p, q] for p in N for q in N)
    Fmix = [[sum(gi_[q, d] * F_[c][d] for d in N) for q in N] for c in N]  # F_c^q
    Y = [[fu_[p][q] + sp.Rational(1, 2) * trh * Fu_[p][q]
          - sum(hu_[c, p] * Fmix[c][q] - hu_[c, q] * Fmix[c][p] for c in N)
          for q in N] for p in N]
    return Y
Y0 = Yup(g, gi, hj, F, f)
dY = [[[None] * 4 for _ in N] for _ in N]
for k in N:
    gik = gi + t * sp.Matrix(4, 4, lambda i, j: dgi[i][j][k])
    gk = g + t * sp.Matrix(4, 4, lambda i, j: dg[i][j][k])
    hk = hj + t * sp.Matrix(4, 4, lambda i, j: dh[i][j][k])
    Fk = [[F[i][j] + t * dF[i][j][k] for j in N] for i in N]
    fk = [[f[i][j] + t * df[i][j][k] for j in N] for i in N]
    Yk = Yup(gk, gik, hk, Fk, fk)
    for i in N:
        for j in N:
            dY[i][j][k] = sp.diff(Yk[i][j], t).subs(t, 0)
trG = [sum(Gam[b][b][k] for b in N) for k in N]
trh0 = sum(gi[p, q] * hj[p, q] for p in N for q in N)
rhsb = [sum(dY[p][b][b] + trG[b] * Y0[p][b] for b in N) - sp.Rational(1, 2) * trh0 * Mu[p] for p in N]  # Y 反对称
print("(b) δ(∇_bF^ab) - ∇_bY^ab + ½hM^a =", [sp.simplify(lhs[i] - rhsb[i]) for i in N])

# ---------- (c) ----------
TE = T4pi(gE_, gE_.inv(), FE)
lhsc = [[sp.diff(TE[i][j], eps).subs(eps, 0) for j in N] for i in N]
hu = gi * hj * gi
Fm = [[sum(gi[q, d] * F[c][d] for d in N) for q in N] for c in N]   # F_c^q
fFsq = sum(f[c][d] * Fu[c][d] for c in N for d in N)
hFF = sum(hu[c, e] * F[c][d] * F[e][dd] * gi[d, dd] for c in N for e in N for d in N for dd in N)
rhsc = [[sum(F[i][c] * gi[c, d] * f[j][d] + F[j][c] * gi[c, d] * f[i][d] for c in N for d in N)
         - sum(hu[c, d] * F[i][c] * F[j][d] for c in N for d in N)
         - sp.Rational(1, 4) * hj[i, j] * Fsq
         - sp.Rational(1, 2) * g[i, j] * (fFsq - hFF)
         for j in N] for i in N]
print("(c) max|δ(4πT) - 公式| =", max(abs(sp.simplify(lhsc[i][j] - rhsc[i][j])) for i in N for j in N))
