"""Static l = 1 response of a Reissner-Nordstrom black hole to a uniform external
field, exact in M and q (calc.tex Sec. 7, the pilot).

Background (Schwarzschild coordinates):  f = 1 - 2M/r + q^2/r^2,  A = -q/r dt.
Perturbation, static, m = 0:
  polar : h_tt = f H0 cos, h_rr = H2 cos / f, h_AB = r^2 K cos Omega_AB, h_r theta = -h1 sin,
          A_t = a cos, A_r = b cos, A_theta = -c sin
  axial : h_t phi = h0 sin^2,  A_phi = w sin^2
Residuals: E_ab = R_ab - 2 F_ac F_b^c + 1/2 g_ab F^2 ;  Mx^a = d_b(sqrt(-g) F^ab).
"""
import sympy as sp

t, r, th, ph = sp.symbols('t r theta phi')
X = (t, r, th, ph)
M, q, eps = sp.symbols('M q epsilon', positive=True)
f = 1 - 2 * M / r + q**2 / r**2
g0 = sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(th)**2)
A0 = [-q / r, 0, 0, 0]
H0, H2, K, h1, a, b, c, h0, w = [sp.Function(n)(r) for n in ('H0', 'H2', 'K', 'h1', 'a', 'b', 'c', 'h0', 'w')]
C, S = sp.cos(th), sp.sin(th)


def residuals(h, da):
    """Linearized residuals, using explicit first-order formulas (no series)."""
    gi = sp.simplify(g0.inv())
    d = lambda e, k: sp.diff(e, X[k])
    G0 = [[[sp.simplify(sum(gi[i, l] * (d(g0[l, j], k) + d(g0[l, k], j) - d(g0[j, k], l)) for l in range(4)) / 2)
            for k in range(4)] for j in range(4)] for i in range(4)]
    # delta Gamma^i_jk = 1/2 g^il (d_j h_lk + d_k h_lj - d_l h_jk) - g^il h_lm Gamma^m_jk
    dG = [[[sp.expand(sum(gi[i, l] * (d(h[l, k], j) + d(h[l, j], k) - d(h[j, k], l)) for l in range(4)) / 2
                      - sum(gi[i, l] * h[l, m] * G0[m][j][k] for l in range(4) for m in range(4)))
            for k in range(4)] for j in range(4)] for i in range(4)]
    dR = sp.zeros(4, 4)
    for j in range(4):
        for k in range(j, 4):
            e = sum(d(dG[i][j][k], i) - d(dG[i][j][i], k) for i in range(4))
            e += sum(dG[i][i][l] * G0[l][j][k] + G0[i][i][l] * dG[l][j][k]
                     - dG[i][k][l] * G0[l][j][i] - G0[i][k][l] * dG[l][j][i] for i in range(4) for l in range(4))
            dR[j, k] = dR[k, j] = e
    F0 = sp.Matrix(4, 4, lambda i, j: d(A0[j], i) - d(A0[i], j))
    df = sp.Matrix(4, 4, lambda i, j: d(da[j], i) - d(da[i], j))
    dgi = -gi * h * gi
    # 8 pi T = 2 F_ik F_jl g^kl - 1/2 g_ij F_kl F_mn g^km g^ln ; linearize
    def T8(F, Fb, ginv, g):
        return sp.Matrix(4, 4, lambda i, j: 2 * sum(F[i, k] * Fb[j, l] * ginv[k, l] for k in range(4) for l in range(4)))
    F2 = lambda F1, F2_, gA, gB: sum(F1[k, l] * F2_[m, n] * gA[k, m] * gB[l, n] for k in range(4) for l in range(4) for m in range(4) for n in range(4))
    s0 = F2(F0, F0, gi, gi)
    ds = 2 * F2(df, F0, gi, gi) + 2 * F2(F0, F0, dgi, gi)
    dT = sp.Matrix(4, 4, lambda i, j: 2 * sum((df[i, k] * F0[j, l] + F0[i, k] * df[j, l]) * gi[k, l]
                                            + F0[i, k] * F0[j, l] * dgi[k, l] for k in range(4) for l in range(4))
                   - h[i, j] * s0 / 2 - g0[i, j] * ds / 2)
    E = (dR - dT).applyfunc(lambda e: sp.simplify(sp.expand(e)))
    sq0 = r**2 * sp.sin(th)
    dsq = sq0 * sp.Rational(1, 2) * sum((gi * h)[i, i] for i in range(4))
    Fup0 = gi * F0 * gi
    dFup = gi * df * gi + dgi * F0 * gi + gi * F0 * dgi
    Mx = [sp.simplify(sp.expand(sum(d(sq0 * dFup[i, j] + dsq * Fup0[i, j], j) for j in range(4)))) for i in range(4)]
    return E, Mx


if __name__ == '__main__':
    import pickle
    hp = sp.zeros(4, 4)
    hp[0, 0] = f * H0 * C
    hp[1, 1] = H2 * C / f
    hp[2, 2] = r**2 * K * C
    hp[3, 3] = r**2 * S**2 * K * C
    hp[1, 2] = hp[2, 1] = -h1 * S
    Ep, Mp = residuals(hp, [a * C, b * C, -c * S, 0])
    ha = sp.zeros(4, 4)
    ha[0, 3] = ha[3, 0] = h0 * S**2
    Ea, Ma_ = residuals(ha, [0, 0, 0, w * S**2])
    pickle.dump({'Ep': Ep, 'Mp': Mp, 'Ea': Ea, 'Ma': Ma_}, open('l1_eqs.pkl', 'wb'))
    for name, (E, Mx) in {'polar': (Ep, Mp), 'axial': (Ea, Ma_)}.items():
        print('=====', name)
        for i in range(4):
            for j in range(i, 4):
                if E[i, j] != 0:
                    print(f'E{i}{j}:', sp.factor(E[i, j]))
        for i in range(4):
            if Mx[i] != 0:
                print(f'M{i}:', sp.factor(Mx[i]))
