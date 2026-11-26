"""First-order test with the AB exact solution (numerical, high precision).
Claim to test: the exact equilibrium satisfies the test-body static balance in the regular
fields to first order:   e/mu = d_z sqrt(H^R) / d_z Phi^R   + O(eps^2).
H, Phi: AB scalars (eqs. (3),(5),(6),(8)); units m2 = 1; small hole at z = 0.
Regular part at first order: value = l=0 part of the r^0 coefficient of H1, gradient = l=1 part
of the r^1 coefficient (first-order singular field has no smooth l=1 r^1 piece)."""
import mpmath as mp
mp.mp.dps = 50


def solution(eps, e, Q, l):
    e1 = eps * e
    gam = lambda m1: (e1 - m1 * Q) / (l + m1 + 1)
    f = lambda m1: (e1 - gam(m1)) * (Q + gam(m1)) - m1
    m1 = mp.findroot(f, eps * Q * e * l / (1 + l - Q**2)) if eps != 0 else mp.mpf(0)
    return m1, e1, gam(m1)


def fields(eps, e, Q, l, rho, z):
    if eps == 0:                      # source 1 absent: RN of source 2 (common factor x1^2 cancelled)
        s2sq = 1 - Q**2
        x2 = mp.sqrt((rho**2 + (z - l)**2 + s2sq + mp.sqrt((rho**2 + (z - l)**2 + s2sq)**2
                                                            - 4 * s2sq * (z - l)**2)) / 2)
        r2 = x2 + 1
        return (x2**2 - s2sq) / r2**2, Q / r2
    m1, e1, g = solution(eps, e, Q, l)
    s1sq = m1**2 - e1**2 + 2 * e1 * g
    s2sq = 1 - Q**2 - 2 * Q * g
    def bip(zk, sksq):
        w2 = rho**2 + (z - zk)**2
        x2 = (w2 + sksq + mp.sqrt((w2 + sksq)**2 - 4 * sksq * (z - zk)**2)) / 2
        x = mp.sqrt(x2)
        return x, ((z - zk) / x if x != 0 else mp.mpf(0))
    x1, c1 = bip(0, s1sq); x2, c2 = bip(l, s2sq)
    r1, r2 = x1 + m1, x2 + 1
    D = r1 * r2 - (e1 - g - g * c2) * (Q + g - g * c1)
    H = ((x1**2 - s1sq + g**2 * (1 - c2**2)) * (x2**2 - s2sq + g**2 * (1 - c1**2))) / D**2
    Phi = ((e1 - g) * (r2 - 1) + (Q + g) * (r1 - m1) + g * (m1 * c1 + c2)) / D
    return H, Phi


def local_mu(eps, e, Q, l):
    m1, e1, g = solution(eps, e, Q, l)
    s1 = mp.sqrt(m1**2 - e1**2 + 2 * e1 * g); s2sq = 1 - Q**2 - 2 * Q * g
    D = (m1 + s1) * (1 + l) - e1 * (Q + g)
    rp = D / mp.sqrt(l**2 - s2sq + g**2)
    return (rp**2 + e1**2) / (2 * rp)


def first_order_coeffs(e, Q, l, which):
    """r-expansion coefficients of F1 = dF/deps|0 near the small hole, projected on
    Legendre P0 (r^0) and P1 (r^1)."""
    def F1(r, c):
        rho, z = r * mp.sqrt(1 - c**2), r * c
        return mp.diff(lambda ep: fields(ep, e, Q, l, rho, z)[which], 0)
    # Gauss-Legendre in c, polynomial fit in r on small radii
    nodes = [(-mp.sqrt(mp.mpf(3) / 5), mp.mpf(5) / 9), (0, mp.mpf(8) / 9), (mp.sqrt(mp.mpf(3) / 5), mp.mpf(5) / 9)]
    xs, ws = mp.gauss_legendre if False else (None, None)
    cs = [mp.cos(mp.pi * (k + mp.mpf(1) / 2) / 12) for k in range(12)]   # Chebyshev nodes in c
    rs = [mp.mpf(10)**(-k) for k in (3, 3.3, 3.6, 3.9, 4.2, 4.5)]
    out0, out1 = [], []
    for c in cs:
        vals = [F1(r, c) * r for r in rs]     # r F1 = A + B r + C r^2 + D r^3 ...
        A = mp.matrix([[r**j for j in range(5)] for r in rs[:5]])
        coef = mp.lu_solve(A, mp.matrix(vals[:5]))
        out0.append(coef[1]); out1.append(coef[2])
    # Legendre projections by fitting polynomials in c (B(c) and C(c) are polynomials of low degree)
    deg = 6
    Ac = mp.matrix([[c**j for j in range(deg + 1)] for c in cs[:deg + 1]])
    b0 = mp.lu_solve(Ac, mp.matrix(out0[:deg + 1])); b1 = mp.lu_solve(Ac, mp.matrix(out1[:deg + 1]))
    # P0 part of B: integral average; P1 part of C: (3/2) int c C dc
    avg = lambda b: sum(b[j] * (1 + (-1)**j) / (j + 1) for j in range(deg + 1)) / 2
    p1 = lambda b: mp.mpf(3) / 2 * sum(b[j] * (1 + (-1)**(j + 1)) / (j + 2) for j in range(deg + 1))
    return avg(b0), p1(b1)


if __name__ == '__main__':
    e, Q, l = mp.mpf(1), mp.mpf('1.5'), mp.mpf(10)
    # background values and gradients at the hole (eps = 0)
    h = mp.mpf('1e-12')
    H0 = fields(0, e, Q, l, 0, 0)[0]; P0 = fields(0, e, Q, l, 0, 0)[1]
    dH0 = mp.diff(lambda z: fields(0, e, Q, l, 0, z)[0], 0)
    dP0 = mp.diff(lambda z: fields(0, e, Q, l, 0, z)[1], 0)
    vH, gH = first_order_coeffs(e, Q, l, 0)
    vP, gP = first_order_coeffs(e, Q, l, 1)
    print('H0, dH0, dP0 =', mp.nstr(H0, 12), mp.nstr(dH0, 12), mp.nstr(dP0, 12))
    print('H1 regular value, gradient =', mp.nstr(vH, 12), mp.nstr(gH, 12))
    print('Phi1 regular value, gradient =', mp.nstr(vP, 12), mp.nstr(gP, 12))
    for eps in [mp.mpf('1e-4'), mp.mpf('2e-4'), mp.mpf('4e-4')]:
        mu = local_mu(eps, e, Q, l)
        HR = H0 + eps * vH; dHR = dH0 + eps * gH; dPR = dP0 + eps * gP
        pred = -(dHR / (2 * mp.sqrt(HR))) / dPR
        T = (eps * e / mu) / pred - 1
        T0 = (eps * e / mu) / (-(dH0 / (2 * mp.sqrt(H0))) / dP0) - 1
        print(f'eps={mp.nstr(eps,3)}:  T_R = {mp.nstr(T, 12)}  T_R/eps^2 = {mp.nstr(T/eps**2, 10)}   T_background/eps = {mp.nstr(T0/eps, 10)}')
