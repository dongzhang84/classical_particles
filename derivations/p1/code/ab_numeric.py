"""Numerical check of the AB invariant balance: B(eps) = (e/mu)/(Bonnor value)."""
import mpmath as mp
mp.mp.dps = 60


def B_of(eps, e, M2, Q2, l):
    e1 = eps * e
    def gam(m1):
        return (M2 * e1 - m1 * Q2) / (l + m1 + M2)
    f = lambda m1: (e1 - gam(m1)) * (Q2 + gam(m1)) - m1 * M2
    m1 = mp.findroot(f, eps * Q2 * e * l / (M2**2 + M2 * l - Q2**2))
    g = gam(m1)
    s1 = mp.sqrt(m1**2 - e1**2 + 2 * e1 * g)
    s2sq = M2**2 - Q2**2 - 2 * Q2 * g
    D = (m1 + s1) * (M2 + l) - e1 * (Q2 + g)
    rp = D / mp.sqrt(l**2 - s2sq + g**2)
    mu = (rp**2 + e1**2) / (2 * rp)
    N = 2 * s1 * rp / (rp**2 - e1**2)
    Mt, Qt = m1 + M2, e1 + Q2
    R = (Mt + mp.sqrt(Mt**2 - (1 - N**2) * Qt**2)) / (1 - N**2)
    return Qt * R * N * e1 / (mu * (Mt * R - Qt**2)), mu / eps, s1


def coeffs(e, M2, Q2, l, n=4):
    hs = [mp.mpf(10)**(-k) for k in range(8, 8 + n + 2)]
    # polynomial fit in eps using small eps values
    xs = [mp.mpf(2)**(-k) * mp.mpf('1e-6') for k in range(n + 3)]
    ys = [B_of(x, e, M2, Q2, l)[0] for x in xs]
    A = mp.matrix([[x**j for j in range(n + 1)] for x in xs[:n + 1]])
    c = mp.lu_solve(A, mp.matrix(ys[:n + 1]))
    return [c[j] for j in range(n + 1)]


if __name__ == '__main__':
    # extremal check: M2 = Q2 = 1 is excluded (needs naked big source), so check
    # the test-limit coefficient and the eps^1, eps^2 coefficients for a few cases
    for (e, M2, Q2, l) in [(mp.mpf('0.5'), 1, mp.mpf('1.5'), 10), (mp.mpf('0.5'), 1, mp.mpf('1.5'), 20),
                           (mp.mpf('0.9'), 1, mp.mpf('1.2'), 10)]:
        B, mu1, s1 = B_of(mp.mpf('1e-8'), e, M2, Q2, l)
        print('e,M2,Q2,l =', e, M2, Q2, l, ' mu/eps =', mp.nstr(mu1, 12), ' sigma1 real:', mp.im(s1) == 0)
        c = coeffs(e, M2, Q2, l)
        print('   B = ', [mp.nstr(ci * mp.mpf('1e-6')**0, 12) for ci in c])
