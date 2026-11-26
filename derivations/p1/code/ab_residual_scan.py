"""Scan the second-order residual  Delta = lim T(2nd-order R)/eps^2  over (Q, l)  (m2 = 1, e = 1)."""
import mpmath as mp
from ab_first_num import fields, local_mu
from ab_second_num import coeffs
mp.mp.dps = 60


def residual(Q, l, e=mp.mpf(1)):
    H0, P0 = fields(0, e, Q, l, 0, 0)
    dH0 = mp.diff(lambda z: fields(0, e, Q, l, 0, z)[0], 0)
    dP0 = mp.diff(lambda z: fields(0, e, Q, l, 0, z)[1], 0)
    v = {}
    for which in (0, 1):
        for k in (1, 2):
            v[(which, k)] = coeffs(e, Q, l, which, k)[:2]
    out = []
    epss = [mp.mpf(k) * mp.mpf('1e-6') for k in (1, 2, 3, 4, 5)]
    for eps in epss:
        mu = local_mu(eps, e, Q, l)
        HR = H0 + eps * v[(0, 1)][0] + eps**2 * v[(0, 2)][0]
        dHR = dH0 + eps * v[(0, 1)][1] + eps**2 * v[(0, 2)][1]
        dPR = dP0 + eps * v[(1, 1)][1] + eps**2 * v[(1, 2)][1]
        out.append(((eps * e / mu) / (-(dHR / (2 * mp.sqrt(HR))) / dPR) - 1) / eps**2)
    # polynomial extrapolation to eps -> 0
    A = mp.matrix([[x**j for j in range(len(epss))] for x in epss])
    cf = mp.lu_solve(A, mp.matrix(out))
    return cf[0], H0, mu / epss[-1]


if __name__ == '__main__':
    for Q in [mp.mpf('1.2'), mp.mpf('1.5'), mp.mpf(2)]:
        for l in [mp.mpf(10), mp.mpf(20), mp.mpf(40)]:
            if l + 1 <= Q**2:
                continue
            d, H0, mu1 = residual(Q, l)
            R = l + 1
            print(f'Q={mp.nstr(Q,3)} R={mp.nstr(R,4)}  Delta={mp.nstr(d, 12)}  Delta*R^2={mp.nstr(d*R**2, 10)}  '
                  f'Delta*R^3={mp.nstr(d*R**3, 10)}  mu/eps={mp.nstr(mu1, 8)}', flush=True)
