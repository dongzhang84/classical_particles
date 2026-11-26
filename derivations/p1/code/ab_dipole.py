"""Is the horizon-rod midpoint the mass centre?  Laurent coefficients of the AB scalars about
the rod midpoint: k-th order fields F_k (k = 1, 2, 3), the r^{-2} l = 1 parts (dipoles)."""
import mpmath as mp
from ab_first_num import fields
mp.mp.dps = 60
e, Q, l = mp.mpf(1), mp.mpf('1.5'), mp.mpf(10)


def laurent(which, order, powers):
    def Fk(r, c):
        rho, z = r * mp.sqrt(1 - c**2), r * c
        return mp.diff(lambda ep: fields(ep, e, Q, l, rho, z)[which], 0, order) / mp.factorial(order)
    xs, ws = zip(*mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(4, mp.mp.prec))
    rs = [mp.mpf(10)**(-mp.mpf(k) / 4) * mp.mpf('1e-3') for k in range(10)]
    npar, p = 9, order + 1
    res = {pw: [0, 0] for pw in powers}
    for c, w in zip(xs, ws):
        vals = [Fk(r, c) * r**p for r in rs]
        A = mp.matrix([[r**j for j in range(npar)] for r in rs[:npar]])
        cf = mp.lu_solve(A, mp.matrix(vals[:npar]))
        for pw in powers:
            co = cf[pw + p]
            res[pw][0] += w * co / 2                       # l = 0
            res[pw][1] += w * co * c * mp.mpf(3) / 2       # l = 1
    return res


for which, nm in [(0, 'H'), (1, 'Phi')]:
    for k in (1, 2, 3):
        res = laurent(which, k, [-3, -2, -1, 0])
        print(nm, f'eps^{k}:', {pw: (mp.nstr(v[0], 12), mp.nstr(v[1], 12)) for pw, v in res.items()}, flush=True)
