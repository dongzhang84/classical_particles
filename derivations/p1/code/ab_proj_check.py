import mpmath as mp
import ab_second_num as S
from ab_first_num import fields
mp.mp.dps = 60
e, Q, l = mp.mpf(1), mp.mpf('1.5'), mp.mpf(10)
base = S.coeffs(e, Q, l, 0, 2)
print('default:', [mp.nstr(x, 20) for x in base[:2]])
# exact Legendre projections via Gauss-Legendre quadrature in c, Laurent via fit in r
def coeffs_gl(which, order, ncq=30):
    def Fk(r, c):
        rho, z = r * mp.sqrt(1 - c**2), r * c
        return mp.diff(lambda ep: fields(ep, e, Q, l, rho, z)[which], 0, order) / mp.factorial(order)
    xs, ws = zip(*mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(5, mp.mp.prec))
    rs = [mp.mpf(10)**(-mp.mpf(k) / 4) * mp.mpf('1e-3') for k in range(9)]
    npar = 8; p = order
    s0 = s1 = 0
    for c, w in zip(xs, ws):
        vals = [Fk(r, c) * r**p for r in rs]
        A = mp.matrix([[r**j for j in range(npar)] for r in rs[:npar]])
        cf = mp.lu_solve(A, mp.matrix(vals[:npar]))
        s0 += w * cf[p] / 2
        s1 += w * cf[p + 1] * c * mp.mpf(3) / 2
    return s0, s1, len(xs)
gl = coeffs_gl(0, 2)
print('Gauss-Legendre (', gl[2], 'nodes):', [mp.nstr(x, 20) for x in gl[:2]])
