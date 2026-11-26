"""AB closure test with rule 3': regular parts taken in the MP-linear variables
psi1 = 1/sqrt(H), psi2 = Phi/sqrt(H); then sqrt(H^R) = 1/psi1^R, Phi^R = psi2^R/psi1^R.
pred = -d sqrt(H^R) / d Phi^R = d psi1 / (psi1 d psi2 - psi2 d psi1)   (at gamma)
T = (e/mu)/pred - 1 ;  Delta' = lim T/eps^2."""
import sys, mpmath as mp
from ab_first_num import fields, local_mu
mp.mp.dps = 60


def psis(eps, e, Q, l, rho, z):
    H, P = fields(eps, e, Q, l, rho, z)
    s = mp.sqrt(H)
    return 1 / s, P / s


def coeffs(e, Q, l, which, order):
    def Fk(r, c):
        rho, z = r * mp.sqrt(1 - c**2), r * c
        return mp.diff(lambda ep: psis(ep, e, Q, l, rho, z)[which], 0, order) / mp.factorial(order)
    xs, ws = zip(*mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(4, mp.mp.prec))
    rs = [mp.mpf(10)**(-mp.mpf(k) / 4) * mp.mpf('1e-3') for k in range(9)]
    npar, p = 8, order
    s0 = s1 = 0
    for c, w in zip(xs, ws):
        vals = [Fk(r, c) * r**p for r in rs]
        A = mp.matrix([[r**j for j in range(npar)] for r in rs[:npar]])
        cf = mp.lu_solve(A, mp.matrix(vals[:npar]))
        s0 += w * cf[p] / 2
        s1 += w * cf[p + 1] * c * mp.mpf(3) / 2
    return s0, s1


def residual(Q, l, e=mp.mpf(1)):
    v0 = [psis(0, e, Q, l, 0, 0)[k] for k in (0, 1)]
    g0 = [mp.diff(lambda z: psis(0, e, Q, l, 0, z)[k], 0) for k in (0, 1)]
    c = {(k, o): coeffs(e, Q, l, k, o) for k in (0, 1) for o in (1, 2)}
    epss = [mp.mpf(j) * mp.mpf('1e-6') for j in (1, 2, 3, 4, 5)]
    out = []
    for eps in epss:
        mu = local_mu(eps, e, Q, l)
        v = [v0[k] + eps * c[(k, 1)][0] + eps**2 * c[(k, 2)][0] for k in (0, 1)]
        g = [g0[k] + eps * c[(k, 1)][1] + eps**2 * c[(k, 2)][1] for k in (0, 1)]
        pred = g[0] / (v[0] * g[1] - v[1] * g[0])
        out.append(((eps * e / mu) / pred - 1) / eps**2)
    A = mp.matrix([[x**j for j in range(len(epss))] for x in epss])
    return mp.lu_solve(A, mp.matrix(out))[0], out[0]


if __name__ == '__main__':
    for q, R in [('1.0001', 11), ('1.001', 11), ('1.01', 11), ('1.2', 11), ('1.5', 11), ('2', 11), ('1.5', 21), ('1.5', 41)]:
        D, T1 = residual(mp.mpf(q), mp.mpf(R) - 1)
        print(f'Q={q} R={R}:  Delta\' = lim T/eps^2 = {mp.nstr(D, 12)}   (T/eps^2 at eps=1e-6: {mp.nstr(T1, 10)})', flush=True)
