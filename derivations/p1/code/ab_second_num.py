"""Second-order test with the AB exact solution, using the rules written in calc.tex
before this computation: gamma = rod midpoint; mass = area-charge mass mu; regular part of
the k-th order field = l=0 part of its r^0 coefficient and l=1 part of its r^1 coefficient."""
import sys
import mpmath as mp
from ab_first_num import fields, local_mu
mp.mp.dps = 60


def coeffs(e, Q, l, which, order):
    """Laurent coefficients (r^0 -> l=0 part, r^1 -> l=1 part) of F_k = (1/k!) d^k F/deps^k at 0."""
    def Fk(r, c):
        rho, z = r * mp.sqrt(1 - c**2), r * c
        return mp.diff(lambda ep: fields(ep, e, Q, l, rho, z)[which], 0, order) / mp.factorial(order)
    cs = [mp.cos(mp.pi * (k + mp.mpf(1) / 2) / 14) for k in range(14)]
    rs = [mp.mpf(10)**(-mp.mpf(k) / 4) * mp.mpf('1e-3') for k in range(9)]
    p = order                                   # most singular power r^-order
    npar = 8
    b0, b1, resid = [], [], 0
    for c in cs:
        vals = [Fk(r, c) * r**p for r in rs]
        A = mp.matrix([[r**j for j in range(npar)] for r in rs[:npar]])
        cf = mp.lu_solve(A, mp.matrix(vals[:npar]))
        # check with the extra radius
        r9 = rs[npar]
        resid = max(resid, abs(sum(cf[j] * r9**j for j in range(npar)) - vals[npar]))
        b0.append(cf[p]); b1.append(cf[p + 1])
    deg = 9
    Ac = mp.matrix([[c**j for j in range(deg + 1)] for c in cs[:deg + 1]])
    B0 = mp.lu_solve(Ac, mp.matrix(b0[:deg + 1])); B1 = mp.lu_solve(Ac, mp.matrix(b1[:deg + 1]))
    avg = sum(B0[j] * (1 + (-1)**j) / (j + 1) for j in range(deg + 1)) / 2
    p1 = mp.mpf(3) / 2 * sum(B1[j] * (1 + (-1)**(j + 1)) / (j + 2) for j in range(deg + 1))
    return avg, p1, resid


if __name__ == '__main__':
    e, Q, l = mp.mpf(1), mp.mpf(sys.argv[1]) if len(sys.argv) > 1 else mp.mpf('1.5'), \
        mp.mpf(sys.argv[2]) if len(sys.argv) > 2 else mp.mpf(10)
    H0, P0 = fields(0, e, Q, l, 0, 0)
    dH0 = mp.diff(lambda z: fields(0, e, Q, l, 0, z)[0], 0)
    dP0 = mp.diff(lambda z: fields(0, e, Q, l, 0, z)[1], 0)
    vH1, gH1, r1 = coeffs(e, Q, l, 0, 1); vP1, gP1, r2 = coeffs(e, Q, l, 1, 1)
    vH2, gH2, r3 = coeffs(e, Q, l, 0, 2); vP2, gP2, r4 = coeffs(e, Q, l, 1, 2)
    print('fit residuals:', [mp.nstr(x, 3) for x in (r1, r2, r3, r4)])
    print('first order  H:', mp.nstr(vH1, 12), mp.nstr(gH1, 12), ' Phi:', mp.nstr(vP1, 12), mp.nstr(gP1, 12))
    print('second order H:', mp.nstr(vH2, 12), mp.nstr(gH2, 12), ' Phi:', mp.nstr(vP2, 12), mp.nstr(gP2, 12))
    for eps in [mp.mpf('1e-4'), mp.mpf('2e-4'), mp.mpf('4e-4')]:
        mu = local_mu(eps, e, Q, l)
        def T(order):
            HR = H0 + eps * vH1 + (eps**2 * vH2 if order >= 2 else 0)
            dHR = dH0 + eps * gH1 + (eps**2 * gH2 if order >= 2 else 0)
            dPR = dP0 + eps * gP1 + (eps**2 * gP2 if order >= 2 else 0)
            return (eps * e / mu) / (-(dHR / (2 * mp.sqrt(HR))) / dPR) - 1
        T1, T2 = T(1), T(2)
        print(f'eps={mp.nstr(eps,3)}:  T(1st-order R)/eps^2 = {mp.nstr(T1/eps**2, 10)}   '
              f'T(2nd-order R)/eps^2 = {mp.nstr(T2/eps**2, 10)}   /eps^3 = {mp.nstr(T2/eps**3, 8)}')
