import mpmath as mp
import ab_residual_scan as S
from ab_first_num import fields, local_mu
from ab_second_num import coeffs
mp.mp.dps = 80
e, Q, l = mp.mpf(1), mp.mpf('1.5'), mp.mpf(10)
H0, P0 = fields(0, e, Q, l, 0, 0)
dH0 = mp.diff(lambda z: fields(0, e, Q, l, 0, z)[0], 0)
dP0 = mp.diff(lambda z: fields(0, e, Q, l, 0, z)[1], 0)
v = {(w, k): coeffs(e, Q, l, w, k)[:2] for w in (0, 1) for k in (1, 2)}
def Tov(eps):
    mu = local_mu(eps, e, Q, l)
    HR = H0 + eps * v[(0, 1)][0] + eps**2 * v[(0, 2)][0]
    dHR = dH0 + eps * v[(0, 1)][1] + eps**2 * v[(0, 2)][1]
    dPR = dP0 + eps * v[(1, 1)][1] + eps**2 * v[(1, 2)][1]
    return ((eps * e / mu) / (-(dHR / (2 * mp.sqrt(HR))) / dPR) - 1) / eps**2
for scale in (mp.mpf('1e-6'), mp.mpf('1e-8'), mp.mpf('1e-10')):
    epss = [k * scale for k in (1, 2, 3, 4, 5)]
    A = mp.matrix([[x**j for j in range(5)] for x in epss])
    cf = mp.lu_solve(A, mp.matrix([Tov(x) for x in epss]))
    print(mp.nstr(scale, 2), mp.nstr(cf[0], 30))
