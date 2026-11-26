import mpmath as mp, sympy as sp
from ab_residual_scan import residual
mp.mp.dps = 60
Rs = [mp.mpf(x) for x in (200, 300, 400, 600, 800, 1200, 1600)]
def c3(Q):
    ys = [residual(Q, R - 1)[0] * R**3 for R in Rs]
    A = mp.matrix([[(1 / R)**j for j in range(len(Rs))] for R in Rs])
    return mp.lu_solve(A, mp.matrix(ys))[0]
Qs = [sp.Rational(k, 10) for k in (11, 12, 13, 15, 17, 20, 25, 30, 40)]
data = [(q, c3(mp.mpf(q.p) / q.q)) for q in Qs]
x = sp.Symbol('x')   # x = Q^2
for deg in (2, 3, 4):
    pts = data[:deg + 1]
    # fit c3*Q = polynomial in Q^2, rationalized
    vals = [sp.nsimplify(float(v * mp.mpf(q.p) / q.q), rational=True, tolerance=1e-9) for q, v in pts]
    poly = sp.interpolate([(q**2, vv) for (q, _), vv in zip(pts, vals)], x)
    err = max(abs(float(poly.subs(x, q**2)) - float(v * mp.mpf(q.p) / q.q)) for q, v in data)
    print('deg', deg, ' c3*Q =', sp.expand(poly), '  max err over all points:', err)
