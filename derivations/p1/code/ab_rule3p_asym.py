import mpmath as mp
from ab_rule3p import residual
mp.mp.dps = 60
Rs = [mp.mpf(x) for x in (400, 600, 800, 1200, 1600, 2400)]
for q in ['1.2', '1.5', '2']:
    Q = mp.mpf(q)
    for p in (3, 4):
        ys = [residual(Q, R - 1)[0] * R**p for R in Rs]
        A = mp.matrix([[(1 / R)**j for j in range(len(Rs))] for R in Rs])
        cf = mp.lu_solve(A, mp.matrix(ys))
        print(f'Q={q}  lim R^{p} Delta\' = {mp.nstr(cf[0], 14)}   /(Q^2-1)^2 = {mp.nstr(cf[0]/(Q**2-1)**2, 14)}', flush=True)
