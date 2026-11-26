import mpmath as mp, sympy as sp
from ab_residual_scan import residual
mp.mp.dps = 70
def c3(Q, Rs):
    ys = [residual(Q, R - 1)[0] * R**3 for R in Rs]
    A = mp.matrix([[(1 / R)**j for j in range(len(Rs))] for R in Rs])
    return mp.lu_solve(A, mp.matrix(ys))[0]
R1 = [mp.mpf(x) for x in (1000, 1500, 2000, 3000, 4000, 6000, 8000, 12000)]
R2 = [mp.mpf(x) for x in (2000, 3000, 4000, 6000, 8000, 12000, 16000, 24000)]
for q in ['1.1', '1.2', '1.5', '2', '2.5', '3', '4']:
    Q = mp.mpf(q)
    a, b = c3(Q, R1), c3(Q, R2)
    print(f'Q={q}  c3(set1)={mp.nstr(a, 18)}  c3(set2)={mp.nstr(b, 18)}  c3*Q={mp.nstr(b*Q, 15)}', flush=True)
