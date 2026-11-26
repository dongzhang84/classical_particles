"""Large-R coefficients of Delta(Q, R) = c3/R^3 + c4/R^4 + ...  (m2 = 1, e = 1),
by polynomial extrapolation in 1/R."""
import sys, mpmath as mp
from ab_residual_scan import residual
mp.mp.dps = 60
Rs = [mp.mpf(x) for x in (200, 300, 400, 600, 800, 1200, 1600)]
for Q in [mp.mpf(x) for x in sys.argv[1:]] or [mp.mpf('1.2'), mp.mpf('1.5'), mp.mpf(2), mp.mpf(3)]:
    ys = [residual(Q, R - 1)[0] * R**3 for R in Rs]
    n = len(Rs)
    A = mp.matrix([[(1 / R)**j for j in range(n)] for R in Rs])
    cf = mp.lu_solve(A, mp.matrix(ys))
    print(f'Q={mp.nstr(Q,4)}  c3={mp.nstr(cf[0], 15)}  c4={mp.nstr(cf[1], 12)}  '
          f'c3*Q={mp.nstr(cf[0]*Q, 12)}  c3/Q={mp.nstr(cf[0]/Q, 12)}', flush=True)
