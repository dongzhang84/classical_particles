import mpmath as mp
from ab_residual_scan import residual
mp.mp.dps = 60
for Q in [mp.mpf('1.2'), mp.mpf('1.5'), mp.mpf(2), mp.mpf(3)]:
    row = []
    for R in [mp.mpf(50), mp.mpf(100), mp.mpf(200), mp.mpf(400)]:
        d = residual(Q, R - 1)[0]
        row.append(d * R**3 / Q)
    print('Q =', mp.nstr(Q, 3), ' Delta R^3 / Q at R = 50,100,200,400:', [mp.nstr(x, 10) for x in row], flush=True)
