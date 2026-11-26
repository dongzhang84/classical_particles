"""Majumdar-Papapetrou limit check: big source Q -> 1+ (m2 = 1), small hole e/mu -> 1."""
import mpmath as mp
from ab_residual_scan import residual
from ab_numeric import coeffs as Bcoeffs
mp.mp.dps = 60
for q in ['1.1', '1.01', '1.001', '1.0001']:
    Q = mp.mpf(q)
    for R in [mp.mpf(11), mp.mpf(41)]:
        D, H0, mu1 = residual(Q, R - 1)
        B = Bcoeffs(mp.mpf(1), mp.mpf(1), Q, R - 1, n=3)
        print(f'Q={q} R={mp.nstr(R,3)}  Delta={mp.nstr(D, 10)}  Delta*R^3={mp.nstr(D*R**3, 8)}  '
              f'B1={mp.nstr(B[1], 8)}  B2={mp.nstr(B[2], 8)}  e/mu={mp.nstr(1/mu1, 8)}', flush=True)
