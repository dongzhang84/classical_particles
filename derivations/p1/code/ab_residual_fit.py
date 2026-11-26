"""Fit the residual Delta(Q, R) (m2 = 1, e = 1) to a rational form.
Try:  Delta = P(Q, R) / ( Q^a R^b (R - Q^2)^c (R^2 - 2R + Q^2)^d )  with P polynomial."""
import itertools, mpmath as mp
from ab_residual_scan import residual
mp.mp.dps = 60
pts = []
for Q in [mp.mpf(x) / 10 for x in (11, 12, 13, 14, 15, 17, 19, 21, 23, 25)]:
    for R in [mp.mpf(x) for x in (9, 11, 13, 16, 19, 23, 27, 33)]:
        if R > Q**2 + 1:
            d = residual(Q, R - 1)[0]
            pts.append((Q, R, d))
print(len(pts), 'points', flush=True)
best = None
for a, b, c, dd in itertools.product(range(0, 3), range(0, 4), range(0, 5), range(0, 4)):
    for degQ, degR in [(3, 3), (4, 3), (4, 4)]:
        mons = [(i, j, w) for i in range(degQ + 1) for j in range(degR + 1) for w in (0, 1)]
        if len(mons) + 4 > len(pts):
            continue
        A = mp.matrix([[Q**i * R**j * mp.sqrt(R**2 - 2 * R + Q**2)**w for (i, j, w) in mons] for (Q, R, _) in pts])
        y = mp.matrix([dv * Q**a * R**b * (R - Q**2)**c * mp.sqrt(R**2 - 2 * R + Q**2)**dd for (Q, R, dv) in pts])
        sol, res = mp.qr_solve(A, y)
        rel = res / mp.norm(y)
        if best is None or rel < best[0]:
            best = (rel, a, b, c, dd, degQ, degR, sol, mons)
        if rel < mp.mpf('1e-30'):
            print('FIT', a, b, c, dd, degQ, degR, 'rel residual', mp.nstr(rel, 3))
            terms = [(mons[k], sol[k]) for k in range(len(mons)) if abs(sol[k]) > mp.mpf('1e-20')]
            for (i, j, w), v in terms:
                print(f'   Q^{i} R^{j} W^{w}: {mp.nstr(v, 20)}  ~ {mp.identify(v)}')
            raise SystemExit
print('best', mp.nstr(best[0], 3), best[1:7])
