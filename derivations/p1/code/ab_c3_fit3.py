import mpmath as mp, sympy as sp
from ab_residual_scan import residual
mp.mp.dps = 70
Rs = [mp.mpf(x) for x in (1000, 1500, 2000, 3000, 4000, 6000, 8000, 12000)]
def c3(Q):
    ys = [residual(Q, R - 1)[0] * R**3 for R in Rs]
    A = mp.matrix([[(1 / R)**j for j in range(len(Rs))] for R in Rs])
    return mp.lu_solve(A, mp.matrix(ys))[0]
qs = [sp.Rational(k, 10) for k in (11, 12, 13, 14, 15, 17, 20, 25, 30, 35, 40)]
vals = []
for q in qs:
    v = c3(mp.mpf(q.p) / q.q)
    from fractions import Fraction
    fr = Fraction(mp.nstr(v, 30)).limit_denominator(10**7)
    vals.append(sp.Rational(fr.numerator, fr.denominator))
print('exact c3 values:', list(zip(qs, vals)))
Qs = sp.Symbol('Q')
for k in range(0, 5):
    for deg in range(2, 9):
        pts = list(zip(qs, vals))[:deg + 1]
        P = sp.interpolate([(q, v * q**k) for q, v in pts], Qs)
        ok = all(sp.simplify(P.subs(Qs, q) - v * q**k) == 0 for q, v in zip(qs, vals))
        if ok:
            print(f'c3 = ({sp.factor(P)}) / Q^{k}')
            raise SystemExit
print('no polynomial/Q^k form found')
