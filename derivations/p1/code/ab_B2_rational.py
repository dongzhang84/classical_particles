"""Reconstruct B2(Q, R) as a rational function (units M2 = 1, e = 1) from exact values at
points where all radicals are rational (Q = (1+u^2)/2u, R - 1 = (p^2/t - t)/2)."""
import re, itertools, sympy as sp
txt = open('/private/tmp/claude-501/-Users-dong-Projects-classical-particles/30be5932-e884-41b6-b8db-ace70328ed01/scratchpad/abc2.txt').read()
v, d = sp.symbols('v d', positive=True)
B2 = sp.sympify(re.search(r'^B eps\^2: (.*)$', txt, re.M).group(1), locals={'v': v, 'd': d})
B2f = sp.lambdify((v, d), B2, 'sympy')
Qs, Rs = sp.symbols('Q R')


def point(u, t):
    Q = (1 + u**2) / (2 * u); p = (1 - u**2) / (2 * u)
    L = (p**2 / t - t) / 2
    vv = 1 / u - 1; ww = p / t - 1
    return Q, L + 1, vv, ww - vv


pts = []
for u in [sp.Rational(1, k) for k in (2, 3, 4, 5, 6, 7)] + [sp.Rational(2, 5), sp.Rational(3, 7), sp.Rational(2, 7), sp.Rational(3, 8)]:
    p = (1 - u**2) / (2 * u)
    for frac in [sp.Rational(1, k) for k in range(2, 14)]:
        t = p * frac
        Q, R, vv, dd = point(u, t)
        if dd > 0 and R > Q**2:
            pts.append((Q, R, B2f(vv, dd)))
print(len(pts), 'points')
# ansatz: B2 = N(Q,R) / (Q^2 (Q^2 - R)^6), N polynomial
for degQ, degR in [(8, 6), (10, 7), (12, 8), (14, 9)]:
    mons = [(i, j) for i in range(degQ + 1) for j in range(degR + 1)]
    if len(mons) > len(pts) - 5:
        print('not enough points for', degQ, degR); continue
    A = sp.Matrix([[Q**i * R**j for (i, j) in mons] for (Q, R, _) in pts])
    b = sp.Matrix([val * Q**2 * (Q**2 - R)**6 for (Q, R, val) in pts])
    try:
        sol = A.solve_least_squares(b)
    except Exception as ex:
        print('fail', ex); continue
    res = (A * sol - b).norm()
    print(degQ, degR, 'residual', sp.N(res))
    if res == 0:
        Np = sum(sol[k] * Qs**i * Rs**j for k, (i, j) in enumerate(mons))
        print('B2 =', sp.factor(Np / (Qs**2 * (Qs**2 - Rs)**6)))
        break
