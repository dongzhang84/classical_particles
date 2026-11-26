"""Rewrite B1 (and check B2 numerically) from ab_closed2.py in physical variables:
units M2 = 1, e = 1; Q = e2, R = l + 1 (Bonnor radius), f = 1 - 2/R + Q^2/R^2."""
import re, sympy as sp
txt = open('/private/tmp/claude-501/-Users-dong-Projects-classical-particles/30be5932-e884-41b6-b8db-ace70328ed01/scratchpad/abc2.txt').read()
v, d = sp.symbols('v d', positive=True)
get = lambda key: sp.sympify(re.search(rf'^{re.escape(key)} (.*)$', txt, re.M).group(1), locals={'v': v, 'd': d})
B1, B2, mu1 = get('B eps^1:'), get('B eps^2:'), get('mu eps^1:')
Q, R = sp.symbols('Q R', positive=True)
L = R - 1
p = sp.sqrt(Q**2 - 1); W = sp.sqrt(L**2 + Q**2 - 1)
u = Q - p; t = W - L
vv = 1 / u - 1; ww = p / t - 1; dd = ww - vv
num = {Q: sp.Rational(3, 2), R: 11}
for name, ex in [('B1', B1), ('B2', B2), ('mu1', mu1)]:
    val = ex.subs({v: vv, d: dd}).subs(num)
    print(name, 'at Q=3/2, R=11:', sp.N(val, 15))
# closed form of B1 in physical variables
B1p = sp.simplify(sp.radsimp(B1.subs({v: vv, d: dd})))
mu1p = sp.simplify(sp.radsimp(mu1.subs({v: vv, d: dd})))
print('mu1 =', sp.factor(mu1p))
print('B1  =', sp.factor(B1p))
f = 1 - 2 / R + Q**2 / R**2
print('B1 * mu1 ... try B1/mu1 =', sp.factor(sp.simplify(B1p / mu1p)))
