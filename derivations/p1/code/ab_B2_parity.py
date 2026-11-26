import re, sympy as sp
txt = open('/private/tmp/claude-501/-Users-dong-Projects-classical-particles/30be5932-e884-41b6-b8db-ace70328ed01/scratchpad/abc2.txt').read()
v, d = sp.symbols('v d')
B2 = sp.sympify(re.search(r'^B eps\^2: (.*)$', txt, re.M).group(1), locals={'v': v, 'd': d})
B1 = sp.sympify(re.search(r'^B eps\^1: (.*)$', txt, re.M).group(1), locals={'v': v, 'd': d})
f2 = sp.lambdify((v, d), B2, 'sympy'); f1 = sp.lambdify((v, d), B1, 'sympy')
def vd(u, t):
    p = (1 - u**2) / (2 * u); vv = 1 / u - 1; ww = p / t - 1
    return vv, ww - vv, (1 + u**2) / (2 * u), (p**2 / t - t) / 2 + 1
u = sp.Rational(1, 3); p = (1 - u**2) / (2 * u)
for t in [p / 5, p / 7]:
    a = vd(u, t); b = vd(u, -p**2 / t)
    print('R equal:', a[3] == b[3], ' B1 equal:', f1(a[0], a[1]) == f1(b[0], b[1]),
          ' B2(W) =', sp.N(f2(a[0], a[1]), 12), ' B2(-W) =', sp.N(f2(b[0], b[1]), 12))
