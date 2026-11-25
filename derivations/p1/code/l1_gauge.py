"""Gauge transformations of the static l = 1 polar fields (H0, H2, K, h1, a)
on RN, and gauge-invariant combinations."""
import sympy as sp
from l1_response import t, r, th, ph, X, M, q, f, g0, A0
Xf, Yf = sp.Function('X')(r), sp.Function('Y')(r)
C, S = sp.cos(th), sp.sin(th)
xi = [0, Xf * C, -Yf * S / r**2, 0]          # xi^mu (contravariant)
d = lambda e, k: sp.diff(e, X[k])
Lg = sp.Matrix(4, 4, lambda m, n: sp.simplify(sum(xi[k] * d(g0[m, n], k) + g0[k, n] * d(xi[k], m) + g0[m, k] * d(xi[k], n) for k in range(4))))
LA = [sp.simplify(sum(xi[k] * d(A0[n], k) + A0[k] * d(xi[k], n) for k in range(4))) for n in range(4)]
dH0 = sp.simplify(Lg[0, 0] / (f * C))
dH2 = sp.simplify(Lg[1, 1] * f / C)
dK = sp.simplify(Lg[2, 2] / (r**2 * C))
dK3 = sp.simplify(Lg[3, 3] / (r**2 * S**2 * C))
dh1 = sp.simplify(-Lg[1, 2] / S)
da = sp.simplify(LA[0] / C)
print('dH0 =', dH0); print('dH2 =', dH2); print('dK =', dK, ' (33 check:', sp.simplify(dK - dK3), ')')
print('dh1 =', dh1); print('da =', da)
print('other comps:', [sp.simplify(Lg[i, j]) for i, j in [(0, 1), (0, 2), (0, 3), (1, 3), (2, 3)]], LA[1:])

# invariants: X_eff = r^2 a / q, Y_eff = r X_eff - r^2 K / 2
import pickle
from l1_response import H0, H2, K, h1, a
ans, s, free = pickle.load(open('l1_polar_sol.pkl', 'rb'))
sol = {fld: sp.expand(ans[fld].subs(s)) for fld in (H0, H2, K, h1, a)}
Xe = r**2 * sol[a] / q
Ye = r * Xe - r**2 * sol[K] / 2
fr = r**2 - 2 * M * r + q**2
I1 = sp.factor(sol[H0] - 2 * (q**2 - M * r) * Xe / (r * fr))
I2 = sp.factor(sol[H2] - 2 * sp.diff(Xe, r) - 2 * (q**2 - M * r) * Xe / (r * fr))
I3 = sp.factor(sol[h1] - (r**3 * Xe + (r * sp.diff(Ye, r) - 2 * Ye) * fr) / (r * fr))
print('I1 =', I1); print('I2 =', I2); print('I3 =', I3)
print('surviving free parameters:', sorted((I1.free_symbols | I2.free_symbols | I3.free_symbols) - {r, M, q}, key=str))
