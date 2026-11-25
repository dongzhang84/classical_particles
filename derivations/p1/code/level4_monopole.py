"""Level 4: first-order monopole corrections dM/r, dq/r in the near zone
contribute W.Q proportional to (dq E - dM a): the zeroth-order force with
corrected mass and charge."""
import sympy as sp
import trunc as T
from trunc import R
from level4_l1 import build, M, q, E
import level4_Q as LQ

dM, dq = sp.symbols('dM dq')
T.set_limits(1, 1, 1)
h2, A2 = LQ.field2('EB')
h1 = T.mat(lambda i, j: T.S({(1, 1, 0): 2 * dM / R if i == j else 0}))
A1 = [T.S({(1, 1, 0): -dq / R})] + [{} for _ in range(3)]
g, A = build(T.mat(lambda i, j: T.add(h1[i][j], h2[i][j])), [T.add(A1[i], A2[i]) for i in range(4)])
Ein, Mx = T.einstein_maxwell(g, A)
out = LQ.project(Ein, Mx, (1, 1, 1), -2)
print('W.Q =', [sp.factor(sp.Rational(2, 3) * out['D'][k] + out['E'][k]) for k in range(3)])
print('identity =', [sp.factor(-out['A'][k] + out['D'][k] + out['F'][k]) for k in range(3)])
