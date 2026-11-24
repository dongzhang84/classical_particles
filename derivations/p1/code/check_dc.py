"""Independent check of the D^C column: brute-force Ricci without any truncation."""
import sympy as sp
from series import t, x, y, z, r, X, christoffel, ricci, DIRECTIONS, radial_coeff
M, e = sp.symbols('M epsilon')
C = sp.symbols('C_x C_y C_z')
psi = 2*M/r
g = sp.diag(-1+psi, 1+psi, 1+psi, 1+psi)
for i in range(3):
    g[0, i+1] = g[i+1, 0] = e*t*C[i]
ginv = g.inv()
G = christoffel(g, ginv)
Rfull = ricci(G)
R00 = sp.diff(Rfull[0, 0], e).subs(e, 0)
R00 = sp.series(R00, M, 0, 2).removeO()
for n in DIRECTIONS[:3]:
    print(n, sp.factor(radial_coeff(R00, n, 0, -2)), ' expected ours: -M C.n =', sp.factor(-M*sum(c*ni for c, ni in zip(C, n))))
