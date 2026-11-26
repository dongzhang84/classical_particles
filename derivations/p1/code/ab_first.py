"""First-order AB fields: H1 = dH/deps, Phi1 = dPhi/deps at eps = 0, near the small hole.
Units M2 = 1. Small source at z = 0, big at z = l. Point rho = r s, z = r c (s^2 = 1 - c^2)."""
import sympy as sp
eps, r, c = sp.symbols('epsilon r c', real=True)
Q, l, e, a1 = sp.symbols('Q ell e a1', positive=True)
m1, e1 = eps * a1, eps * e
g = (e1 - m1 * Q) / (l + m1 + 1)
s1sq = m1**2 - e1**2 + 2 * e1 * g          # O(eps^2): drops out at first order
s2sq = 1 - Q**2 - 2 * Q * g
rho2, zz = r**2 * (1 - c**2), r * c
def bip(zk, sksq):
    w2 = rho2 + (zz - zk)**2
    x2 = (w2 + sksq + sp.sqrt((w2 + sksq)**2 - 4 * sksq * (zz - zk)**2)) / 2
    x = sp.sqrt(x2)
    return x, (zz - zk) / x
x1, c1 = bip(0, 0)                         # sigma1^2 = O(eps^2)
x2, c2 = bip(l, s2sq)
r1, r2 = x1 + m1, x2 + 1
D = r1 * r2 - (e1 - g - g * c2) * (Q + g - g * c1)
H = ((x1**2 + g**2 * (1 - c2**2)) * (x2**2 - s2sq + g**2 * (1 - c1**2))) / D**2
Phi = ((e1 - g) * (r2 - 1) + (Q + g) * (r1 - m1) + g * (m1 * c1 + c2)) / D
H0 = sp.simplify(H.subs(eps, 0)); P0 = sp.simplify(Phi.subs(eps, 0))
H1 = sp.diff(H, eps).subs(eps, 0); P1 = sp.diff(Phi, eps).subs(eps, 0)
import pickle
pickle.dump({'H0': H0, 'P0': P0, 'H1': H1, 'P1': P1}, open('ab_first.pkl', 'wb'))
for nm, F in [('H0', H0), ('Phi0', P0), ('H1', H1), ('Phi1', P1)]:
    ser = sp.series(F, r, 0, 2).removeO()
    for p in (-1, 0, 1):
        cf = sp.factor(sp.simplify(sp.expand(ser).coeff(r, p)))
        if cf != 0:
            print(f'{nm} r^{p}:', cf, flush=True)
