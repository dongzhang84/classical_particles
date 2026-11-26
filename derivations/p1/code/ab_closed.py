"""Alekseev-Belinski equilibrium in invariant form (closed-form version).
On the small hole's horizon (c = 0 point):
    r+^2 = sigma1/kappa1 = D^2/(l^2 - sigma2^2 + gamma^2),  D = (m1 + sigma1)(m2 + l) - e1 (e2 + gamma),
    mu = (r+^2 + e1^2)/(2 r+),   N = kappa1/kappa_loc = 2 sigma1 r+/(r+^2 - e1^2).
Bonnor (test body) in RN(M_tot = m1+m2, Q_tot = e1+e2) at the radius with redshift N:
    e/mu = (M_tot R - Q_tot^2)/(Q_tot R N),   N^2 = 1 - 2 M_tot/R + Q_tot^2/R^2.
Instead of inverting N(R), compare at fixed invariants by computing
    B := Q_tot R N e / (mu (M_tot R - Q_tot^2))  with R defined from N.
"""
import sys, sympy as sp
eps = sp.Symbol('epsilon', positive=True)
e, M2, Q2, l = sp.symbols('e M_2 Q_2 ell', positive=True)
K = int(sys.argv[1]) if len(sys.argv) > 1 else 3
mm = sp.symbols(f'mm0:{K+2}')
m1 = eps * sum(mm[k] * eps**k for k in range(K + 2))
e1 = eps * e
gam = (M2 * e1 - m1 * Q2) / (l + m1 + M2)
S = lambda x, n: sp.expand(sp.series(x, eps, 0, n).removeO())
F = S((e1 - gam) * (Q2 + gam) - m1 * M2, K + 3)
sol = {}
for k in range(K + 1):
    sol[mm[k]] = sp.factor(sp.solve(sp.expand(F.coeff(eps, k + 1).subs(sol)), mm[k])[0])
    print(f'm1 eps^{k+1}:', sol[mm[k]], flush=True)
m1s = sp.expand(m1.subs(sol))
gs = (M2 * e1 - m1s * Q2) / (l + m1s + M2)
s1 = sp.sqrt(m1s**2 - e1**2 + 2 * e1 * gs)
s2sq = M2**2 - Q2**2 - 2 * Q2 * gs
D = (m1s + s1) * (M2 + l) - e1 * (Q2 + gs)
rp = D / sp.sqrt(l**2 - s2sq + gs**2)
mu = (rp**2 + e1**2) / (2 * rp)
N = 2 * s1 * rp / (rp**2 - e1**2)
Mt, Qt = m1s + M2, e1 + Q2
R = (Mt + sp.sqrt(Mt**2 - (1 - N**2) * Qt**2)) / (1 - N**2)
B = Qt * R * N * e1 / (mu * (Mt * R - Qt**2))
Bs = S(B, K)
print('mu leading:', sp.factor(S(mu / eps, 1)), flush=True)
for k in range(K):
    print(f'B eps^{k}:', sp.factor(Bs.coeff(eps, k)), flush=True)
