"""Effective-mass check, efficient version. In the gauge of mass2_gauge.py the
first-order metric is diagonal (h1 = 0), so the full metric to eps^2 is
  ds^2 = -A dt^2 + B dr^2 + C (dth^2 + sin^2 th dph^2),  A_t = Phi,
with A, B, C, Phi functions of (r, th). Ricci/Maxwell are computed once for
generic A, B, C, Phi; then the eps-expansion is substituted and the eps^2
coefficient projected on l = 0 (exact 3-point Gauss-Legendre in cos th).
Usage: python mass2_diag.py q   (M = 1, E = 1)"""
import pickle, sys, sympy as sp
from l1_response import r, th, M as Ms, q as qs, H0, H2, K, h1, a
t, ph, eps = sp.symbols('t phi epsilon')
X = [t, r, th, ph]
A, B, Cf, Phi = [sp.Function(n)(r, th) for n in ('A', 'B', 'Cf', 'Phi')]
g = sp.diag(-A, B, Cf, Cf*sp.sin(th)**2)
gi = sp.diag(-1/A, 1/B, 1/Cf, 1/(Cf*sp.sin(th)**2))
d = lambda e, k: sp.diff(e, X[k])
Gam = [[[sp.simplify(sum(gi[i, l]*(d(g[l, j], k) + d(g[l, k], j) - d(g[j, k], l)) for l in range(4))/2)
         for k in range(4)] for j in range(4)] for i in range(4)]
def Ric(j, k):
    e = sum(d(Gam[i][j][k], i) - d(Gam[i][j][i], k) for i in range(4))
    e += sum(Gam[i][i][l]*Gam[l][j][k] - Gam[i][k][l]*Gam[l][j][i] for i in range(4) for l in range(4))
    return sp.simplify(e)
Av = [Phi, 0, 0, 0]
F = sp.Matrix(4, 4, lambda i, j: d(Av[j], i) - d(Av[i], j))
Fup = gi*F*gi
F2 = sum(F[i, j]*Fup[i, j] for i in range(4) for j in range(4))
T8 = lambda i, j: 2*sum(F[i, k]*gi[k, l]*F[j, l] for k in range(4) for l in range(4)) - g[i, j]*F2/2
sq = sp.sqrt(A*B)*Cf*sp.sin(th)
eqs = {(0, 0): Ric(0, 0) - T8(0, 0), (1, 1): Ric(1, 1) - T8(1, 1), (2, 2): (Ric(2, 2) - T8(2, 2) + (Ric(3, 3) - T8(3, 3))/sp.sin(th)**2)/(2*r**2),
       'M': sum(d(sq*Fup[0, m], m) for m in range(4))/(r**2*sp.sin(th))}
print('generic done', flush=True)

qv = sp.Rational(sys.argv[1]); Mv = 1
G, rest = pickle.load(open('mass2_gauge.pkl', 'rb'))
Esym = sp.Symbol('E')
G = {k: sp.factor(v.subs({u: 0 for u in rest}).subs({Ms: Mv, qs: qv, Esym: 1})) for k, v in G.items()}
assert G[h1] == 0
f = 1 - 2*Mv/r + qv**2/r**2
P, R, b = [sp.Function(n)(r) for n in ('P', 'R', 'b')]
c = sp.cos(th)
sub = {A: f*(1 - eps*G[H0]*c - eps**2*P), B: (1 + eps*G[H2]*c + eps**2*R)/f,
       Cf: r**2*(1 + eps*G[K]*c), Phi: -qv/r + eps*G[a]*c + eps**2*b}
NODES = [(sp.Integer(0), sp.Rational(8, 9)), (sp.sqrt(sp.Rational(3, 5)), sp.Rational(5, 9)), (-sp.sqrt(sp.Rational(3, 5)), sp.Rational(5, 9))]
out = {}
for key, e in eqs.items():
    e2 = sp.diff(e.subs(sub).doit(), eps, 2).subs(eps, 0)/2
    e0 = e.subs(sub).doit().subs(eps, 0)
    tot = 0
    for cv, w in NODES:
        tot += w*e2.subs(th, sp.acos(cv))
    out[key] = sp.simplify(tot/2)
    print(key, 'done; zeroth-order residual =', sp.simplify(e0), flush=True)
pickle.dump(out, open(f'mass2_diag_q{str(qv).replace("/", "_")}.pkl', 'wb'))
for k, v in out.items():
    print(k, ':', sp.collect(sp.expand(sp.numer(sp.together(v))), [P, R, b], evaluate=True))
print('ALL DONE', flush=True)
