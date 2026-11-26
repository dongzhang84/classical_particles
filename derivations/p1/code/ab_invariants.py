"""Alekseev-Belinski equilibrium in invariant form, expanded in the small black hole.

Exact solution (Alekseev & Belinski, PRD 76, 021501 (2007), eqs. (1)-(10)):
source 1 = small black hole at z1, source 2 = big source at z2 = z1 + l.
On the horizon rod of source 1 (rho -> 0, z = z1 + sigma1 c, theta2 = pi,
r2 - m2 = l - sigma1 c), from kappa^2 = lim H/(rho^2 f) (checked on Schwarzschild):
    kappa1 = sigma1 [ (y^2 - sigma2^2) + gamma^2 (1 - c^2) ] / D^2 ,
    D = r1 r2 - e1 (e2 + gamma - gamma c),  r1 = m1 + sigma1,  y = l - sigma1 c ;
and kappa1 A1 = 4 pi sigma1 (Weyl rod).
Invariants: e1 (Gauss), A1, kappa1, M_tot = m1 + m2, Q_tot = e1 + e2.
Small hole: r+^2 = A1/4pi, mu = (r+^2 + e1^2)/(2 r+), kappa_loc = (r+^2 - e1^2)/(2 r+^3),
redshift N = kappa1/kappa_loc.  Test body (Bonnor) in RN(M_tot, Q_tot) at radius R with
N^2 = 1 - 2 M_tot/R + Q_tot^2/R^2 :   e/mu = (M_tot R - Q_tot^2)/(Q_tot R N).
Output: (e/mu)/(Bonnor value) = 1 + c1 + c2 + ...  (c_k of order eps^k).
Units: Q2 = 1 is NOT imposed; all symbols kept.
"""
import sys
import sympy as sp
eps = sp.Symbol('epsilon')
e, M2, Q2, l, c = sp.symbols('e M_2 Q_2 ell c', real=True)
K = int(sys.argv[1]) if len(sys.argv) > 1 else 3      # relative orders eps^0 .. eps^(K-1)


def tr(x, n):
    """truncate a polynomial in eps to degrees < n (after expansion)."""
    x = sp.expand(x)
    return sum(x.coeff(eps, k) * eps**k for k in range(0, n))


def inv(x, n):
    """1/x as a series in eps to degree < n; x = x0 + eps x1 + ..., x0 != 0."""
    x = sp.expand(x)
    x0 = x.coeff(eps, 0)
    d = sp.expand((x - x0) / x0)
    out, term = 1, 1
    for k in range(1, n):
        term = tr(-term * d, n)
        out += term
    return tr(out / x0, n)


def sqrt_s(x, n):
    x = sp.expand(x)
    x0 = x.coeff(eps, 0)
    d = sp.expand((x - x0) / x0)
    out, term = 1, 1
    for k in range(1, n):
        term = tr(term * d, n)
        out += sp.binomial(sp.Rational(1, 2), k) * term
    return tr(sp.sqrt(x0) * out, n)


N_ = K + 2
mm = sp.symbols(f'mm0:{N_}')
m1 = eps * sum(mm[k] * eps**k for k in range(N_))
e1 = eps * e
# gamma = (M2 e1 - m1 Q2)/(l + M2 + m1)
gam = tr((M2 * e1 - m1 * Q2) * inv(l + M2 + m1, N_ + 1), N_ + 1)
F = tr((e1 - gam) * (Q2 + gam) - m1 * M2, N_ + 1)
sol = {}
for k in range(N_ - 1):
    eq = sp.expand(F.coeff(eps, k + 1).subs(sol))
    sol[mm[k]] = sp.factor(sp.solve(eq, mm[k])[0])
    print(f'm1, eps^{k+1}:', sol[mm[k]], flush=True)
sub = lambda x: tr(sp.expand(sp.sympify(x).subs(sol)), N_)
m1s, gams = sub(m1), sub(gam)
s1sq_over = sp.expand(sub(m1s**2 - e1**2 + 2 * e1 * gams) / eps**2)    # sigma1^2 / eps^2
sigma1 = eps * sqrt_s(s1sq_over, N_ - 1)
s2sq = sub(M2**2 - Q2**2 - 2 * Q2 * gams)
y = l - sigma1 * c
r1 = m1s + sigma1
r2 = M2 + y
D = tr(r1 * r2 - e1 * (Q2 + gams - gams * c), N_)
kap1 = tr(sigma1 * ((y**2 - s2sq) + gams**2 * (1 - c**2)) * inv(tr(D**2, N_), N_), N_)
dc = [sp.simplify(sp.diff(kap1.coeff(eps, k), c)) for k in range(1, N_)]
print('kappa1 independent of c on the horizon:', all(v == 0 for v in dc), flush=True)
kap1 = sp.expand(kap1.subs(c, 0))
# r+^2 = sigma1 / kappa1 : kappa1 = eps^0 ... ; sigma1 = eps (...)
kap_over = sp.expand(kap1)                                 # starts at eps^0? check
k0 = kap_over.coeff(eps, 0)
print('kappa1 leading order eps^0 coefficient zero?', sp.simplify(k0) == 0, flush=True)
# kappa1 ~ 1/eps : multiply by eps
kap_e = sp.expand(kap1 * eps)                              # = eps * kappa1, starts at eps^0
rp2 = tr(sp.expand(sigma1 * eps) * inv(kap_e, N_), N_ + 1) # r+^2 ~ eps^2
rp = eps * sqrt_s(sp.expand(rp2 / eps**2), N_ - 1)
rp_over = sp.expand(rp / eps)
mu = tr(eps * (sp.expand(rp2 / eps**2) + e**2) * inv(2 * rp_over, N_ - 1), N_)
kloc_e = tr((sp.expand(rp2 / eps**2) - e**2) * inv(2 * rp_over**3, N_ - 1), N_ - 1)   # eps * kappa_loc
N = tr(kap_e * inv(kloc_e, N_ - 1), N_ - 1)
Mt, Qt = M2 + m1s, Q2 + e1
one_m = tr(1 - N**2, K + 1)
R = tr((Mt + sqrt_s(tr(Mt**2 - one_m * Qt**2, K + 1), K + 1)) * inv(one_m, K + 1), K + 1)
xi_test = tr((Mt * R - Qt**2) * inv(tr(Qt * R * N, K + 1), K + 1), K + 1)
mu_over = sp.expand(mu / eps)
ratio = tr(e * inv(tr(mu_over * xi_test, K), K), K)
print('mu leading:', sp.factor(mu_over.coeff(eps, 0)), flush=True)
for k in range(K):
    print(f'(e/mu)/Bonnor, eps^{k}:', sp.factor(ratio.coeff(eps, k)), flush=True)
