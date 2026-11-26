"""Route 2 (a): test body (m, q) moving slowly in RN(M, Q) near the MP point.
Q = M(1 + eta), q = m(1 - delta).  L = -m sqrt(f - rdot^2/f - r^2 Omega^2) + q Q/r  (A_t = -Q/r convention).
Slow motion: L = -V + (m/2) g_ij v^i v^j + O(v^4)."""
import sympy as sp
r, M, m, eta, delta, rd, Om, rho = sp.symbols('r M m eta delta rdot Omega rho', positive=True)
Q, q = M * (1 + eta), m * (1 - delta)
f = 1 - 2 * M / r + Q**2 / r**2
V = m * sp.sqrt(f) + q * Q / r            # like charges repel: V = m sqrt(f) + q Q / r
Vs = sp.series(sp.series(V, eta, 0, 2).removeO(), delta, 0, 2).removeO()
Vs = sp.expand(Vs).subs(eta * delta, 0)
Vr = sp.simplify(Vs.subs(r, rho + M).subs(sp.sqrt(rho**2), rho))
print('V at eta = delta = 0:', sp.simplify(Vr.subs({eta: 0, delta: 0})))
print('V - m to first order (isotropic rho = r - M):', sp.factor(sp.simplify(Vr - m)))
# kinetic metric at the extremal point, in isotropic rho = r - M
fe = (1 - M / r)**2
grr, gang = 1 / fe**sp.Rational(3, 2), r**2 / sp.sqrt(fe)
U = 1 + M / rho
print('g_rr * (drho)^2 at extremality  -> U^3 ?', sp.simplify(grr.subs(r, rho + M) - U**3))
print('angular: r^2/sqrt(f) vs rho^2 U^3 :', sp.simplify(gang.subs(r, rho + M) - rho**2 * U**3))
