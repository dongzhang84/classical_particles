"""MP check of X_(iii): local data of an extremal RN background (M_b = Q_b = 1) in the Fermi
(radial) coordinates of the static observer at distance L, mapped to the ansatz parameters of
ckpt2_axi (order-1 slot 0: h00 = -2 a z, A0 = E z; order-2 slot 0: h00, A0, h_ij quadratic in
radial gauge), then X evaluated with an extremal body (q = M)."""
import sympy as sp, pickle
import sys
L = sp.Rational(sys.argv[1]) if len(sys.argv) > 1 else sp.Integer(3)
X1, X2, X3 = sp.symbols('X1 X2 X3')        # Fermi (Riemann normal) spatial coordinates
y = sp.symbols('y1:4')
Xv = (X1, X2, X3)
Uf = lambda yy: 1 + 1 / sp.sqrt(yy[0]**2 + yy[1]**2 + (yy[2] - L)**2)
U0 = Uf((0, 0, 0))
dU = [sp.diff(Uf(y), y[i]).subs({y[0]: 0, y[1]: 0, y[2]: 0}) for i in range(3)]
# Christoffels of h = U^2 delta at the origin
G = [[[((1 if i == j else 0) * dU[k] + (1 if i == k else 0) * dU[j] - (1 if j == k else 0) * dU[i]) / U0
       for k in range(3)] for j in range(3)] for i in range(3)]
Gy = lambda i, j, k: ((1 if i == j else 0) * sp.diff(Uf(y), y[k]) + (1 if i == k else 0) * sp.diff(Uf(y), y[j])
                      - (1 if j == k else 0) * sp.diff(Uf(y), y[i])) / Uf(y)
at0 = {y[0]: 0, y[1]: 0, y[2]: 0}
dG = [[[[sp.simplify(sp.diff(Gy(i, j, k), y[l]).subs(at0)) for l in range(3)] for k in range(3)] for j in range(3)] for i in range(3)]
v0 = [Xv[i] / U0 for i in range(3)]
ymap = [v0[i] - sp.Rational(1, 2) * sum(G[i][j][k] * v0[j] * v0[k] for j in range(3) for k in range(3))
        + sp.Rational(1, 6) * sum((-dG[i][j][k][l] + 2 * sum(G[i][m][k] * G[m][j][l] for m in range(3)))
                                  * v0[l] * v0[j] * v0[k] for j in range(3) for k in range(3) for l in range(3))
        for i in range(3)]
eps = sp.Symbol('eps')
def ser2(expr):
    e = expr.subs({X1: eps * X1, X2: eps * X2, X3: eps * X3})
    return sp.expand(sp.series(e, eps, 0, 3).removeO().subs(eps, 1))
g00 = ser2(-U0**2 / Uf(ymap)**2)          # t_F = t / U0
for sgn in (1, -1):
    At = lambda yy: sgn / Uf(yy)
    A0 = ser2(U0 * (At(ymap) - At((0, 0, 0))))
    a = -sp.diff(g00, X3).subs({X1: 0, X2: 0, X3: 0}) / 2        # h00 = -2 a z at linear order
    E = sp.diff(A0, X3).subs({X1: 0, X2: 0, X3: 0})
    print('sign', sgn, ': a =', sp.simplify(a), ' E =', sp.simplify(E), ' a/E =', sp.simplify(a / E))
    if sp.simplify(a / E - 1) == 0:
        break
z_, R_ = sp.symbols('z R')
quad = lambda e: sp.expand(sum(term for term in sp.Add.make_args(sp.expand(e)) if sp.Poly(term, X1, X2, X3).total_degree() == 2))
h00q, A0q = quad(g00 + 1), quad(A0)
# axisymmetric: with X1 = x, X2 = 0 (y = 0), X3 = z, R^2 = x^2 + z^2 -> coefficients of R^2 and z^2
def axi(e):
    e = sp.expand(e.subs({X2: 0}))
    cx = e.coeff(X1, 2); cz = e.subs(X1, 0).coeff(X3, 2)
    return sp.simplify(cx), sp.simplify(cz - cx)       # e = cx R^2 + (cz - cx) z^2
print('order-2 slot-0 h00: R^2, z^2 coefficients:', axi(h00q))
print('order-2 slot-0 A0 : R^2, z^2 coefficients (Ao2s0_0, Ao2s0_2):', axi(A0q))
# spatial metric pulled back, quadratic part
J = [[sp.diff(ymap[i], Xv[a_]) for a_ in range(3)] for i in range(3)]
gij = [[ser2(sum(Uf(ymap)**2 * J[i][a_] * J[i][b] for i in range(3))) for b in range(3)] for a_ in range(3)]
hq = [[quad(gij[a_][b] - (1 if a_ == b else 0)) for b in range(3)] for a_ in range(3)]
print('radial gauge check x^j h_ij:', [sp.simplify(sum(Xv[j] * hq[i][j] for j in range(3))) for i in range(3)])
# fit h_ij = delta S1 + x x S2 + (x zh + zh x) S3 + zh zh S4 on the plane X2 = 0
s = sp.symbols('s1R s1z s2 s3 s4R s4z')
Xa, Za = X1, X3
S1 = s[0] * (Xa**2 + Za**2) + s[1] * Za**2
S2 = s[2]; S3 = s[3] * Za; S4 = s[4] * (Xa**2 + Za**2) + s[5] * Za**2
xs = (Xa, 0, Za); zh = (0, 0, 1)
eqs = []
for i in range(3):
    for j in range(3):
        model = (1 if i == j else 0) * S1 + xs[i] * xs[j] * S2 + (xs[i] * zh[j] + zh[i] * xs[j]) * S3 + zh[i] * zh[j] * S4
        eqs += sp.Poly(sp.expand(model - hq[i][j].subs(X2, 0)), Xa, Za).coeffs()
fit = sp.solve(eqs, s, dict=True)
print('h_ij fit:', fit)
pickle.dump({'a': a, 'E': E, 'h00': axi(h00q), 'A0': axi(A0q), 'hij': fit}, open(f'ck2iii_mp_L{L}.pkl', 'wb'))
