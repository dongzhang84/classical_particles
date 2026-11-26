"""MP local data in the isotropic (MP harmonic) coordinates of the static observer:
X = U0 y, t_F = t/U0;  g00 = -U0^2/U^2, A0 = U0 (1/U - 1/U0), g_ij = (U/U0)^2 delta_ij.
Fit to the near-zone order-2 slot-0 solution of ckpt2_axi --light --iso and evaluate Bal_2 (q = M)."""
import pickle, sys, sympy as sp
import __main__
import ckpt2_axi as C
__main__.Field = C.Field
from trunc import x, y, z, R
X1, X3 = sp.symbols('X1 X3')
Lv = sp.Rational([a for a in sys.argv[1:] if not a.startswith('--')][0])
eps = sp.Symbol('eps')
U = lambda xx, zz: 1 + 1 / sp.sqrt(xx**2 + (zz - Lv)**2)
U0 = U(0, 0)
yx, yz = X1 / U0, X3 / U0
def ser2(e):
    return sp.expand(sp.series(e.subs({X1: eps * X1, X3: eps * X3}), eps, 0, 3).removeO().subs(eps, 1))
g00 = ser2(-U0**2 / U(yx, yz)**2)
A0 = ser2(U0 * (1 / U(yx, yz) - 1 / U0))
gxx = ser2((U(yx, yz) / U0)**2)
E = sp.diff(A0, X3).subs({X1: 0, X3: 0}); a = -sp.diff(g00, X3).subs({X1: 0, X3: 0}) / 2
print('a =', a, ' E =', E)
quad = lambda e: sp.expand(sum(t_ for t_ in sp.Add.make_args(e) if sp.Poly(t_, X1, X3).total_degree() == 2))
h00q, A0q, hq = quad(g00 + 1), quad(A0), quad(gxx - 1)
# near-zone order-2 slot-0 fields (q = M), with kept params unknown
d = pickle.load(open('ckpt2_axi_o2_iso_conf' + ('_isoall' if '--isoall' in sys.argv else '') + '.pkl', 'rb'))
F = d['F']; key = (0, 2, 0)
kept = [sp.Symbol('Ao2s0_2'), sp.Symbol('S1o2s0_2')]
M = C.M
sub0 = {sp.Symbol('E'): E, C.q: M, sp.Symbol('Ao2s1_1'): 0}
nz = {'h00': F.h[0][0].get(key, 0), 'A0': F.A[0].get(key, 0), 'hxx': F.h[1][1].get(key, 0), 'hzz': F.h[3][3].get(key, 0)}
toX = lambda e: sp.expand(sp.sympify(e).subs(sub0).subs({y: 0}).subs({R: sp.sqrt(X1**2 + X3**2)}).subs({x: X1, z: X3}))
eqs = []
for nm, mpv in [('h00', h00q), ('A0', A0q), ('hxx', hq), ('hzz', hq)]:
    eqs += sp.Poly(sp.expand(toX(nz[nm]) - mpv), X1, X3).coeffs()
sol = sp.solve(eqs, kept, dict=True)
print('fit of kept parameters (consistency: solvable =', bool(sol), '):', sol)
pickle.dump({'E': E, 'sol': sol[0] if sol else None}, open(f'ck2iii_mpiso_L{Lv}.pkl', 'wb'))
