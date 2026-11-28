"""Compare quasi-local mass definitions at second order (static slice, spheres r = const):
 (a) from the l=0 part of g_rr only:  dm_a = (r/2) f R
 (b) angle-averaged Misner-Sharp with g^rr to eps^2
 (c) Hawking mass m_H = sqrt(A/16pi) (1 - (1/16pi) oint k^2 dA), k = mean curvature in the slice.
Usage: python mass2_defs.py q"""
import pickle, sys, sympy as sp
from l1_response import r, th, M as Ms, q as qs, H0, H2, K, h1, a
qv = sp.Rational(sys.argv[1]); Mv = 1
G, rest = pickle.load(open('mass2_gauge.pkl', 'rb'))
Esym = sp.Symbol('E')
G = {k: sp.factor(v.subs({u: 0 for u in rest}).subs({Ms: Mv, qs: qv, Esym: 1})) for k, v in G.items()}
Rs, Ps, bs, free = pickle.load(open(f'mass2_sol_q{str(qv).replace("/", "_")}.pkl', 'rb'))
Rs = Rs.subs(sp.Symbol('r'), r)
eps = sp.Symbol('epsilon'); c = sp.cos(th)
f = 1 - 2*Mv/r + qv**2/r**2
B = (1 + eps*G[H2]*c + eps**2*Rs)/f
Cf = r**2*(1 + eps*G[K]*c)
avg = lambda e: sp.integrate(sp.expand(e)*sp.sin(th), (th, 0, sp.pi))/2
def e2(x): return sp.diff(x, eps, 2).subs(eps, 0)/2
# (a)
dma = sp.expand(r/2*f*Rs)
# (b) g^rr = 1/B ; m = r/2 (1 - <g^rr>)
grr = 1/B
dmb = sp.expand(sp.simplify(-r/2*avg(e2(grr))))
# (c) Hawking: area A = oint Cf sin = 4 pi r^2 (exact); k = (1/sqrt(B)) d_r ln Cf
k = (1/sp.sqrt(B))*sp.diff(Cf, r)/Cf
integrand = k**2*Cf/(16*sp.pi)*2*sp.pi       # (1/16pi) oint k^2 dA, phi integrated
I2 = sp.simplify(sp.integrate(sp.expand(e2(integrand)*sp.sin(th)), (th, 0, sp.pi)))
dmc = sp.expand(sp.simplify(-r/2*I2))       # sqrt(A/16pi) = r/2
for name, dm in (('a', dma), ('b', dmb), ('c', dmc)):
    ser = sp.series(dm.subs(r, 1/sp.Symbol('u')), sp.Symbol('u'), 0, 2).removeO()
    print(name, 'large-r:', sp.collect(sp.expand(ser.subs(sp.Symbol('u'), 1/r)), r))
pred = -sp.Rational(1, 20)*(Mv**2 - qv**2)**2/Mv
print('AB-fit prediction for constant term:', pred, float(pred))
