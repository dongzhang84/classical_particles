"""Find, in the exact static l=1 family (l1_polar_sol.pkl), a gauge that is
(i) regular on the future horizon: H0 + H2 divisible by D = r^2 - 2Mr + q^2,
(ii) asymptotically standard near-zone form: H2, K, h1 have no growing terms,
     a ~ E r, H0 ~ -2 a r with M a = q E."""
import pickle, sympy as sp
from l1_response import r, M, q, H0, H2, K, h1, a
ans, s, free = pickle.load(open('l1_polar_sol.pkl', 'rb'))
E = sp.Symbol('E')
F = {fld: sp.expand(ans[fld].subs(s)) for fld in (H0, H2, K, h1, a)}
fr = [u for u in free]
conds = [F[a].coeff(r, 1) - E]
for fld, powers in ((H2, (1, 2)), (K, (1, 2)), (h1, (1, 2))):
    conds += [F[fld].coeff(r, p) for p in powers]
D = r**2 - 2 * M * r + q**2
rem = sp.rem(sp.numer(sp.together(F[H0] + F[H2])) , D, r)
conds += sp.Poly(sp.expand(rem), r).coeffs()
sol = sp.solve(conds, fr, dict=True)
print('solutions:', len(sol))
for so in sol:
    rest = [u for u in fr if u not in so]
    print('fixed:', {k: sp.factor(v) for k, v in so.items()})
    print('remaining free:', rest)
    G = {fld: sp.factor(F[fld].subs(so)) for fld in F}
    for k, v in G.items():
        print(' ', k, '=', v)
    pickle.dump((G, rest), open('mass2_gauge.pkl', 'wb'))
