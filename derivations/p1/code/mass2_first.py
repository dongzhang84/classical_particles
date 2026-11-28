"""Step 1 of the effective-mass check: pick a gauge for the exact static l=1
first-order solution (RN in uniform E with acceleration Ma=qE), so that it is
as regular as possible, and print it."""
import pickle, sympy as sp
from l1_response import r, M, q, H0, H2, K, h1, a
ans, s, free = pickle.load(open('l1_polar_sol.pkl', 'rb'))
E = sp.Symbol('E')
sym = lambda n, k: sp.Symbol(f'{n}[{k}]')
choice = {sym('K', 2): 0, sym('a', 1): E, sym('a', 0): 0, sym('a', -1): 0, sym('a', -2): 0, sym('a', -3): 0,
          sym('K', 0): 0, sym('h1', -1): 0, sym('h1', -2): 0, sym('h1', -3): 0, sym('h1', -4): 0, sym('h1', 0): 0}
# remove the r^1 term of K (large translation) with h1[2]
K1 = sp.expand(ans[K].subs(s).coeff(r, 1).subs(choice))
h12 = sp.solve(K1.subs(sym('h1', 2), sp.Symbol('X')), sp.Symbol('X'))[0]
choice[sym('h1', 2)] = h12
sol = {fld: sp.factor(sp.expand(ans[fld].subs(s).subs(choice))) for fld in (H0, H2, K, h1, a)}
for k, v in sol.items():
    print(k, '=', v)
pickle.dump(sol, open('mass2_first.pkl', 'wb'))
