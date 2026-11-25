"""Solve the static l = 1 equations of l1_response.py by finite Laurent
polynomials, and identify the solution space (calc.tex Sec. 7 pilot)."""
import pickle
import sympy as sp
from l1_response import r, th, M, q, H0, H2, K, h1, a, b, c, h0, w

eqs = pickle.load(open('l1_eqs.pkl', 'rb'))
lo, hi = -5, 2


def laurent(name):
    cs = sp.symbols(f'{name}_m5:{hi - lo + 1}')   # placeholder names
    cs = [sp.Symbol(f'{name}[{k}]') for k in range(lo, hi + 1)]
    return sum(ck * r**k for ck, k in zip(cs, range(lo, hi + 1))), cs


def solve_sector(fields, exprs, subs_extra):
    ans, unk = {}, []
    for fld in fields:
        e, cs = laurent(str(fld.func))
        ans[fld] = e
        unk += cs
    sys_eqs = []
    for ex in exprs:
        ex = ex.subs(subs_extra)
        ex = ex.subs({fld: ans[fld] for fld in fields}).doit()
        ex = sp.simplify(ex / (sp.cos(th) if ex.has(sp.cos(th)) else 1))
        num = sp.numer(sp.together(sp.expand(ex.subs(th, sp.pi / 3))))
        num = sp.expand(num)
        sys_eqs += sp.Poly(num, r).coeffs()
    sol = sp.solve(sys_eqs, unk, dict=True)
    return ans, unk, sol


# axial
ans, unk, sol = solve_sector([h0, w], [eqs['Ea'][0, 3], eqs['Ma'][3]], {})
print('axial solutions:', len(sol))
s = sol[0]
for fld in (h0, w):
    print(' ', fld, '=', sp.factor(sp.expand(ans[fld].subs(s))))
print('  free:', [u for u in unk if u not in s and ans[h0].subs(s).has(u) or (u not in s and ans[w].subs(s).has(u))])

# polar, gauge K = h1 = 0, b = c = 0
pol = [eqs['Ep'][i, j] for i in range(4) for j in range(i, 4) if eqs['Ep'][i, j] != 0] + [m for m in eqs['Mp'] if m != 0]
ans, unk, sol = solve_sector([H0, H2, a], pol, {K: 0, h1: 0, b: 0, c: 0})
print('polar solutions:', len(sol))
s = sol[0]
for fld in (H0, H2, a):
    print(' ', fld, '=', sp.factor(sp.expand(ans[fld].subs(s))))
free = sorted({u for fld in (H0, H2, a) for u in ans[fld].subs(s).free_symbols if u in unk}, key=str)
print('  free:', free)

print('--- polar, no gauge fixing (b = c = 0) ---')
ans, unk, sol = solve_sector([H0, H2, K, h1, a], pol, {b: 0, c: 0})
print('polar solutions:', len(sol))
s = sol[0]
for fld in (H0, H2, K, h1, a):
    print(' ', fld, '=', sp.factor(sp.expand(ans[fld].subs(s))))
free = sorted({u for fld in (H0, H2, K, h1, a) for u in ans[fld].subs(s).free_symbols if u in unk}, key=str)
print('  free:', free)
pickle.dump((ans, s, free), open('l1_polar_sol.pkl', 'wb'))
