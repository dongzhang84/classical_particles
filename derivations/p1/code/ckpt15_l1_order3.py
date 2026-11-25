"""Checkpoint 1.5, near-zone third order, l = 1 sector (see ckpt15_l1.py).

Third-order field = tbar * d/dt0 (order-2 solution) + general static l = 1
terms at rbar^2 (eps^1), rbar^1 (eps^2), rbar^0 (eps^3).  d/dt0 acts on the
order-2 data (a, E, B) and on its free parameters; dotted symbols are new
independent symbols.  The solvability conditions at eps^1..eps^3 are printed.
"""
import pickle, time, sys
import sympy as sp
import trunc as T
from trunc import t
from ckpt15_l1 import (vec, unknown_set, l1_form, put, empty, assemble, equations,
                       dot, cross, half, xs, M, q, solve_order2)

T.set_limits(3, 1, 0)
sol2, (a, E, B, U1, U2, U3), P2 = solve_order2()
syms2 = set().union(*[sp.sympify(v).free_symbols for U in (U1, U2, U3) for vv in U.values() for v in vv])
syms2 |= set(a) | set(E) | set(B)
free2 = sorted([s for s in syms2 if s not in sol2], key=str)
print('order-2 free parameters:', free2, flush=True)
dotted = {s: sp.Symbol(str(s) + 'd') for s in free2 + [s for s in a if s in sol2]}
def ddt(expr):
    """d/dt0 of an order-2 expression: replace each free symbol by its dot (linear)."""
    expr = sp.sympify(expr).subs(sol2)
    return sp.expand(sum(sp.diff(expr, s) * dotted[s] for s in free2))

P3 = empty()
for a_ in range(4):
    for slot_dict, tgt in [(P2[1][a_], ('A', a_))] + [(P2[0][a_][b], ('h', a_, b)) for b in range(4)]:
        for key, val in slot_dict.items():
            v = t * ddt(val)
            if v == 0:
                continue
            if tgt[0] == 'A':
                P3[1][a_] = T.add(P3[1][a_], {key: v})
            else:
                P3[0][tgt[1]][tgt[2]] = T.add(P3[0][tgt[1]][tgt[2]], {key: v})
V1, V2, V3 = unknown_set('k1'), unknown_set('k2'), unknown_set('k3')
a2, E2, B2 = vec('A2'), vec('F2'), vec('C2')     # second-order acceleration, uniform E, uniform B
data3 = a2 + E2 + B2
V2['h00'] = [-2 * c for c in a2]; V2['a0'] = E2; V2['aM'] = [half * c for c in B2]
for U, p, k in [(V1, 2, 1), (V2, 1, 2), (V3, 0, 3)]:
    h, A = l1_form(U, p)
    put(P3, k, h, A)
g, A = assemble(P3[0], P3[1])
t0 = time.time()
eqs = equations(g, A, 3)
print(f'order 3: {len(eqs)} equations, {time.time()-t0:.0f}s', flush=True)
from ckpt15_l1 import solve_ordered
pickle.dump({'eqs': eqs, 'V': (V1, V2, V3), 'data3': data3, 'dotted': dotted, 'sol2': sol2},
            open('ckpt15_eqs3.pkl', 'wb'))
dsyms = sorted(set(dotted.values()), key=str)
unk = [s for U in (V3, V2, V1) for v in U.values() for s in v if s not in data3 and isinstance(s, sp.Symbol)]
sol3, cond, free3 = solve_ordered(eqs, unk, data3 + dsyms)
print('order-3 conditions on the data:')
for c in cond:
    print('  ', c)
print('order-3 free:', free3)
sol3 = []
if sol3:
    S = sol3[0]
    print('order 3 solvable; free:', [s for s in unk if s not in S])
    for tag, U in [('r^2', V1), ('r^1', V2), ('r^0', V3)]:
        for k, v in U.items():
            vals = [sp.factor(S.get(s, s)) for s in v]
            if any(vv != 0 for vv in vals):
                print(tag, k, vals)
    cond = []
else:
    A_, b_ = sp.linear_eq_to_matrix(eqs, unk)
    aug = A_.row_join(b_)
    rref, piv = aug.rref()
    cond = [sp.factor(rref[i, -1]) for i in range(rref.rows)
            if all(rref[i, j] == 0 for j in range(len(unk))) and rref[i, -1] != 0]
    print('conditions on the data:')
    for c in cond:
        print('  ', c)
pickle.dump({'sol2': {str(k): v for k, v in sol2.items()}, 'cond': cond,
             'sol3': {str(k): v for k, v in (sol3[0] if sol3 else {}).items()}},
            open('ckpt15_order3.pkl', 'wb'))
