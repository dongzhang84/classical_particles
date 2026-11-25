"""Re-solve the third-order l = 1 equations (saved by ckpt15_l1_order3.py) with
the slot-1 (rbar^2) coefficients -- first-order far-zone data a31 -- treated
as data, to expose the solvability condition on the second-order data."""
import pickle, sympy as sp
from ckpt15_l1 import solve_ordered
d = pickle.load(open('ckpt15_eqs3.pkl', 'rb'))
eqs, (V1, V2, V3), data3 = d['eqs'], d['V'], d['data3']
dsyms = sorted(set(d['dotted'].values()), key=str)
isS = lambda s: isinstance(s, sp.Symbol)
resp = [s for U in (V3, V2) for v in U.values() for s in v if isS(s) and s not in data3]
v1 = [s for v in V1.values() for s in v if isS(s)]
sol, cond, free = solve_ordered(eqs, resp + v1, data3 + dsyms)   # v1 first: see what is forced
print('--- V1 as responses: conditions', cond)
sol, cond, free = solve_ordered(eqs, resp, data3 + v1 + dsyms)
print('--- V1 as data: conditions on the data:')
for c in cond:
    print('  ', sp.factor(c))
print('free responses:', free)
