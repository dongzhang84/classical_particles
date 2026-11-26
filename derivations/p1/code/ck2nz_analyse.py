"""Re-solve the class-(ii) projected equations with the slot-1 coefficients (first-order far-zone
r^2 terms) treated as data, to expose the a2 relation (as in ckpt15_analyse3.py)."""
import pickle, sympy as sp
from ckpt15_l1 import solve_ordered, vec, unknown_set, NAMES
d = pickle.load(open('ck2nz_ii.pkl', 'rb'))
eqs = d['eqs']
a2, e2, B2 = vec('A2'), vec('F2'), vec('C2')
syms = set().union(*[e.free_symbols for e in eqs])
k1 = sorted([s for s in syms if str(s).split('_')[0].endswith('k1')], key=str)
k23 = sorted([s for s in syms if str(s).split('_')[0].endswith(('k2', 'k3'))], key=str)
data = k1 + list(e2) + list(B2) + list(a2)
print('slot-1 symbols present:', k1)
sol, cond, free = solve_ordered(eqs, k23, data)
print('conditions:')
for c in cond:
    print('   ', sp.factor(c))
