"""Locate the inconsistency of the order-3 near-zone equations in ckpt2_axi."""
import pickle, sys, sympy as sp
import __main__
import ckpt2_axi as C
__main__.Field = C.Field
import trunc as T
d = pickle.load(open('ckpt2_axi_o2.pkl', 'rb'))
dcols = d['log'][2]['free']
F = d['F']
m = 3
unk = []
for k in range(0, C.EPS + 1):
    h, A, syms = C.ansatz(m, k)
    F.add((k, m, 0), h, A)
    unk.extend(syms)
T.set_limits(C.EPS, m, 0)
g00, gss, A0 = C.rn_background(C.EPS)
G = T.mat(lambda a_, b: T.add(g00 if a_ == b == 0 else (gss if a_ == b else {}), F.h[a_][b]))
AA = [T.add(A0 if a_ == 0 else {}, F.A[a_]) for a_ in range(4)]
Ein, Mx = T.einstein_maxwell(G, AA)
names = [f'E{a_}{b}' for a_ in range(4) for b in range(a_, 4)] + [f'Mx{a_}' for a_ in range(4)]
comps = [Ein[a_][b] for a_ in range(4) for b in range(a_, 4)] + Mx
byN = {}
for nm, cpt in zip(names, comps):
    for N in range(C.EPS + 1):
        ex = cpt.get((N, m, 0), 0)
        if ex != 0:
            byN.setdefault(N, []).append((nm, C.reduce_eq(ex)))
acc = []
for N in sorted(byN):
    acc_new = acc + [e for nm, eqs in byN[N] for e in eqs]
    sol, cond, free = C.solve_ordered(acc_new, unk, dcols)
    print(f'N<={N}: {len(acc_new)} eqs, conditions:', cond, flush=True)
    if 1 in cond:
        sol0, _, _ = C.solve_ordered(acc, unk, dcols) if acc else ({}, [], [])
        for nm, eqs in byN[N]:
            res = [sp.factor(sp.sympify(e).subs(sol0)) for e in eqs]
            print('  ', nm, len(eqs), 'eqs; sample residuals:', [r for r in res if r != 0][:2])
        break
    acc = acc_new
