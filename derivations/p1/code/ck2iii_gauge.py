"""Class (iii), static axisymmetric electric case: gauge orbit of the order-3 slot-1 coefficients.
For static axisymmetric gauge vectors at near-zone order 3, slot 1 (length power 3):
   xi_i = x_i f(z, R) + zhat_i g(z, R),   f ~ z^a R^(2-a),  g ~ z^a R^(3-a),
the change of the order-3 fields is L_xi(RN) (keys: slots 1, 2, 3), expanded in the ansatz basis
of ckpt2_axi.  We (1) check that every a2-relation found by ckpt2_axi --light is invariant,
(2) list which slot-1 coefficients can be gauged to zero."""
import pickle, sympy as sp
import __main__
import ckpt2_axi as C
__main__.Field = C.Field
import trunc as T
from trunc import x, y, z, R
d = pickle.load(open('ckpt2_axi.pkl', 'rb'))
cond3 = d['log'][3]['cond']
m = 3
T.set_limits(C.EPS, 1, 0)
g00, gss, A0 = C.rn_background(C.EPS)
gb = [[(g00 if a == b == 0 else (gss if a == b else {})) for b in range(4)] for a in range(4)]
xs = (x, y, z); zh = (0, 0, 1)
gens = []
for a_ in range(0, 6):
    for kind in ('f', 'g'):
        P = 2 if kind == 'f' else 3
        if a_ > P + 3:
            continue
        mono = z**a_ * R**(P - a_)
        xi = [0] + [(xs[i] * mono if kind == 'f' else zh[i] * mono) for i in range(3)]
        gens.append((f'{kind}{a_}', xi))
results = {}
for name, xi in gens:
    # L_xi g = xi^c d_c g_ab + g_cb d_a xi^c + g_ac d_b xi^c ; xi carries slot 0, marker 1 -> keys (k,1,0)
    xiq = [T.S({(1, 1, 0): sp.sympify(v)}) if v != 0 else {} for v in xi]
    Lg = T.mat(lambda a, b: T.add(*[T.mul(xiq[c], T.d(gb[a][b], c)) for c in range(4)],
                                   *[T.mul(gb[c][b], T.d(xiq[c], a)) for c in range(4)],
                                   *[T.mul(gb[a][c], T.d(xiq[c], b)) for c in range(4)]))
    LA = [T.add(*[T.mul(xiq[c], T.d(A0 if b == 0 else {}, c)) for c in range(4)],
                *[T.mul(A0 if c == 0 else {}, T.d(xiq[c], b)) for c in range(4)]) for b in range(4)]
    # fit to the ansatz at each slot k (order 3): ansatz(m, k) - L = 0
    shift = {}
    for k in (1, 2, 3):
        h, A, syms = C.ansatz(m, k)
        eqs = []
        for a in range(4):
            for b in range(a, 4):
                eqs += C.reduce_eq(sp.expand(h[a][b] - Lg[a][b].get((k, 1, 0), 0)))
            eqs += C.reduce_eq(sp.expand(A[a] - LA[a].get((k, 1, 0), 0)))
        sol = sp.solve(eqs, syms, dict=True)
        if not sol:
            print(name, 'slot', k, ': L_xi(RN) not representable in the ansatz'); continue
        sol = sol[0]
        for s_ in syms:
            v = sp.simplify(sol.get(s_, 0))
            if v != 0:
                shift[s_] = v
    results[name] = shift
    print(name, 'shifts:', {str(k): v for k, v in shift.items()}, flush=True)
# invariance of the conditions that contain a2 (ho3s2_1)
a2s = sp.Symbol('ho3s2_1')
for ci, c in enumerate(cond3):
    num = sp.numer(sp.together(c))
    lin = [s_ for s_ in num.free_symbols if str(s_).startswith(('h', 'A', 'S')) and 'o3' in str(s_)]
    for name, sh in results.items():
        delta = sp.simplify(sum(sp.diff(num, s_) * sh.get(s_, 0) for s_ in lin))
        if delta != 0:
            print(f'condition {ci} NOT invariant under {name}: {delta}')
print('invariance check done')
pickle.dump(results, open('ck2iii_gauge.pkl', 'wb'))
