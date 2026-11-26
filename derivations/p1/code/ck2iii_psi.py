"""Class (iii), static electric: the second-order balance with the MP-linear split (rule 3'),
computed from the NEAR-ZONE solution alone (no exact solution used).

Regular parts are taken in psi1 = (-g00)^(-1/2), psi2 = A0 psi1:
  value of the far-zone order-k field at gamma  = l=0 part of the R^0 coefficient at near-zone key (k, k),
  gradient (d_z) of the order-k field at gamma   = l=1 part of the R^1 coefficient at key (k, k+1).
Test-body static balance (body potential A0 = -q/r, so Phi = -A0):  M d psi1 + q (psi1 d psi2 - psi2 d psi1) = 0.
Residual Bal_k at each order k; Bal_0 = 0 is the zeroth-order balance.  Expected: Bal_2 = 0 on MP data
(extremal body in extremal RN background); otherwise Bal_2 is the residual force in this convention."""
import pickle, sys, sympy as sp
import __main__
import ckpt2_axi as C
__main__.Field = C.Field
import trunc as T
from trunc import x, y, z, R, nn
M, q = C.M, C.q
SUF = ('_iso' if '--iso' in sys.argv else '_radial') + ('_conf' if '--conf' in sys.argv else '') + ('_isoall' if '--isoall' in sys.argv else '')
d = pickle.load(open(f'ck2iii_struct{SUF}.pkl', 'rb'))
F, sol, cond, data, free = d['F'], d['sol'], d['cond'], d['data'], d['free']
A2, F2 = sp.symbols('A2 F2')
hE, hF = sp.Symbol('hEk1'), sp.Symbol('hFk1')
gauge = {hE: 0, hF: 0}
# solve the conditions for the determined slot-1 coefficients and a2 (data columns), keep F2 (e2)
conds = [sp.numer(sp.together(c)).subs(gauge) for c in cond]
dsyms = [s_ for s_ in data if s_ not in (F2,) and s_ not in gauge]
csol = sp.solve(conds, [s_ for s_ in dsyms if any(c.has(s_) for c in conds)], dict=True)[0]
print('conditions solved for:', list(csol.keys()))
full = {}
full.update(gauge); full.update(csol)
for s_ in dsyms:
    if s_ not in full:
        full[s_] = 0                       # remaining data (vector/magnetic slot-1 structures): zero in this class
for s_ in free:
    full[s_] = 0                           # free (gauge/homogeneous) responses: zero
resp = {k: sp.sympify(v).subs(full) for k, v in sol.items()}
full.update(resp)
full[sp.Symbol('Ao2s1_1')] = 0
Fs = F.subs({k: sp.sympify(v).subs(full) for k, v in full.items()})
# total g00 and A0 as trunc series (eps, lambda)
T.set_limits(C.EPS, 3, 0)
g00, gss, A0 = C.rn_background(C.EPS)
G00 = T.add(g00, Fs.h[0][0])
AA0 = T.add(A0, Fs.A[0])
u = T.add(T.scal(-1, G00), T.const(-1))          # -g00 - 1
psi1, term = T.const(1), T.const(1)
coef = 1
for n in range(1, 7):
    coef = coef * (-sp.Rational(1, 2) - (n - 1)) / n     # binomial(-1/2, n)
    term = T.mul(term, u)
    psi1 = T.add(psi1, T.scal(coef, term))
psi2 = T.mul(AA0, psi1)


def l0(expr):
    return T.sphere_avg(T.radial_part(expr, 0))


def l1z(expr):
    return 3 * T.sphere_avg(T.radial_part(expr, 1) * nn[2])


vals = {}
for name, P in (('psi1', psi1), ('psi2', psi2)):
    for k in range(0, 3):
        vals[(name, 'v', k)] = sp.expand(l0(P.get((k, k, 0), 0)))
        vals[(name, 'g', k)] = sp.expand(l1z(P.get((k, k + 1, 0), 0)))
for k_, v in vals.items():
    print(k_, sp.factor(v))
Bal = {}
for k in range(0, 3):
    s_ = 0
    for i in range(0, k + 1):
        j = k - i
        s_ += q * (vals[('psi1', 'v', i)] * vals[('psi2', 'g', j)] - vals[('psi2', 'v', i)] * vals[('psi1', 'g', j)])
    s_ += M * vals[('psi1', 'g', k)]          # balance: M d psi1 + q (psi1 d psi2 - psi2 d psi1) = 0  (Phi = -A0)
    Bal[k] = sp.factor(sp.expand(s_))
    print(f'Bal_{k} =', Bal[k])
pickle.dump({'vals': vals, 'Bal': Bal}, open('ck2iii_psi.pkl', 'wb'))

for Lv in [a for a in sys.argv[1:] if not a.startswith('--')]:
    if '--conf' in sys.argv:
        mp = pickle.load(open(f'ck2iii_mpiso_L{Lv}.pkl', 'rb'))
        sub = {sp.Symbol('E'): mp['E']}; sub.update(mp['sol'])
    else:
        mp = pickle.load(open(f'ck2iii_mp_L{Lv}.pkl', 'rb'))
        sub = {sp.Symbol('E'): mp['E'], sp.Symbol('Ao2s0_2'): mp['A0'][1], sp.Symbol('S4o2s0_0'): mp['hij'][0][sp.Symbol('s4R')]}
    print(f'L = {Lv}: Bal_2 (MP, q = M) =', sp.simplify(Bal[2].subs(sub).subs(q, M)))
