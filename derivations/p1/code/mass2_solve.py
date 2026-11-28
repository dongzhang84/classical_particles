"""Solve the projected eps^2 l=0 equations (mass2_diag_q*.pkl) with a Laurent
ansatz (plus log terms if needed); report free parameters and the mass function
dm(r) = (r/2) f R(r) at large r. Usage: python mass2_solve.py q"""
import pickle, sys, sympy as sp
r = sp.Symbol('r')
qv = sp.Rational(sys.argv[1]); Mv = 1
f = 1 - 2*Mv/r + qv**2/r**2
out = pickle.load(open(f'mass2_diag_q{str(qv).replace("/", "_")}.pkl', 'rb'))
P, R, b = [sp.Function(n)(r) for n in ('P', 'R', 'b')]
L = sp.log(r)
lo, hi = -12, 4
def ans(name, withlog):
    cs = [sp.Symbol(f'{name}{k}') for k in range(lo, hi + 1)]
    e = sum(ck*r**k for ck, k in zip(cs, range(lo, hi + 1)))
    if withlog:
        ls = [sp.Symbol(f'{name}L{k}') for k in range(lo, hi + 1)]
        e += L*sum(ck*r**k for ck, k in zip(ls, range(lo, hi + 1)))
        cs += ls
    return e, cs
for withlog in (False, True):
    eP, cP = ans('p', withlog); eR, cR = ans('R', withlog); eb, cb = ans('b', withlog)
    unk = cP + cR + cb
    eqs = []
    for k, v in out.items():
        num = sp.numer(sp.together(v))
        e = sp.expand(num.subs({P: eP, R: eR, b: eb}).doit())
        e = sp.expand(e)
        if withlog:
            e0 = e.subs(L, 0); e1 = sp.expand((e - e0)/L).subs(L, 0) if e.has(L) else 0
            parts = [sp.expand(e.coeff(L, 0)), sp.expand(e.coeff(L, 1)), sp.expand(e.coeff(L, 2))]
            parts = [sp.expand(e.subs(sp.log(r), sp.Symbol('LL'))).coeff(sp.Symbol('LL'), j) for j in range(3)]
        else:
            parts = [e]
        for pt in parts:
            pt = sp.expand(pt)
            if pt != 0:
                eqs += sp.Poly(pt, r, 1/r).coeffs() if pt.has(1/r) else sp.Poly(pt, r).coeffs()
    sol = sp.solve(eqs, unk, dict=True)
    print('withlog =', withlog, 'solutions:', len(sol))
    if sol:
        s = sol[0]
        Rs = sp.expand(eR.subs(s)); Ps = sp.expand(eP.subs(s)); bs = sp.expand(eb.subs(s))
        free = sorted({x for x in unk if x not in s and (Rs.has(x) or Ps.has(x) or bs.has(x))}, key=str)
        print('free:', free)
        dm = sp.expand(r/2*f*Rs)
        print('R =', Rs)
        print('P =', Ps)
        print('b =', bs)
        print('dm(r) =', sp.collect(dm, r))
        pickle.dump((Rs, Ps, bs, free), open(f'mass2_sol_q{str(qv).replace("/", "_")}.pkl', 'wb'))
        break
