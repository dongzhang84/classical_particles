"""Checkpoint 2, near-zone verification, class (ii): the part of the third-order
near-zone l = 1 solvability condition that is BILINEAR in
    background uniform data  (E_b, B_b; acceleration a_0 = q E_b / M)   [near-zone order 1, slot 0]
and first-order far-zone data (f, b; acceleration a_1 = q f / M)        [near-zone order 2, slot 1].

Structure (P gauge, see calc.tex sec:L6ck2plan):
  order-1 near-zone field U_b : RN in the background uniform field + acceleration (l = 1),
  order-2 near-zone field U_1 : the checkpoint-1.5 order-2 solution (l = 1), independent of the background,
  order-3 bilinear field  U_3 : general l = 1 ansatz (both parities) at slots 1, 2, 3 (powers 2, 1, 0),
                                with the second-order data a2 (h00 = -2 a2.x), e2 (A0 = e2.x), B2 at slot 2.
trunc keys: (eps, mu1, mu2) = (powers of M,q incl. data, degree in background data, degree in first-order data).
The order-3 bilinear equations sit at keys (N, 1, 1); they are projected on l = 1 (electric and magnetic parts of
every component) and solved with the data columns last.  Test-body prediction (Z = 0 in H fields): no bilinear
term in the a2 relation, i.e. X_(ii) = 0.
"""
import sys, time, pickle
import sympy as sp
import trunc as T
from trunc import t, x, y, z, R, nn
from ckpt15_l1 import (vec, unknown_set, l1_form, put, empty, dot, cross, half, xs, M, q,
                       solve_ordered, rn_background, NAMES)
from series import DIRECTIONS

EPS = 3


def keyed(field, marker):
    """Re-key an 'empty()'-type field built with keys (slot, 1, 0) to the given marker."""
    out = empty()
    for a in range(4):
        out[1][a] = {(k[0],) + marker: v for k, v in field[1][a].items()}
        for b in range(4):
            out[0][a][b] = {(k[0],) + marker: v for k, v in field[0][a][b].items()}
    return out


def eqs_points(cpt_list, key):
    """Equations from components at a given trunc key, by evaluating on rational directions."""
    eqs = []
    for ex in cpt_list:
        ex = ex.get(key, 0)
        if ex == 0:
            continue
        for nv in DIRECTIONS:
            e = sp.expand(sp.sympify(ex).subs({x: R * nv[0], y: R * nv[1], z: R * nv[2]}, simultaneous=True))
            e = sp.expand(e * R**12)
            eqs += [cf for cf in sp.Poly(e, R, t).coeffs() if cf != 0]
    return eqs


def uniform_solution(slot0, tag, Eh, Bh, ah):
    """RN in a uniform field + acceleration, l = 1, data at slot slot0, responses at slot0+1..EPS.
    Returns (solution dict, field) with field keys (slot, 1, 0)."""
    P = empty()
    dh = [[0] * 4 for _ in range(4)]
    dh[0][0] = -2 * dot(ah, xs)
    put(P, slot0, dh, [dot(Eh, xs)] + [half * cross(Bh, xs)[i] for i in range(3)])
    Us = []
    for k in range(slot0, EPS + 1):
        U = unknown_set(f'{tag}{k}')
        if k == slot0:
            for nm in ('h00', 'a0', 'aM'):
                U[nm] = [0, 0, 0]
        h, A = l1_form(U, 1 + slot0 - k)
        put(P, k, h, A)
        Us.append(U)
    T.set_limits(EPS, 1, 0)
    g00, gss, A0 = rn_background(EPS)
    G = T.mat(lambda a, b: T.add(g00 if a == b == 0 else (gss if a == b else {}), P[0][a][b]))
    AA = [T.add(A0 if a == 0 else {}, P[1][a]) for a in range(4)]
    Ein, Mx = T.einstein_maxwell(G, AA)
    comps = [Ein[a][b] for a in range(4) for b in range(a, 4)] + Mx
    eqs = []
    for N in range(EPS + 1):
        eqs += eqs_points(comps, (N, 1, 0))
    unk = [s for U in reversed(Us) for v in U.values() for s in v if s != 0]
    sol, cond, free = solve_ordered(eqs, unk, list(ah) + list(Eh) + list(Bh))
    print(f'  [{tag}] conditions:', cond, flush=True)
    ca = sp.solve(cond, ah, dict=True)[0] if cond else {}
    zero = {s: 0 for s in free}                     # P gauge: free (gauge / homogeneous) parameters = 0
    sol = {k: sp.factor(sp.sympify(v).subs(ca).subs(zero)) for k, v in sol.items()}
    sol.update(ca); sol.update(zero)
    F = empty()
    for a in range(4):
        F[1][a] = {k: sp.expand(sp.sympify(v).subs(sol)) for k, v in P[1][a].items()}
        F[1][a] = {k: v for k, v in F[1][a].items() if v != 0}
        for b in range(4):
            F[0][a][b] = {k: sp.expand(sp.sympify(v).subs(sol)) for k, v in P[0][a][b].items()}
            F[0][a][b] = {k: v for k, v in F[0][a][b].items() if v != 0}
    return sol, F


# ------------------------------------------------------------------ l = 1 projections
def proj_scalar(coef):
    return T.l1_scalar(coef)


def proj_vector(coef_vec):
    P, C = T.l1_vector(coef_vec)
    W = T.l1_magnetic_vector(coef_vec)
    return P + C + W


def proj_tensor(coef_mat):
    D, E, F = T.l1_tensor(coef_mat)
    nT = [sum(nn[i] * coef_mat[i][j] for i in range(3)) for j in range(3)]
    W = [2 * w for w in T.l1_magnetic_vector(nT)]
    return D + E + F + W


def to_n(expr, p, d=0):
    """coefficient of R^p t^d after x -> R n."""
    e = sp.expand(expr)
    e = e.coeff(t, d) if d else e.subs(t, 0)
    e = sp.expand(e.subs({x: R * nn[0], y: R * nn[1], z: R * nn[2]}, simultaneous=True))
    out = 0
    for term in sp.Add.make_args(e):
        c, pw = term.as_coeff_exponent(R)
        if pw == p:
            out += c
    return sp.expand(out)


def projected_equations(Ein, Mx, key, powers):
    eqs, tags = [], []
    for p in powers:
        comps = []
        comps.append(('E00', proj_scalar(to_n(Ein[0][0].get(key, 0), p))))
        comps.append(('Mx0', proj_scalar(to_n(Mx[0].get(key, 0), p))))
        comps.append(('E0i', proj_vector([to_n(Ein[0][i + 1].get(key, 0), p) for i in range(3)])))
        comps.append(('Mxi', proj_vector([to_n(Mx[i + 1].get(key, 0), p) for i in range(3)])))
        comps.append(('Eij', proj_tensor([[to_n(Ein[i + 1][j + 1].get(key, 0), p) for j in range(3)] for i in range(3)])))
        for nm, lst in comps:
            for e in lst:
                e = sp.expand(e)
                if e != 0:
                    eqs.append(e); tags.append((nm, key[0], p))
    return eqs, tags


if __name__ == '__main__':
    t0 = time.time()
    Eb, Bb, ab = vec('Eb'), vec('Bb'), vec('ab')
    f1, b1, a1 = vec('f'), vec('b'), vec('aa')
    print('background solution (order 1, data at slot 0):', flush=True)
    solb, Fb = uniform_solution(0, 'bg', Eb, Bb, ab)
    print('first-order-data solution (order 2, data at slot 1):', flush=True)
    sol1, F1 = uniform_solution(1, 'fo', f1, b1, a1)
    print(f'  done, {time.time()-t0:.0f}s', flush=True)
    pickle.dump({'solb': solb, 'sol1': sol1}, open('ck2nz_lower.pkl', 'wb'))
    # order-3 bilinear unknowns
    P3 = empty()
    a2, e2, B2 = vec('A2'), vec('F2'), vec('C2')
    V = {}
    for k, p in [(1, 2), (2, 1), (3, 0)]:
        U = unknown_set(f'k{k}')
        if k == 2:
            U['h00'] = [-2 * c for c in a2]; U['a0'] = e2; U['aM'] = [half * c for c in B2]
        h, A = l1_form(U, p)
        put(P3, k, h, A)
        V[k] = U
    T.set_limits(EPS, 1, 1)
    g00, gss, A0 = rn_background(EPS)
    Fb_k, F1_k, P3_k = keyed(Fb, (1, 0)), keyed(F1, (0, 1)), keyed(P3, (1, 1))
    G = T.mat(lambda a, b: T.add(g00 if a == b == 0 else (gss if a == b else {}),
                                 Fb_k[0][a][b], F1_k[0][a][b], P3_k[0][a][b]))
    AA = [T.add(A0 if a == 0 else {}, Fb_k[1][a], F1_k[1][a], P3_k[1][a]) for a in range(4)]
    t1 = time.time()
    Ein, Mx = T.einstein_maxwell(G, AA)
    print(f'order-3 Einstein-Maxwell: {time.time()-t1:.0f}s', flush=True)
    eqs, tags = [], []
    for N in range(EPS + 1):
        e_, t_ = projected_equations(Ein, Mx, (N, 1, 1), range(-6, 3))
        eqs += e_; tags += t_
    print(f'{len(eqs)} projected equations', flush=True)
    data = list(e2) + list(B2) + list(a2)
    unk = [s for k in (3, 2, 1) for v in V[k].values() for s in v
           if isinstance(s, sp.Symbol) and s not in data]
    sol3, cond, free3 = solve_ordered(eqs, unk, data)
    print('class (ii) solvability conditions on the second-order data:')
    for c in cond:
        print('   ', c)
    pickle.dump({'eqs': eqs, 'tags': tags, 'cond': cond, 'sol3': sol3, 'free3': free3},
                open('ck2nz_ii.pkl', 'wb'))
    print(f'total {time.time()-t0:.0f}s')
