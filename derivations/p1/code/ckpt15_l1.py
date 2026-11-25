"""Checkpoint 1.5 (vacuum background, charged black hole): the l = 1 sector of
the near zone through third order, in P gauge.

P gauge: the far-zone fields are expressed about the background geodesic gamma
with the body mass-centred on it; the near-zone first-order field then vanishes
(up to monopole redefinitions), so the second- and third-order near-zone fields
obey the LINEARIZED Einstein-Maxwell equations on RN.  The motion is encoded in
the l = 1 data: the acceleration potential h_00 = -2 a.x and the uniform fields.

Bookkeeping: slot eps counts powers of (M, q), including those carried by the
far-zone data (first-order far-zone fields are proportional to M, q).  At
near-zone order m every term has (length power) + (eps power) = m, where tbar
counts as a length.  The equations at total eps-order N then sit at
rbar^{-N} (order 2) or rbar^{1-N} (order 3).

Order 2 (static): data at eps^1, rbar^1:
    h_00 = -2 a.x,   A_0 = E.x,   A_i = (B x x)_i / 2,
plus the most general l = 1 terms at rbar^0 (eps^2) and rbar^-1 (eps^3), and
general l = 1 terms at rbar^1 (eps^1) besides the data (gauge).
Order 3: tbar * d/dt0 (order-2 solution) + general static l = 1 terms at
rbar^2 (eps^1), rbar^1 (eps^2), rbar^0 (eps^3).
"""
import pickle, sys, time
import sympy as sp
import trunc as T
from trunc import t, x, y, z, R
from series import DIRECTIONS

M, q = sp.symbols('M q')
xs = (x, y, z)
n = [xi / R for xi in xs]
dot = lambda u, v: sum(ui * vi for ui, vi in zip(u, v))
cross = lambda u, v: [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
half = sp.Rational(1, 2)
dl = lambda i, j: 1 if i == j else 0
vec = lambda name: list(sp.symbols(f'{name}_1:4'))
NAMES = ['h00', 'hB', 'hC', 'hM', 'hD', 'hE', 'hF', 'hMij', 'a0', 'aB', 'aC', 'aM']


def l1_form(U, p):
    """Most general l = 1 field (both parities) times R^p, from 12 vectors."""
    h = [[0] * 4 for _ in range(4)]
    h[0][0] = dot(U['h00'], n) * R**p
    w = cross(U['hMij'], n)
    for i in range(3):
        v = dot(U['hB'], n) * n[i] + U['hC'][i] + cross(U['hM'], n)[i]
        h[0][i+1] = h[i+1][0] = v * R**p
        for j in range(3):
            h[i+1][j+1] = (dot(U['hD'], n) * n[i] * n[j] + half * (U['hE'][i] * n[j] + U['hE'][j] * n[i])
                           + dot(U['hF'], n) * dl(i, j) + half * (n[i] * w[j] + n[j] * w[i])) * R**p
    A = [dot(U['a0'], n) * R**p] + [(dot(U['aB'], n) * n[i] + U['aC'][i] + cross(U['aM'], n)[i]) * R**p
                                    for i in range(3)]
    return h, A


def unknown_set(tag):
    return {k: vec(f'{k}{tag}') for k in NAMES}


# ------------------------------------------------------------ RN background
def rn_background(order):
    u, w = M / (2 * R), q / (2 * R)
    Bf = (1 + u)**2 - w**2
    s = sp.Symbol('s')
    sub = {M: s * M, q: s * q}
    def ser(e):
        e = sp.series(sp.sympify(e).subs(sub, simultaneous=True), s, 0, order + 1).removeO()
        return {(k, 0, 0): sp.expand(sp.expand(e).coeff(s, k)) for k in range(order + 1)
                if sp.expand(e).coeff(s, k) != 0}
    g00 = ser(-((1 - u**2 + w**2) / Bf)**2)
    gss = ser(Bf**2)
    A0 = ser(-q / (R * Bf))
    return g00, gss, A0


def assemble(pert_h, pert_A, order=3):
    g00, gss, A0 = rn_background(order - 1)
    g = T.mat(lambda a, b: T.add(g00 if a == b == 0 else (gss if a == b else {}), pert_h[a][b]))
    A = [T.add(A0 if a == 0 else {}, pert_A[a]) for a in range(4)]
    return g, A


def equations(g, A, Nmax):
    Ein, Mx = T.einstein_maxwell(g, A)
    comps = [Ein[a][b] for a in range(4) for b in range(a, 4)] + Mx
    eqs = []
    for N in range(0, Nmax + 1):
        for cpt in comps:
            ex = cpt.get((N, 1, 0), 0)
            if ex == 0:
                continue
            for nv in DIRECTIONS:
                e = sp.expand(sp.sympify(ex).subs({x: R * nv[0], y: R * nv[1], z: R * nv[2]},
                                                  simultaneous=True))
                e = sp.expand(e * R**8)
                eqs += [cf for cf in sp.Poly(e, R, t).coeffs() if cf != 0]
    return eqs


def put(dst, slot, h, A):
    for a in range(4):
        A_ = A[a]
        if A_ != 0:
            dst[1][a] = T.add(dst[1][a], {(slot, 1, 0): A_})
        for b in range(4):
            if h[a][b] != 0:
                dst[0][a][b] = T.add(dst[0][a][b], {(slot, 1, 0): h[a][b]})


def empty():
    return [[[{} for _ in range(4)] for _ in range(4)], [{} for _ in range(4)]]


def solve_ordered(eqs, responses, data):
    """Row-reduce with the data columns last.  Returns (solution of the pivot
    response columns in terms of free columns, list of conditions on the data,
    free response columns)."""
    from sympy.polys.matrices import DomainMatrix
    one = sp.Symbol('ONE')                  # inhomogeneous terms ride on a unit column
    cols = list(responses) + list(data) + [one]
    Am, bm = sp.linear_eq_to_matrix(eqs, cols[:-1])
    Am = Am.row_join(-bm)
    dm = DomainMatrix.from_Matrix(Am).to_field()
    rr, piv = dm.rref()
    rr = rr.to_Matrix()
    nresp = len(responses)
    sol, cond = {}, []
    for i, pc in enumerate(piv):
        row = rr.row(i)
        expr = sum(row[j] * cols[j] for j in range(len(cols)) if j != pc)
        if pc < nresp:
            sol[cols[pc]] = sp.factor(-expr)
        else:
            cond.append(sp.factor(cols[pc] + expr))
    free = [c for j, c in enumerate(responses) if j not in piv]
    sol = {k: sp.factor(v.subs(one, 1)) for k, v in sol.items()}
    cond = [sp.factor(c.subs(one, 1)) for c in cond]
    return sol, cond, free


def solve_order2():
    a, E, B = vec('a'), vec('E'), vec('B')
    U1, U2, U3 = unknown_set('g1'), unknown_set('g2'), unknown_set('g3')
    P = empty()
    data_h = [[0] * 4 for _ in range(4)]
    data_h[0][0] = -2 * dot(a, xs)
    data_A = [dot(E, xs)] + [half * cross(B, xs)[i] for i in range(3)]
    put(P, 1, data_h, data_A)
    for nm in ('h00', 'a0', 'aM'):          # these r^1 structures are the data
        U1[nm] = [0, 0, 0]
    for U, p, k in [(U1, 1, 1), (U2, 0, 2), (U3, -1, 3)]:
        h, A = l1_form(U, p)
        put(P, k, h, A)
    g, A = assemble(P[0], P[1])
    t0 = time.time()
    eqs = equations(g, A, 3)
    print(f'order 2: {len(eqs)} equations, {time.time()-t0:.0f}s', flush=True)
    unk = [s for U in (U3, U2, U1) for v in U.values() for s in v if s != 0]
    sol, cond, free = solve_ordered(eqs, unk, a + E + B)
    print('order-2 conditions on the data:', cond)
    print('order-2 free response parameters:', free)
    # impose the conditions by solving them for a
    ca = sp.solve(cond, a, dict=True)[0] if cond else {}
    sol = {k: sp.factor(v.subs(ca)) for k, v in sol.items()}
    sol.update(ca)
    return sol, (a, E, B, U1, U2, U3), P


if __name__ == '__main__':
    T.set_limits(3, 1, 0)
    sol, syms, P = solve_order2()
    a, E, B, U1, U2, U3 = syms
    print('a =', [sp.factor(sol.get(s, s)) for s in a])
    for tag, U in [('r^1', U1), ('r^0', U2), ('r^-1', U3)]:
        for k, v in U.items():
            vals = [sp.factor(sp.sympify(sol.get(s, s))) for s in v]
            if any(vv != 0 for vv in vals):
                print(tag, k, vals)
    pickle.dump({str(k): v for k, v in sol.items()}, open('ckpt15_order2.pkl', 'wb'))
