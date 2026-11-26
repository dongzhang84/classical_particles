"""Checkpoint 2, part 1: static, axisymmetric, purely electric electrovac
background.  Near zone through third order, fully nonlinear, in P gauge.

Fields are T-even: h_00, h_ij, A_0 (h_0i = A_i = 0).  Axisymmetry about z:
    scalar   S(R, z) ,
    h_ij  =  delta_ij S1 + x_i x_j S2 + (x_i zhat_j + zhat_i x_j) S3 + zhat_i zhat_j S4 .
Near-zone order m, slot k (power of M, q, including the far-zone data): every
term has length power m - k.  Slot 0 = background (smooth polynomials only).

Order 1 data: h_00 = -2 a z, A_0 = E z (uniform field and acceleration).
All free parameters are kept as symbols (data or gauge); the mass dipole and
charge dipole (l = 1 decaying modes) are set to zero (P gauge).

Output: the l = 1 (z) solvability condition at each order, compared with the
test-body condition  M d_z sqrt(-g_00) = q d_z A_0  at x = 0, evaluated with the
smooth parts of the fields (far-zone order = slot).
"""
import pickle, sys, time
import sympy as sp
import trunc as T
from trunc import t, x, y, z, R

M, q = sp.symbols('M q', positive=True)
if '--q0' in sys.argv:
    q = sp.Integer(0)
xs = (x, y, z)
zh = (0, 0, 1)
LMAX = 6
LMAX_ORDER = {1: 4, 2: 6, 3: 8}
import os
EPS = int(os.environ.get('CK2_EPS', '3'))
USELOG = '--log' in sys.argv
RADIAL = '--radial' in sys.argv
ISO = '--iso' in sys.argv
ISOALL = '--isoall' in sys.argv   # isotropic (conformally flat) spatial metric in ALL slots
CONF = '--conf' in sys.argv     # conformastatic background coordinates: order-1 h_ij = -2 a z delta_ij
KEEPALL = '--keepall' in sys.argv
TIME = '--time' in sys.argv
DIPOLE = '--dipole' in sys.argv
LORENZ = '--lorenz' in sys.argv
LIGHT = '--light' in sys.argv
d1free = []


def rn_background(order):
    u, w = M / (2 * R), q / (2 * R)
    Bf = (1 + u)**2 - w**2
    s = sp.Symbol('s')
    sub = {M: s * M, q: s * q}
    def ser(e):
        e = sp.expand(sp.series(sp.sympify(e).subs(sub, simultaneous=True), s, 0, order + 1).removeO())
        return {(k, 0, 0): sp.expand(e.coeff(s, k)) for k in range(order + 1) if e.coeff(s, k) != 0}
    return ser(-((1 - u**2 + w**2) / Bf)**2), ser(Bf**2), ser(-q / (R * Bf))


def scalar(name, P, smooth, lmax=None):
    """General axisymmetric scalar of length power P; returns (expr, symbols)."""
    syms, e = [], 0
    for a_ in range(0, (lmax or LMAX) + 1):
        if smooth and (P - a_ < 0 or (P - a_) % 2):
            continue
        c = sp.Symbol(f'{name}_{a_}')
        syms.append(c)
        e += c * z**a_ * R**(P - a_)
    return e, syms


LOG = sp.log(R)


def ansatz(m, k, withlog=False):
    """General static axisymmetric T-even field at order m, slot k
    (times log R if withlog)."""
    P = m - k
    sm = (k == 0)
    tag = f'o{m}s{k}' + ('L' if withlog else '')
    h00, s0 = scalar(f'h{tag}', P, sm, LMAX_ORDER.get(m))
    a0, s1 = scalar(f'A{tag}', P, sm, LMAX_ORDER.get(m))
    S1, s2 = scalar(f'S1{tag}', P, sm, LMAX_ORDER.get(m))
    S2, s3 = scalar(f'S2{tag}', P - 2, sm, LMAX_ORDER.get(m))
    S3, s4 = scalar(f'S3{tag}', P - 1, sm, LMAX_ORDER.get(m))
    S4, s5 = scalar(f'S4{tag}', P, sm, LMAX_ORDER.get(m))
    h = [[0] * 4 for _ in range(4)]
    h[0][0] = h00
    for i in range(3):
        for j in range(3):
            h[i+1][j+1] = ((1 if i == j else 0) * S1 + xs[i] * xs[j] * S2
                           + (xs[i] * zh[j] + zh[i] * xs[j]) * S3 + zh[i] * zh[j] * S4)
    A = [a0, 0, 0, 0]
    if withlog:
        h = [[e * LOG for e in row] for row in h]
        A = [e * LOG for e in A]
    return h, A, s0 + s1 + s2 + s3 + s4 + s5


def ansatz_t(m, k, d):
    """Time-dependent pieces tbar^d at order m, slot k (tbar counts as a length).
    d even: T-even fields h00, h_ij, A0; d odd: T-odd fields h_0i, A_i with
    the axisymmetric polar form V_i = x_i f + zhat_i g."""
    P = m - k
    lm = 6
    tag = f'o{m}s{k}t{d}'
    h = [[0] * 4 for _ in range(4)]
    A = [0, 0, 0, 0]
    syms = []
    if d % 2 == 0:
        # rebuild with shifted power and new names
        h00, s0 = scalar(f'h{tag}', P - d, False, lm)
        a0, s1 = scalar(f'A{tag}', P - d, False, lm)
        S1, s2 = scalar(f'S1{tag}', P - d, False, lm)
        S2, s3 = scalar(f'S2{tag}', P - d - 2, False, lm)
        S3, s4 = scalar(f'S3{tag}', P - d - 1, False, lm)
        S4, s5 = scalar(f'S4{tag}', P - d, False, lm)
        h[0][0] = h00 * t**d
        for i in range(3):
            for j in range(3):
                h[i+1][j+1] = t**d * ((1 if i == j else 0) * S1 + xs[i] * xs[j] * S2
                                      + (xs[i] * zh[j] + zh[i] * xs[j]) * S3 + zh[i] * zh[j] * S4)
        A[0] = a0 * t**d
        syms = s0 + s1 + s2 + s3 + s4 + s5
    else:
        f1, s0 = scalar(f'Vf{tag}', P - d - 1, False, lm)
        g1, s1 = scalar(f'Vg{tag}', P - d, False, lm)
        f2, s2 = scalar(f'Wf{tag}', P - d - 1, False, lm)
        g2, s3 = scalar(f'Wg{tag}', P - d, False, lm)
        for i in range(3):
            h[0][i+1] = h[i+1][0] = t**d * (xs[i] * f1 + zh[i] * g1)
            A[i+1] = t**d * (xs[i] * f2 + zh[i] * g2)
        syms = s0 + s1 + s2 + s3
    return h, A, syms


class Field:
    def __init__(self):
        self.h = [[{} for _ in range(4)] for _ in range(4)]
        self.A = [{} for _ in range(4)]

    def add(self, key, h, A):
        for a_ in range(4):
            if A[a_] != 0:
                self.A[a_] = T.add(self.A[a_], {key: A[a_]})
            for b in range(4):
                if h[a_][b] != 0:
                    self.h[a_][b] = T.add(self.h[a_][b], {key: h[a_][b]})

    def subs(self, sol):
        f = Field()
        f.h = [[{k: sp.expand(sp.sympify(v).subs(sol)) for k, v in c.items()} for c in row] for row in self.h]
        f.A = [{k: sp.expand(sp.sympify(v).subs(sol)) for k, v in c.items()} for c in self.A]
        f.h = [[{k: v for k, v in c.items() if v != 0} for c in row] for row in f.h]
        f.A = [{k: v for k, v in c.items() if v != 0} for c in f.A]
        return f


def reduce_eq(e):
    """Coefficients of an equation component after y -> 0, x^2 -> R^2 - z^2."""
    Ls = sp.Symbol('LOGR')
    e = sp.expand(sp.sympify(e).subs(y, 0).subs(sp.log(R), Ls))
    if e == 0:
        return []
    if e == 0:
        return []
    lo = max([-term.as_powers_dict().get(R, 0) for term in sp.Add.make_args(e)] + [0])
    num = sp.expand(e * R**lo)
    assert not num.has(1 / R)
    num = sp.expand(sp.rem(sp.Poly(num, x), sp.Poly(x**2 + z**2 - R**2, x)).as_expr())
    return [c for c in sp.Poly(num, x, z, R, Ls, t).coeffs() if c != 0]


def equations(field, m, tagged=False):
    g00, gss, A0 = rn_background(EPS)
    G = T.mat(lambda a_, b: T.add(g00 if a_ == b == 0 else (gss if a_ == b else {}), field.h[a_][b]))
    AA = [T.add(A0 if a_ == 0 else {}, field.A[a_]) for a_ in range(4)]
    Ein, Mx = T.einstein_maxwell(G, AA)
    eqs, tags = [], []
    names = [f'E{a_}{b}' for a_ in range(4) for b in range(a_, 4)] + [f'Mx{a_}' for a_ in range(4)]
    for nm, cpt in zip(names, [Ein[a_][b] for a_ in range(4) for b in range(a_, 4)] + Mx):
        for N in range(EPS + 1):
            ex = cpt.get((N, m, 0), 0)
            if ex != 0:
                r = reduce_eq(ex)
                eqs += r
                tags += [(nm, N)] * len(r)
    return (eqs, tags) if tagged else eqs


def slot_of_g(s):
    return int(str(s).rsplit('_', 1)[0].split('s')[-1].split('t')[0].rstrip('L'))


def lorenz_eqs(field, m):
    """Lorenz conditions on the far-zone FIRST-order fields (slot 1) about the
    far-zone background (slot 0), at near-zone order m:
        g^{ac} nabla_c (h_ab - g_ab h/2) = 0,   g^{ab} nabla_a A_b = 0 ."""
    g00, gss, A0 = rn_background(1)
    # background: eta + slot-0 parts ; first-order field: slot-1 parts (incl. RN eps^1)
    keep = lambda dct, k: {kk: v for kk, v in dct.items() if kk[0] == k}
    gb = T.mat(lambda a_, b: T.add(T.const(T.eta[a_, b]), keep(field.h[a_][b], 0)))
    h1 = T.mat(lambda a_, b: T.add(keep(g00 if a_ == b == 0 else (gss if a_ == b else {}), 1),
                                   keep(field.h[a_][b], 1)))
    A1 = [T.add(keep(A0, 1) if a_ == 0 else {}, keep(field.A[a_], 1)) for a_ in range(4)]
    gi = T.inverse(gb)
    G = T.christoffel(gb, gi)
    tr = T.add(*[T.mul(gi[c][d], h1[c][d]) for c in range(4) for d in range(4)])
    hb = T.mat(lambda a_, b: T.add(h1[a_][b], T.scal(-sp.Rational(1, 2), T.mul(gb[a_][b], tr))))
    out = []
    for b in range(4):
        terms = []
        for a_ in range(4):
            for c in range(4):
                cov = T.add(T.d(hb[a_][b], c),
                            T.scal(-1, T.add(*[T.mul(G[d][c][a_], hb[d][b]) for d in range(4)])),
                            T.scal(-1, T.add(*[T.mul(G[d][c][b], hb[a_][d]) for d in range(4)])))
                terms.append(T.mul(gi[a_][c], cov))
        out.append(T.add(*terms).get((1, m, 0), 0))
    div = []
    for a_ in range(4):
        for b in range(4):
            div.append(T.mul(gi[a_][b], T.add(T.d(A1[b], a_),
                                              T.scal(-1, T.add(*[T.mul(G[c][a_][b], A1[c]) for c in range(4)])))))
    out.append(T.add(*div).get((1, m, 0), 0))
    eqs = []
    for e_ in out:
        if e_ != 0:
            eqs += reduce_eq(e_)
    return eqs


def solve_ordered(eqs, responses, data):
    from sympy.polys.matrices import DomainMatrix
    one = sp.Symbol('ONE')
    cols = list(responses) + list(data)
    Am, bm = sp.linear_eq_to_matrix(eqs, cols)
    Am = Am.row_join(-bm)
    cols = cols + [one]
    gens = sorted(Am.free_symbols, key=str)
    if gens:
        dom = sp.QQ.frac_field(*gens)
        dm = DomainMatrix([[dom.from_sympy(Am[i, j]) for j in range(Am.cols)] for i in range(Am.rows)],
                          Am.shape, dom)
    else:
        dm = DomainMatrix.from_Matrix(Am).to_field()
    print('   rref domain:', dm.domain, flush=True)
    rr, piv = dm.rref()
    rr = rr.to_Matrix()
    nresp = len(responses)
    sol, cond = {}, []
    for i, pc in enumerate(piv):
        row = rr.row(i)
        expr = sum(row[j] * cols[j] for j in range(len(cols)) if j != pc)
        if pc < nresp:
            sol[cols[pc]] = sp.factor(sp.sympify(-expr).subs(one, 1))
        else:
            cond.append(sp.factor(sp.sympify(cols[pc] + expr).subs(one, 1)))
    free = [c for j, c in enumerate(responses) if j not in piv]
    return sol, cond, free


if __name__ == '__main__':
    T.set_limits(EPS, 3, 0)
    a, E = sp.symbols('a E')
    if '--vacuum' in sys.argv:
        E = sp.Integer(0)
    F = Field()
    hd = [[0] * 4 for _ in range(4)]; hd[0][0] = -2 * a * z
    if CONF:
        for i_ in range(1, 4):
            hd[i_][i_] = -2 * a * z
    F.add((0, 1, 0), hd, [E * z, 0, 0, 0])
    log = {}
    data_syms = [a, E]
    orders = (1, 2, 3)
    if '--resume' in sys.argv:
        dd = pickle.load(open('ckpt2_axi_o2.pkl', 'rb'))
        dd['log'] = {k: v for k, v in dd['log'].items() if k < 3}
        F, log, orders = dd['F'], dd['log'], (3,)
    for m in orders:
        unk_resp, unk_data = [], []
        gauge_eqs = []
        for k in range(0, EPS + 1):
            if m == 1 and k == 0:
                continue                      # order-1 background = the data above
            h, A, syms = ansatz(m, k)
            if RADIAL and k == 0:
                gauge_eqs.extend(sum([reduce_eq(sum(xs[j] * h[i+1][j+1] for j in range(3))) for i in range(3)], []))
            if (ISO and k == 0) or ISOALL:      # conformally flat: traceless part of h_ij vanishes
                trh = sum(h[i+1][i+1] for i in range(3))
                for i in range(3):
                    for j in range(i, 3):
                        gauge_eqs.extend(reduce_eq(sp.expand(h[i+1][j+1] - (1 if i == j else 0) * trh / 3)))
            if m == 3 and k >= 1 and TIME:
                for d_ in (1, 2):
                    ht, At, st = ansatz_t(m, k, d_)
                    F.add((k, m, 0), ht, At)
                    unk_resp.extend(st)
            if m == 3 and k >= 1 and USELOG:
                hL, AL, symsL = ansatz(m, k, withlog=True)
                F.add((k, m, 0), hL, AL)
                unk_resp.extend(symsL)
            # P gauge: no decaying l = 1 modes (mass / charge dipole) at order 1, slot 3
            F.add((k, m, 0), h, A)
            (unk_data if k == 0 else unk_resp).extend(syms)
        # order responses: highest slot (most singular) first, so data stay free
        unk_resp.sort(key=lambda s: -int(str(s).rsplit('_', 1)[0].split('s')[-1].split('t')[0].rstrip('L')))
        t0 = time.time()
        T.set_limits(EPS, m, 0)
        eqs, tags = equations(F, m, tagged=True)
        eqs += gauge_eqs; tags += [('radial', -1)] * len(gauge_eqs)
        if LORENZ:
            le = lorenz_eqs(F, m)
            print(f'   Lorenz conditions: {len(le)} equations', flush=True)
            eqs += le; tags += [('lorenz', -1)] * len(le)
        print(f'order {m}: {len(eqs)} equations, {len(unk_resp)+len(unk_data)} unknowns, '
              f'{time.time()-t0:.0f}s', flush=True)
        dcols = data_syms if m == 1 else (log[2]['free'] if m == 3 else [])
        if DIPOLE and m == 3:
            dcols = list(d1free)
        if m == 3 and not DIPOLE:
            a2s, e2s = sp.Symbol('ho3s2_1'), sp.Symbol('Ao3s2_1')
            v1 = [] if (LORENZ and not LIGHT) else [u for u in unk_resp if slot_of_g(u) == 1]
            unk_resp = [u for u in unk_resp if u not in (a2s, e2s) and u not in v1]
            dcols = list(dcols) + v1 + [e2s, a2s]
            if LIGHT:
                dcols = v1 + [e2s, a2s]
        allsym = set().union(*[sp.sympify(e_).free_symbols for e_ in eqs])
        extra = allsym - set(unk_resp + unk_data + list(dcols)) - {M, q, E}
        print('   symbols outside the column list:', sorted(extra, key=str), flush=True)
        sol, cond, free = solve_ordered(eqs, unk_resp + unk_data, dcols)
        print(f'order {m} conditions on data:', cond, flush=True)
        if m == 1 and (E == 0 or q == 0):
            cond = [a]
        if cond and m == 1:
            assert len(cond) == 1, cond
            csol = {a: sp.solve(cond[0], a)[0]}
            print('   solved as', csol)
            sol = {k: sp.factor(sp.sympify(v).subs(csol)) for k, v in sol.items()}
            sol.update(csol)
        # P gauge and the definition of M, q: drop free modes decaying faster than r^0
        slot_of = lambda s: int(str(s).rsplit('_', 1)[0].split('s')[-1].split('t')[0].rstrip('L'))
        # gauge choice: keep only the smooth constants/linear terms as free; zero the rest
        def smooth_free(s_):
            name, a_ = str(s_).rsplit('_', 1)
            pre = name.split('o')[0]
            off = {'h': 0, 'A': 0, 'S1': 0, 'S4': 0, 'S2': 2, 'S3': 1}[pre]
            mm = int(name.split('o')[1].split('s')[0]); kk = slot_of(s_)
            b = mm - kk - off - int(a_)
            return kk <= mm and b >= 0 and b % 2 == 0
        if m == 1 and DIPOLE:
            zero = {s_: 0 for s_ in free if slot_of(s_) <= 2}   # h^H(gamma) = 0, monopoles 0; keep dipoles
            d1free = [s_ for s_ in free if slot_of(s_) >= 3]
            print('   kept order-1 parameters:', d1free, flush=True)
        elif m == 1:
            zero = {s_: 0 for s_ in free}          # P gauge: h^H(gamma) = 0, no monopole/dipole change
        elif m == 2 and DIPOLE:
            import random
            random.seed(1)
            zero = {s_: sp.Rational(random.randint(-9, 9), random.randint(1, 9)) for s_ in free}
        elif m == 2 and LIGHT:
            # P gauge: j(gamma) = 0, no second-order monopole corrections, residual gauge = 0;
            # keep only the background data (slot 0) and the first-order uniform field f1
            zero = {s_: 0 for s_ in free if slot_of(s_) != 0 and str(s_) != 'Ao2s1_1'}
            print('   order-2 parameters kept:', [s_ for s_ in free if s_ not in zero], flush=True)
        else:
            # keep the order-2 free parameters symbolic, except pure high-degree gauge modes
            zero = {} if KEEPALL else {s_: 0 for s_ in free if int(str(s_).rsplit('_', 1)[1]) >= 4}
        sol = {k: sp.factor(sp.sympify(v).subs(zero)) for k, v in sol.items()}
        sol.update(zero)
        if cond and m > 1 and not DIPOLE:
            zz = {s_: 0 for s_ in free}
            zz.update({s_: 0 for s_ in dcols})
            full = {k: sp.sympify(v).subs(zz) for k, v in sol.items()}
            print('order', m, 'residual equations with all free parameters = 0:')
            for e_, tg in zip(eqs, tags):
                r_ = sp.factor(sp.sympify(e_).subs(full).subs(zz))
                if r_ != 0:
                    print('    ', tg, r_)
            print('order', m, 'conditions:')
            for c_ in cond:
                print('    ', c_)
        F = F.subs(sol)
        newfree = [s for s in free if s not in sol]
        if m == 1:
            newfree = []
        print(f'order {m} free parameters ({len(newfree)}):', newfree, flush=True)
        data_syms += newfree
        log[m] = dict(sol=sol, cond=cond, free=newfree)
        pickle.dump({'log': log, 'F': F}, open('ckpt2_axi.pkl', 'wb'))
        tag_ = ("_eps" + str(EPS) if EPS != 3 else "") + ("_radial" if RADIAL else "") + ("_iso" if ISO else "") + ("_conf" if CONF else "") + ("_isoall" if ISOALL else "")
        pickle.dump({'log': log, 'F': F}, open(f'ckpt2_axi_o{m}{tag_}.pkl', 'wb'))
        if m == 2 and '--stop2' in sys.argv:
            break
