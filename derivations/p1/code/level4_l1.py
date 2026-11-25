"""Level 4, near-zone first order: RN (to O(M, q)) in an accelerated frame with a
uniform background field. Solve for the O(eps) r^0 l = 1 response."""
import sympy as sp
import trunc as T
from trunc import t, x, y, z, R

M, q, eps = sp.symbols('M q epsilon')
xs = (x, y, z)
n = [xi / R for xi in xs]
a = sp.symbols('a_x a_y a_z')
E = sp.symbols('E_x E_y E_z')
B = sp.symbols('B_x B_y B_z')
dot = lambda u, v: sum(ui * vi for ui, vi in zip(u, v))
cross = lambda u, v: [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
half = sp.Rational(1, 2)
dl = lambda i, j: 1 if i == j else 0

# general l = 1 r^0 forms with unknown vectors
U = {k: sp.symbols(f'{k}_1:4') for k in
     ['h00', 'hB', 'hC', 'hM', 'hD', 'hE', 'hF', 'hMij', 'a0', 'aB', 'aC', 'aM']}
unknowns = [s for v in U.values() for s in v]


def response():
    h = [[0] * 4 for _ in range(4)]
    h[0][0] = dot(U['h00'], n)
    for i in range(3):
        v = dot(U['hB'], n) * n[i] + U['hC'][i] + cross(U['hM'], n)[i]
        h[0][i+1] = h[i+1][0] = v
        for j in range(3):
            w = cross(U['hMij'], n)
            h[i+1][j+1] = (dot(U['hD'], n) * n[i] * n[j] + half * (U['hE'][i] * n[j] + U['hE'][j] * n[i])
                           + dot(U['hF'], n) * dl(i, j) + half * (n[i] * w[j] + n[j] * w[i]))
    A = [dot(U['a0'], n)] + [dot(U['aB'], n) * n[i] + U['aC'][i] + cross(U['aM'], n)[i] for i in range(3)]
    return h, A


def build(h1, A1):
    """g = eta + eps b + mu1 h1 ; A = eps A0 + mu1 A1  (h1, A1 given as {(i,j,k): expr} entries)."""
    g = T.mat(lambda i, j: T.add(T.const(T.eta[i, j]),
                                 {(1, 0, 0): 2 * M / R} if i == j else {},
                                 h1[i][j]))
    A = [T.add({(1, 0, 0): -q / R} if i == 0 else {}, A1[i]) for i in range(4)]
    return g, A


if __name__ == '__main__':
    import time
    T.set_limits(1, 1, 0)
    hr, Ar = response()
    h1 = T.mat(lambda i, j: T.S({(0, 1, 0): (-2 * dot(a, xs) if i == j == 0 else 0),
                                 (1, 1, 0): hr[i][j]}))
    A0 = [dot(E, xs)] + [half * cross(B, xs)[i] for i in range(3)]
    A1 = [T.S({(0, 1, 0): A0[i], (1, 1, 0): Ar[i]}) for i in range(4)]
    g, A = build(h1, A1)
    t0 = time.time()
    Ein, Mx = T.einstein_maxwell(g, A)
    print('computed', time.time() - t0)
    comps = [Ein[i][j].get((1, 1, 0), 0) for i in range(4) for j in range(i, 4)] + [Mx[i].get((1, 1, 0), 0) for i in range(4)]
    assert all(Ein[i][j].get((0, 1, 0), 0) == 0 for i in range(4) for j in range(4))
    from series import DIRECTIONS
    eqs = []
    for c in comps:
        for nv in DIRECTIONS:
            e = sp.expand(c.subs({x: R * nv[0], y: R * nv[1], z: R * nv[2]}, simultaneous=True))
            e = sp.expand(e * R**2)
            assert not e.has(R), e
            eqs.append(e)
    sol = sp.solve(eqs, unknowns + list(a), dict=True)
    print(len(sol))
    sol = sol[0]
    print('a =', [sol[s] for s in a])
    for k, v in U.items():
        print(k, [sp.factor(sol.get(s, s)) for s in v])
    free = [s for s in unknowns if s not in sol]
    print('free:', free)
    # which constraints on a, E, B ?
    res = [sp.factor(sp.expand(e.subs(sol))) for e in eqs]
    print('residual constraints:', set(r for r in res if r != 0))
