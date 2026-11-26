"""B1, B2 of the Alekseev-Belinski invariant balance (see calc.tex), analytic.
Units M2 = 1; e = 1 (B_k is proportional to e^k: restore by B_k -> B_k e^k).
Symbols: Q (= e2 > 1 for the naked big source), L (= ell).  Truncated power series in eps
with every coefficient simplified on the fly."""
import sys, sympy as sp
eps = sp.Symbol('epsilon')
v, d = sp.symbols('v d', positive=True)
w = v + d                                          # w > v <=> beyond the zero-force radius
uu = 1 / (1 + v)
Q = sp.factor((1 + uu**2) / (2 * uu))             # Q > 1 (naked big source)
pp = sp.factor((1 - uu**2) / (2 * uu))            # sqrt(Q^2 - 1)
tt = pp / (1 + w)
L = sp.factor((pp**2 / tt - tt) / 2)              # L > 0
Wsq = sp.factor(L**2 + Q**2 - 1)                  # = ((p^2/t + t)/2)^2
K = 3                                   # B to eps^(K-1)
NN = K + 1


def simp(x):
    return sp.factor(sp.cancel(sp.sympify(x)))


def psqrt(x):
    x = sp.factor(sp.cancel(x))
    r = sp.sqrt(x)
    r = sp.factor(sp.powsimp(sp.powdenest(r, force=True), force=True))
    assert not r.has(sp.Abs), r
    return r


class S:
    """power series in eps: dict k -> coefficient (k >= lo)."""
    def __init__(self, c):
        self.c = {k: simp(v) for k, v in c.items() if simp(v) != 0}

    @staticmethod
    def of(x):
        return S({0: x})

    def __add__(a, b):
        b = b if isinstance(b, S) else S.of(b)
        ks = set(a.c) | set(b.c)
        return S({k: a.c.get(k, 0) + b.c.get(k, 0) for k in ks if k < NN})

    __radd__ = __add__

    def __neg__(a):
        return S({k: -v for k, v in a.c.items()})

    def __sub__(a, b):
        return a + (-(b if isinstance(b, S) else S.of(b)))

    def __rsub__(a, b):
        return S.of(b) - a

    def __mul__(a, b):
        b = b if isinstance(b, S) else S.of(b)
        out = {}
        for i, x in a.c.items():
            for j, y in b.c.items():
                if i + j < NN:
                    out[i + j] = out.get(i + j, 0) + x * y
        return S(out)

    __rmul__ = __mul__

    def lo(a):
        return min(a.c)

    def shift(a, n):
        return S({k + n: v for k, v in a.c.items()})

    def inv(a):
        l0 = a.lo(); a0 = a.c[l0]
        d = S({k - l0: v / a0 for k, v in a.c.items() if k != l0})
        out, term = S.of(1), S.of(1)
        for _ in range(NN):
            term = -(term * d)
            out = out + term
        return (out * (1 / a0)).shift(-l0)

    def __truediv__(a, b):
        b = b if isinstance(b, S) else S.of(b)
        return a * b.inv()

    def sqrt(a):
        l0 = a.lo(); assert l0 % 2 == 0
        a0 = a.c[l0]
        d = S({k - l0: v / a0 for k, v in a.c.items() if k != l0})
        out, term = S.of(1), S.of(1)
        for n in range(1, NN):
            term = term * d
            out = out + sp.binomial(sp.Rational(1, 2), n) * term
        return (out * psqrt(a0)).shift(l0 // 2)


M2, e = 1, 1
# m1 solved order by order (from ab_closed.py output / re-derived here)
a = sp.symbols('a1:6')
m1 = S({k + 1: a[k] for k in range(NN - 1)})
e1 = S({1: e})
gam = (M2 * e1 - m1 * Q) / (L + M2 + m1)
F = (e1 - gam) * (Q + gam) - m1 * M2
sol = {}
for k in range(1, NN):
    eq = sp.cancel(F.c.get(k, 0).subs(sol))
    if eq == 0:
        continue
    sol[a[k - 1]] = simp(sp.solve(sp.numer(sp.together(eq)), a[k - 1])[0])
    print(f'm1 eps^{k}:', sol[a[k - 1]], flush=True)
sub = lambda s_: S({k: v.subs(sol) for k, v in s_.c.items()})
m1 = sub(m1); gam = sub(gam)
s1 = (m1 * m1 - e1 * e1 + 2 * e1 * gam).sqrt()
s2sq = M2**2 - Q**2 - 2 * Q * gam
D = (m1 + s1) * (M2 + L) - e1 * (Q + gam)
rp = D / (L**2 - s2sq + gam * gam).sqrt()
mu = (rp * rp + e1 * e1) / (2 * rp)
N = 2 * s1 * rp / (rp * rp - e1 * e1)
print('N leading:', N.c.get(0), flush=True)
Mt, Qt = m1 + M2, e1 + Q
one_m = 1 - N * N
R = (Mt + (Mt * Mt - one_m * Qt * Qt).sqrt()) / one_m
B = Qt * R * N * e1 / (mu * (Mt * R - Qt * Qt))
for k in range(K):
    print(f'B eps^{k}:', simp(B.c.get(k, 0)), flush=True)
mu1 = mu.c[1]
print('mu eps^1:', mu1)
print('Q =', Q, '  L =', L)
