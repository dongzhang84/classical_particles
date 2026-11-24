"""Shared tools for the first-order motion calculation (levels 1-4).

Conventions follow Gralla & Wald 2008 (GW), Sec. VI: near-zone Cartesian
coordinates (t, x, y, z) = (tbar, xbar^i), signature (-,+,+,+), n^i = x^i/r.

Two jobs:
  * linearized Ricci tensor of a perturbation about a background that is
    flat plus a term linear in a small parameter (M), kept to first order in M;
  * extraction of the coefficient of t^k r^-p, and its l = 1 tensor-harmonic
    pieces, by evaluating on rational unit vectors and solving a linear fit.
"""
import itertools
import sympy as sp

t, x, y, z = sp.symbols('t x y z', real=True)
X = (t, x, y, z)
r = sp.sqrt(x**2 + y**2 + z**2)
eta = sp.diag(-1, 1, 1, 1)
eps = sp.Symbol('epsilon')


def christoffel(g, ginv):
    G = [[[0] * 4 for _ in range(4)] for _ in range(4)]
    dg = [[[sp.diff(g[a, b], X[c]) for c in range(4)] for b in range(4)] for a in range(4)]
    for a, b, c in itertools.product(range(4), repeat=3):
        G[a][b][c] = sum(ginv[a, d] * (dg[d][b][c] + dg[d][c][b] - dg[b][c][d])
                         for d in range(4)) / 2
    return G


def ricci(G):
    R = sp.zeros(4, 4)
    for b, c in itertools.product(range(4), repeat=2):
        R[b, c] = (sum(sp.diff(G[a][b][c], X[a]) for a in range(4))
                   - sum(sp.diff(G[a][b][a], X[c]) for a in range(4))
                   + sum(G[a][a][d] * G[d][b][c] for a in range(4) for d in range(4))
                   - sum(G[a][c][d] * G[d][b][a] for a in range(4) for d in range(4)))
    return R


def truncate(expr, keep):
    """Keep only the monomials eps^1 * lam^k for the (symbol, power) pairs in keep."""
    expr = sp.expand(expr)
    out = 0
    for sym_powers in keep:
        term = expr
        for s, p in sym_powers:
            term = sp.expand(term).coeff(s, p)
        mono = sp.Mul(*[s**p for s, p in sym_powers])
        out += mono * term
    return out


def linearized_ricci(b1, h, M):
    """Ricci of eta + M*b1 + eps*h, coefficient of eps^1, kept through M^1.

    Returns (R0, R1): the eps^1 M^0 and eps^1 M^1 parts (R1 without the M).
    """
    g = eta + M * b1 + eps * h
    # inverse to the needed order: second order only in the cross term M*eps
    ei = eta  # eta is its own inverse
    ginv = (ei - M * ei * b1 * ei - eps * ei * h * ei
            + M * eps * (ei * b1 * ei * h * ei + ei * h * ei * b1 * ei))
    G = christoffel(g, ginv)
    # cut each Christoffel to monomials 1, M, eps, M*eps before squaring
    def cut(e):
        e = sp.expand(e)
        return sum(sp.expand(e).coeff(eps, i).coeff(M, j) * eps**i * M**j
                   for i in (0, 1) for j in (0, 1))
    G = [[[cut(G[a][b][c]) for c in range(4)] for b in range(4)] for a in range(4)]
    R = ricci(G)
    R0 = sp.zeros(4, 4)
    R1 = sp.zeros(4, 4)
    for a, b in itertools.product(range(4), repeat=2):
        e = sp.expand(R[a, b]).coeff(eps, 1)
        R0[a, b] = sp.expand(e).coeff(M, 0)
        R1[a, b] = sp.expand(e).coeff(M, 1)
    return R0, R1


# Rational unit vectors (Pythagorean quadruples), so r is rational.
DIRECTIONS = [
    (sp.Rational(1, 3), sp.Rational(2, 3), sp.Rational(2, 3)),
    (sp.Rational(2, 7), sp.Rational(3, 7), sp.Rational(6, 7)),
    (sp.Rational(3, 13), sp.Rational(4, 13), sp.Rational(12, 13)),
    (sp.Rational(2, 3), sp.Rational(-1, 3), sp.Rational(2, 3)),
    (sp.Rational(-6, 7), sp.Rational(2, 7), sp.Rational(3, 7)),
    (sp.Rational(4, 9), sp.Rational(4, 9), sp.Rational(-7, 9)),
    (sp.Rational(1, 9), sp.Rational(-4, 9), sp.Rational(8, 9)),
    (sp.Rational(-2, 11), sp.Rational(6, 11), sp.Rational(9, 11)),
]

rho = sp.Symbol('rho', positive=True)


def radial_coeff(expr, n, tpow, rpow):
    """Coefficient of t^tpow * rho^rpow after x^i = rho * n^i (n rational unit)."""
    e = expr.subs({x: rho * n[0], y: rho * n[1], z: rho * n[2]}, simultaneous=True)
    e = sp.expand(sp.simplify(e))
    e = sp.expand(e).coeff(t, tpow)
    return _laurent_coeff(e, rpow)


def _laurent_coeff(e, rpow):
    e = sp.expand(e)
    out = 0
    for term in sp.Add.make_args(e):
        c, p = term.as_coeff_exponent(rho)
        if p == rpow:
            out += c
    return sp.expand(out)


def fit_l1(values, ansatz, unknowns):
    """Solve for unknowns such that ansatz(n) == values[n] on all DIRECTIONS.

    values: dict n -> expression; ansatz: function n -> expression linear in
    unknowns. Returns (solution dict, list of residuals that must vanish).
    """
    eqs = [sp.expand(ansatz(n) - values[n]) for n in values]
    sol = sp.solve(eqs, unknowns, dict=True)
    if not sol:
        raise ValueError('l=1 ansatz does not fit: output has other harmonics')
    sol = sol[0]
    resid = [sp.simplify(e.subs(sol)) for e in eqs]
    return sol, resid
