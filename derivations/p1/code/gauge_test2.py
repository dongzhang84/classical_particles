"""Level 2 gauge checks on the matrix written by level2.py.

Pure-gauge perturbations must lie in the null space of the full matrix
(Einstein rows A-F with the EM source, Maxwell rows X, Y, Z).
For a diffeomorphism xi, the potential changes by L_xi A with A = -(q/r) dt;
the metric part is the same as in gauge_test.py (level 1).
For an EM gauge transformation psi, only the potential changes.
"""
import ast
import sympy as sp

M, q = sp.symbols('M q', positive=True)
c = sp.Symbol('c')
lines = open('level2_matrix.txt').read().strip().split('\n')
cols = ast.literal_eval(lines[0])
rows = []
for ln in lines[1:]:
    name, rest = ln.split(' ', 1)
    rows.append((name, list(sp.sympify(rest, locals={'M': M, 'q': q}))))
Mt = sp.Matrix([r for _, r in rows])

def vec(**kw):
    return sp.Matrix([kw.get(k, 0) for k in cols])

tests = {
    # diffeos: metric content from gauge_test.py, EM content from L_xi(-(q/r)dt)
    'xi_0 = t (c.x)':   vec(DA=2*c, DC=c, KB=-2*M*c, HP=q*c, KCp=q*c),
    'xi_i = (c.x) x_i': vec(DE=2*c, DF=2*c, HA=-2*M*c, HD=4*M*c, HF=-2*M*c, HP=q*c),
    'xi_i = r^2 c_i':   vec(DE=4*c, HA=-2*M*c, HE=4*M*c, HF=-2*M*c, HP=q*c),
    # EM gauge psi = t (c.x): a_0 = c.x, a_i = t c_i
    'psi = t (c.x)':    vec(E0=c, EC=c),
}
ok = True
for name, v in tests.items():
    res = sp.simplify(Mt * v)
    print(f'{name:18s}', [r for r in res.T], '' if all(x == 0 for x in res) else '  <-- NONZERO')
    ok &= all(x == 0 for x in res)
print('all gauge checks pass' if ok else 'SOME GAUGE CHECKS FAIL')
