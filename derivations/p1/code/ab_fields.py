"""Step 2: the Alekseev-Belinski scalars H = -g(xi,xi) and Phi = A.xi expanded in the small
black hole (eps) at a fixed far-zone point, then near the small hole (r -> 0).
Units M2 = 1; small source at z1 = 0, big source at z2 = l.  Point: rho = r sin(th), z = r cos(th).
Exact functions (AB eqs. (3), (5), (6), (8)):
  rho = sqrt((r_k - m_k)^2 - s_k^2) sin th_k,   z = z_k + (r_k - m_k) cos th_k ;
  H = [(r1-m1)^2 - s1^2 + g^2 sin^2 th2][(r2-m2)^2 - s2^2 + g^2 sin^2 th1] / D^2 ,
  Phi = [(e1-g)(r2-m2) + (e2+g)(r1-m1) + g(m1 cos th1 + m2 cos th2)] / D ,
  D = r1 r2 - (e1 - g - g cos th2)(e2 + g - g cos th1).
For a point off the rods: (r_k - m_k) = (P_k+ + P_k-)/2, (r_k - m_k) cos th_k = z - z_k,
P_k+- = sqrt(rho^2 + (z - z_k +- s_k)^2).
Output: coefficients of eps^0, eps^1, eps^2 of H and Phi as Laurent series in r about the
small hole (to O(r^1)), with their angular dependence.
"""
import sys, sympy as sp
eps, r, c = sp.symbols('epsilon r c', real=True)          # c = cos(th) about the small hole
Q, l, e = sp.symbols('Q ell e', positive=True)
a = sp.symbols('a1:4')                                    # m1 = eps a1 + eps^2 a2 + eps^3 a3
NE = 3                                                    # eps^0..eps^2
m1 = sum(a[k] * eps**(k + 1) for k in range(NE))
e1 = eps * e
M2, Q2 = 1, Q
g = (M2 * e1 - m1 * Q2) / (l + m1 + M2)
s1sq = m1**2 - e1**2 + 2 * e1 * g
s2sq = M2**2 - Q2**2 - 2 * Q2 * g
rho2 = r**2 * (1 - c**2)
zz = r * c


def ser(x, n=NE):
    return sp.expand(sp.series(x, eps, 0, n).removeO())


def bipolar(zk, sksq):
    """(r_k - m_k) and cos th_k at the point, as eps series (sksq may be O(eps))."""
    # P+ P- product and sum via x = r_k - m_k : x^2 - sk^2 cos^2 = ... use the exact
    # relations: x^2 = (rho^2 + (z-zk)^2 + sk^2 + sqrt((rho^2+(z-zk)^2+sk^2)^2 - 4 sk^2 (z-zk)^2))/2
    w2 = rho2 + (zz - zk)**2
    x2 = (w2 + sksq + sp.sqrt((w2 + sksq)**2 - 4 * sksq * (zz - zk)**2)) / 2
    x = sp.sqrt(x2)
    return x, (zz - zk) / x


x1, c1 = bipolar(0, s1sq)
x2, c2 = bipolar(l, s2sq)
r1, r2 = x1 + m1, x2 + M2
sin1sq, sin2sq = 1 - c1**2, 1 - c2**2
D = r1 * r2 - (e1 - g - g * c2) * (Q2 + g - g * c1)
H = ((x1**2 - s1sq + g**2 * sin2sq) * (x2**2 - s2sq + g**2 * sin1sq)) / D**2
Phi = ((e1 - g) * (r2 - M2) + (Q2 + g) * (r1 - m1) + g * (m1 * c1 + M2 * c2)) / D
if __name__ == '__main__':
    out = {}
    for name, F in [('H', H), ('Phi', Phi)]:
        Fs = ser(F)
        for k in range(NE):
            ck = Fs.coeff(eps, k)
            # Laurent expansion in r about the small hole, to O(r)
            lo = -2 * k - 1
            sr = sp.expand(sp.series(ck, r, 0, 2).removeO())
            terms = {p: sp.factor(sp.simplify(sr.coeff(r, p))) for p in range(lo, 2)}
            out[(name, k)] = terms
            for p, v in terms.items():
                if v != 0:
                    print(f'{name} eps^{k} r^{p}:', v, flush=True)
    import pickle
    pickle.dump(out, open('ab_fields.pkl', 'wb'))
