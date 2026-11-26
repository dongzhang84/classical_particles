"""Alekseev-Belinski (PRD 76, 021501, 2007) equilibrium of two RN sources,
expanded in the small source (m1 = eps m, e1 = eps e):
    m1 m2 = (e1 - gamma)(e2 + gamma),  gamma = (m2 e1 - m1 e2)/(l + m1 + m2)   (10), (9)
F = (e1-gamma)(e2+gamma) - m1 m2 = eps F1 + eps^2 F2 + eps^3 F3 + ...
"""
import sympy as sp
m, e, M, Q, l, eps = sp.symbols('m e M Q ell epsilon')
m1, e1 = eps * m, eps * e
gam = (M * e1 - m1 * Q) / (l + m1 + M)
F = (e1 - gam) * (Q + gam) - m1 * M
ser = sp.expand(sp.series(F, eps, 0, 4).removeO())
for k in (1, 2, 3):
    print(f'F{k} =', sp.factor(ser.coeff(eps, k)))
# with R = l + M (test-particle radius) and gamma1 = (Me - mQ)/R
R = sp.Symbol('R')
print('\nin terms of R = l + M:')
for k in (1, 2, 3):
    print(f'F{k} =', sp.factor(ser.coeff(eps, k).subs(l, R - M)))
