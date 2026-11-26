"""Compare the order-2 background (slot 0) of ckpt2_axi --light --radial with the MP local data,
then evaluate X_(iii) for an extremal body (q = M)."""
import pickle, sys, sympy as sp
import __main__
import ckpt2_axi as C
__main__.Field = C.Field
from trunc import x, y, z, R
L = sp.Rational(sys.argv[1]) if len(sys.argv) > 1 else sp.Integer(3)
mp = pickle.load(open(f'ck2iii_mp_L{L}.pkl', 'rb'))
E0, A2c, fit = mp['E'], mp['A0'][1], mp['hij'][0]
S4R = fit[sp.Symbol('s4R')]
d = pickle.load(open('ckpt2_axi_o2_radial.pkl', 'rb'))
F = d['F']
sub = {sp.Symbol('E'): E0, sp.Symbol('Ao2s0_2'): A2c, sp.Symbol('S4o2s0_0'): S4R, sp.Symbol('Ao2s1_1'): 0}
key = (0, 2, 0)
h00 = sp.expand(sp.sympify(F.h[0][0].get(key, 0)).subs(sub))
A0 = sp.expand(sp.sympify(F.A[0].get(key, 0)).subs(sub))
print('near-zone order-2 slot-0 h00 :', h00, '   MP:', mp['h00'][0], '* R^2 +', mp['h00'][1], '* z^2')
print('near-zone order-2 slot-0 A0  :', A0, '   MP:', mp['A0'][0], '* R^2 +', mp['A0'][1], '* z^2')
hzz = sp.expand(sp.sympify(F.h[3][3].get(key, 0)).subs(sub)); hxx = sp.expand(sp.sympify(F.h[1][1].get(key, 0)).subs(sub))
print('h_zz:', hzz, ' h_xx:', hxx)
print('MP h_ij fit:', fit)
