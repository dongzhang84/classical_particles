"""Debug: is the order-3 system of ckpt2_axi (radial, resume) consistent?"""
import pickle, sys, sympy as sp
import __main__
import ckpt2_axi as C
__main__.Field = C.Field
import trunc as T
from sympy.polys.matrices import DomainMatrix
dd = pickle.load(open('ckpt2_axi_o2.pkl', 'rb'))
F = dd['F']; dcols = dd['log'][2]['free']
m = 3; unk = []; geqs = []
for k in range(0, C.EPS + 1):
    h, A, syms = C.ansatz(m, k)
    if k == 0:
        geqs += sum([C.reduce_eq(sum(C.xs[j] * h[i+1][j+1] for j in range(3))) for i in range(3)], [])
    F.add((k, m, 0), h, A); unk += syms
T.set_limits(C.EPS, m, 0)
eqs, tags = C.equations(F, m, tagged=True)
eqs += geqs; tags += [('radial', -1)] * len(geqs)
cols = unk + list(dcols)
Am, bm = sp.linear_eq_to_matrix(eqs, cols)
import random
gens = sorted((Am.free_symbols | bm.free_symbols), key=str)
for trial in range(2):
    num = {g: sp.Rational(random.randint(2, 19), random.randint(2, 19)) for g in gens}
    An = DomainMatrix.from_Matrix(Am.subs(num)).convert_to(sp.QQ)
    Abn = DomainMatrix.from_Matrix(Am.row_join(bm).subs(num)).convert_to(sp.QQ)
    print('numeric ranks A, [A|b]:', An.rank(), Abn.rank(), ' shape', Am.shape, flush=True)
    # which equations make it inconsistent: drop each tag group and recheck
    for Nmax in range(0, 4):
        keep = [i for i, tg in enumerate(tags) if tg[1] <= Nmax]
        a1 = DomainMatrix.from_Matrix(Am.extract(keep, list(range(Am.cols))).subs(num)).convert_to(sp.QQ).rank()
        a2 = DomainMatrix.from_Matrix(Am.row_join(bm).extract(keep, list(range(Am.cols + 1))).subs(num)).convert_to(sp.QQ).rank()
        print('   N <=', Nmax, 'consistent' if a1 == a2 else 'INCONSISTENT', flush=True)
    groups = []
    for gname in groups:
        keep = [i for i, tg in enumerate(tags) if tg != gname]
        a1 = DomainMatrix.from_Matrix(Am.extract(keep, list(range(Am.cols))).subs(num)).convert_to(sp.QQ).rank()
        a2 = DomainMatrix.from_Matrix(Am.row_join(bm).extract(keep, list(range(Am.cols + 1))).subs(num)).convert_to(sp.QQ).rank()
        if a1 == a2:
            print('   consistent without', gname, flush=True)
    break
