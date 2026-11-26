"""Displacement of the reference point needed to cancel Delta:
delta_z_req = -Delta eps^2 / (d ln pred / dz),  pred = -(d sqrt(H)/dz)/(d Phi/dz) on the background."""
import mpmath as mp
from ab_first_num import fields
from ab_residual_scan import residual
mp.mp.dps = 50
e = mp.mpf(1)
for Q in [mp.mpf('1.2'), mp.mpf('1.5'), mp.mpf(2)]:
    for R in [mp.mpf(11), mp.mpf(21), mp.mpf(41), mp.mpf(81)]:
        l = R - 1
        if R <= Q**2 + 1:
            continue
        def lnpred(z):
            dH = mp.diff(lambda zz: fields(0, e, Q, l, 0, zz)[0], z)
            dP = mp.diff(lambda zz: fields(0, e, Q, l, 0, zz)[1], z)
            H = fields(0, e, Q, l, 0, z)[0]
            return mp.log(-(dH / (2 * mp.sqrt(H))) / dP)
        dl = mp.diff(lnpred, 0)
        D, H0, mu1 = residual(Q, l)
        dz = -D / dl
        Ebg = mp.diff(lambda zz: fields(0, e, Q, l, 0, zz)[1], 0)   # background field strength (Weyl z)
        print(f'Q={mp.nstr(Q,3)} R={mp.nstr(R,3)}  dlnpred/dz={mp.nstr(dl,8)}  dz_req/eps^2={mp.nstr(dz,10)}'
              f'  dz_req*R^2={mp.nstr(dz*R**2,8)}  mu/eps={mp.nstr(mu1,8)}', flush=True)
