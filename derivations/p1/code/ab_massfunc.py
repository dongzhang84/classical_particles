"""Hypothesis: the residual Delta' (rule 3') is the gradient force of a field-dependent mass,
   m_eff = mu [1 + kappa (mu^2 - e^2) (E_loc^2 - a_loc^2)],
with E_loc = proper electric field and a_loc = acceleration of the static observer at the hole.
Then T = sqrt(H) d_z ln m_eff / d_z sqrt(H)  ->  Delta' = kappa (mu1^2 - e^2) sqrt(H) d_z(E^2-a^2)/d_z sqrt(H)
(mu1 = mu/eps).  Background RN(m2=1, Q): f = 1 - 2/r + Q^2/r^2, E_loc = Q/r^2, a_loc = (r - Q^2)/(r^3 sqrt f).
Test: is Delta'/[ (mu1^2-e^2) sqrt(f) d_r(E^2-a^2)/d_r sqrt(f) ] the same constant for all (Q, R)?"""
import mpmath as mp
from ab_rule3p import residual
from ab_first_num import local_mu
mp.mp.dps = 50
e = mp.mpf(1)
for q, R in [('1.2', 11), ('1.2', 41), ('1.5', 11), ('1.5', 21), ('1.5', 41), ('2', 11), ('2', 41), ('3', 21)]:
    Q, R = mp.mpf(q), mp.mpf(R)
    if R <= Q**2 + 1:
        continue
    D = residual(Q, R - 1)[0]
    mu1 = local_mu(mp.mpf('1e-12'), e, Q, R - 1) / mp.mpf('1e-12')
    f = lambda r: 1 - 2 / r + Q**2 / r**2
    W = lambda r: (Q / r**2)**2 - ((r - Q**2) / (r**3 * mp.sqrt(f(r))))**2
    sf = lambda r: mp.sqrt(f(r))
    pred_shape = (mu1**2 - e**2) * sf(R) * mp.diff(W, R) / mp.diff(sf, R)
    print(f'Q={q} R={mp.nstr(R,3)}  Delta\'={mp.nstr(D,10)}  ratio kappa = {mp.nstr(D / pred_shape, 12)}', flush=True)
