import sympy as sp
t,r,th,ph,eps=sp.symbols('t r theta phi epsilon')
M,Q,E=sp.symbols('M Q E',positive=True)
k,K0,K1,K2,K3=sp.symbols('k K0 K1 K2 K3')
X=[t,r,th,ph]
f=1-2*M/r+Q**2/r**2; D=r**2-2*M*r+Q**2
Y=(3*sp.cos(th)**2-1)/2
H0=-E*D; H2=H0
K=K0*r**2+K1*r+K2+K3/r
g0=sp.diag(-f,1/f,r**2,r**2*sp.sin(th)**2)
h=sp.diag(f*H0*Y, H2*Y/f, r**2*K*Y, r**2*sp.sin(th)**2*K*Y)
A=[-Q/r+eps*k*E*Q*D/r*Y,0,0,0]
g=g0+eps*h
gi=g.inv()
gi=gi.applyfunc(lambda e: sp.series(e,eps,0,2).removeO())
n=4
Gam=[[[sp.expand(sum(gi[a,s]*(sp.diff(g[s,b],X[c])+sp.diff(g[s,c],X[b])-sp.diff(g[b,c],X[s])) for s in range(n))/2) for c in range(n)] for b in range(n)] for a in range(n)]
def trunc(e): return sp.series(sp.expand(e),eps,0,2).removeO()
Gam=[[[trunc(Gam[a][b][c]) for c in range(n)] for b in range(n)] for a in range(n)]
def Ric(b,d):
    return trunc(sum(sp.diff(Gam[a][b][d],X[a])-sp.diff(Gam[a][b][a],X[d])+sum(Gam[a][a][e]*Gam[e][b][d]-Gam[a][d][e]*Gam[e][b][a] for e in range(n)) for a in range(n)))
F=[[sp.diff(A[b],X[a])-sp.diff(A[a],X[b]) for b in range(n)] for a in range(n)]
Fup=[[trunc(sum(gi[a,c]*gi[b,d]*F[c][d] for c in range(n) for d in range(n))) for b in range(n)] for a in range(n)]
F2=trunc(sum(F[a][b]*Fup[a][b] for a in range(n) for b in range(n)))
def T(a,b): # 4pi T_ab with Gaussian units: (F_ac F_b^c - 1/4 g F^2)
    Fmix=sum(F[a][c]*sum(gi[c,d]*F[b][d] for d in range(n)) for c in range(n))
    return trunc(Fmix-g[a,b]*F2/4)
eqs=[]
for (a,b) in [(0,0),(1,1),(2,2),(1,2),(3,3)]:
    e=trunc(Ric(a,b)-2*T(a,b))   # R_ab = 8pi T_ab = 2 (F F - 1/4 g F^2)
    eqs.append(sp.simplify(e.coeff(eps,1)))
sqrtg=sp.sqrt(-g.det()); 
# Maxwell t component: d_mu (sqrt(-g) F^{t mu}) 
sg=trunc(r**2*sp.sin(th)*(1+eps*sp.Rational(1,2)*sum(sp.Matrix(g0).inv()[a,a]*h[a,a] for a in range(4))))
mx=trunc(sum(sp.diff(sg*Fup[0][m],X[m]) for m in range(n)))
eqs.append(sp.simplify(mx.coeff(eps,1)))
unk=[k,K0,K1,K2,K3]
conds=[]
for e in eqs:
    num=sp.numer(sp.together(sp.expand(e.rewrite(sp.cos))))
    num=sp.expand(num.subs(sp.sin(th)**2,1-sp.cos(th)**2))
    P=sp.Poly(num, r, sp.cos(th))
    conds+=P.coeffs()
sol=sp.solve(conds,unk,dict=True)
print(sol)
