import sympy as sp
t,r,th,ph,eps=sp.symbols('t r theta phi epsilon')
M,Q,B=sp.symbols('M Q B',positive=True)
k=sp.symbols('k')
X=[t,r,th,ph]; n=4
f=1-2*M/r+Q**2/r**2
Y=(3*sp.cos(th)**2-1)/2
Sph=sp.sin(th)*sp.diff(Y,th)          # S_phi = sin(theta) Y_,theta
h0=r**3*f*B/3
u4=k*r**2*Q*B*(1-Q**2/r**2)           # ratio fixed by equations
g0=sp.diag(-f,1/f,r**2,r**2*sp.sin(th)**2)
h=sp.zeros(4,4); h[0,3]=h[3,0]=h0*Sph
A=[-Q/r,0,0,eps*u4*Sph]
g=g0+eps*h
gi=g.inv().applyfunc(lambda e: sp.series(e,eps,0,2).removeO())
def trunc(e): return sp.series(sp.expand(e),eps,0,2).removeO()
Gam=[[[trunc(sum(gi[a,s]*(sp.diff(g[s,b],X[c])+sp.diff(g[s,c],X[b])-sp.diff(g[b,c],X[s])) for s in range(n))/2) for c in range(n)] for b in range(n)] for a in range(n)]
def Ric(b,d):
    return trunc(sum(sp.diff(Gam[a][b][d],X[a])-sp.diff(Gam[a][b][a],X[d])+sum(Gam[a][a][e]*Gam[e][b][d]-Gam[a][d][e]*Gam[e][b][a] for e in range(n)) for a in range(n)))
F=[[sp.diff(A[b],X[a])-sp.diff(A[a],X[b]) for b in range(n)] for a in range(n)]
F2=trunc(sum(F[a][b]*gi[a,c]*gi[b,d]*F[c][d] for a in range(n) for b in range(n) for c in range(n) for d in range(n)))
def T(a,b):
    return trunc(sum(F[a][c]*gi[c,d]*F[b][d] for c in range(n) for d in range(n))-g[a,b]*F2/4)
out=[]
for (a,b) in [(0,3),(2,3),(1,3)]:
    e=sp.simplify(trunc(Ric(a,b)-2*T(a,b)).coeff(eps,1)); out.append(e)
sol=[sp.solve(sp.numer(sp.together(e)),k) if e!=0 else 'ok' for e in out]
print(sol)
