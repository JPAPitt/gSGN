from sympy import *
from numpy import matrix
from sympy.solvers.solveset import linsolve

def Polynomial(ujmh,ujph,dujmh,dujph,dx):
    a3 = (dx*dujph + dx*dujmh - 2*ujph + 2*ujmh)/dx**3
    a2 = (dujph- dujmh)/(2*dx)
    a1 = (-dx*dujph - dx*dujmh + 6*ujph - 6*ujmh)/(4*dx)
    a0 = -dx*dujph/8 + dx*dujmh/8 + ujph/2 + ujmh/2
    
    return a3,a2,a1,a0



# symbols
x,t,dx,dt,ga,b1 = symbols('x,t,dx,dt,g,beta1', positive = True, nonzero = True,real=True)

h = Function('h')(x,t)
u = Function('u')(x,t)

# #Symbols
ha3,ha2,ha1,ha0 = symbols('a_h^3 a_h^2 a_h^1 a_h^0',real=True)
ua3,ua2,ua1,ua0 = symbols('a_u^3 a_u^2 a_u^1 a_u^0',real=True)
Ga3,Ga2,Ga1,Ga0 = symbols('a_G^3 a_G^2 a_G^1 a_G^0',real=True)
xmxj = symbols('xmxj', positive = True, nonzero = True,real=True)

ujmh, ujph,dujmh, dujph = symbols('u_{j-1/2} u_{j+1/2} u^d_{j-1/2} u^d_{j+1/2}')
ujp3h,dujp3h = symbols('u_{j+3/2} u^d_{j+3/2}')
ujm3h,dujm3h = symbols('u_{j-3/2} u^d_{j-3/2}')

hjmh, hjph,dhjmh, dhjph = symbols('hjmh hjph dhjmh dhjph')
hjp3h,dhjp3h = symbols('hjp3h dhjp3h ')
hjm3h,dhjm3h = symbols('hjm3h dhjm3h')
Gaj,Gajp1,Gajm1 = symbols('Gaj Gajp1 Gajm1')

qjmh, qjph,dqjmh, dqjph = symbols('q_{j-1/2} q_{j+1/2} q^d_{j-1/2} q^d_{j+1/2}')

Geq = u*h - b1/2*diff(h**3*diff(u,x), x)

DivTerm = h**3*diff(u,x)
DivTerm = DivTerm.doit()
DivTermjph = DivTerm.subs(h,hjph).subs(diff(u,x),dujph).subs(u,ujph)
DivTermjmh = DivTerm.subs(h,hjmh).subs(diff(u,x),dujmh).subs(u,ujmh)


Pua3,Pua2,Pua1,Pua0 = Polynomial(ujmh,ujph,dujmh,dujph,dx)
Pha3,Pha2,Pha1,Pha0 = Polynomial(hjmh,hjph,dhjmh,dhjph,dx)
Pqa3,Pqa2,Pqa1,Pqa0 = Polynomial(qjmh,qjph,dqjmh,dqjph,dx)

PqPoly = Pqa3*x**3 + Pqa2*x**2 + Pqa1*x + Pqa0 
PuPoly = Pua3*x**3 + Pua2*x**2 + Pua1*x + Pua0 
PhPoly = Pha3*x**3 + Pha2*x**2 + Pha1*x + Pha0 
uhA = integrate(PqPoly,(x,-dx/2,dx/2))
uhA = uhA.simplify()
uhA = uhA.subs(qjmh,ujmh*hjmh).subs(qjph,ujph*hjph).subs(dqjmh,dujmh*hjmh + ujmh*dhjmh).subs(dqjph,dujph*hjph + ujph*dhjph)

#Cell average eq
GeqInt = uhA - b1/2*(DivTermjph - DivTermjmh)
GeqInt = GeqInt.expand()
GeqInt = collect(GeqInt,ujmh)
GeqInt = collect(GeqInt,ujph)
GeqInt = collect(GeqInt,dujmh)
GeqInt = collect(GeqInt,dujph)


CGajujph = GeqInt.subs(ujmh,0).subs(dujmh,0).subs(dujph,0).subs(ujph,1)
CGajujmh = GeqInt.subs(ujmh,1).subs(dujmh,0).subs(dujph,0).subs(ujph,0)
CGajdujph = GeqInt.subs(ujmh,0).subs(dujmh,0).subs(dujph,1).subs(ujph,0)
CGajdujmh = GeqInt.subs(ujmh,0).subs(dujmh,1).subs(dujph,0).subs(ujph,0)

#Cell average d/dx (eq)
dUHA = ujph*hjph - ujmh*hjmh
dh3uxjph =( hjph**3*diff(PuPoly,(x,2)).subs(x,dx/2) + 3*hjph**2*dhjph*dujph )
dh3uxjph = dh3uxjph.simplify()
dh3uxjmh =( hjmh**3*diff(PuPoly,(x,2)).subs(x,-dx/2) + 3*hjmh**2*dhjmh*dujmh )
dh3uxjmh = dh3uxjmh.simplify()
dGeqInt = dUHA - b1/2*( dh3uxjph-dh3uxjmh)
dGeqInt = dGeqInt.expand()
dGeqInt = collect(dGeqInt,ujmh)
dGeqInt = collect(dGeqInt,ujph)
dGeqInt = collect(dGeqInt,dujmh)
dGeqInt = collect(dGeqInt,dujph)

CdGajujph = dGeqInt.subs(ujmh,0).subs(dujmh,0).subs(dujph,0).subs(ujph,1)
CdGajujmh = dGeqInt.subs(ujmh,1).subs(dujmh,0).subs(dujph,0).subs(ujph,0)
CdGajdujph = dGeqInt.subs(ujmh,0).subs(dujmh,0).subs(dujph,1).subs(ujph,0)
CdGajdujmh = dGeqInt.subs(ujmh,0).subs(dujmh,1).subs(dujph,0).subs(ujph,0)

