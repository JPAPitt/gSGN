from sympy import *
from numpy import matrix
from sympy.solvers.solveset import linsolve


def  Polynomial(yjmh,yjms,yjps,yjph,dx):
    a3 = (-9*yjmh + 27*yjms - 27*yjps + 9*yjph)/ (2*dx**3)
    a2 = (9*yjmh - 9*yjms - 9*yjps + 9*yjph  )/ (4*dx**2)
    a1 = (yjmh  - 27*yjms  + 27*yjps  - yjph  )/ (8*dx)
    a0 = (-yjmh+ 9*yjms + 9*yjps  - yjph)/ 16
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

ujmh, ujph,ujms, ujps = symbols('u_{j-1/2} u_{j+1/2} u_{j-1/6} u_{j+1/6}')

hjmh, hjph,hjms, hjps = symbols('hjmh hjph hjms hjps')


Geq = u*h - b1/2*diff(h**3*diff(u,x), x)


Pua3,Pua2,Pua1,Pua0 = Polynomial(ujmh,ujms,ujps,ujph,dx)
Pha3,Pha2,Pha1,Pha0 = Polynomial(hjmh,hjms,hjps,hjph,dx)

PuPoly = Pua3*x**3 + Pua2*x**2 + Pua1*x + Pua0 
PhPoly = Pha3*x**3 + Pha2*x**2 + Pha1*x + Pha0 


GeqPoly = Geq.subs(u,PuPoly).subs(h,PhPoly).doit()

GeqPolyjmh = GeqPoly.subs(x,-dx/2)
GeqPolyjmh = GeqPolyjmh.expand()

GeqPolyjmh = collect(GeqPolyjmh,ujmh)
GeqPolyjmh = collect(GeqPolyjmh,ujms)
GeqPolyjmh = collect(GeqPolyjmh,ujps)
GeqPolyjmh = collect(GeqPolyjmh,ujph)

GjmhujmhC = GeqPolyjmh.subs(ujmh,1).subs(ujms,0).subs(ujps,0).subs(ujph,0)
GjmhujmsC = GeqPolyjmh.subs(ujmh,0).subs(ujms,1).subs(ujps,0).subs(ujph,0)
GjmhujpsC = GeqPolyjmh.subs(ujmh,0).subs(ujms,0).subs(ujps,1).subs(ujph,0)
GjmhujphC = GeqPolyjmh.subs(ujmh,0).subs(ujms,0).subs(ujps,0).subs(ujph,1)


GeqPolyjms = GeqPoly.subs(x,-dx/6)
GeqPolyjms = GeqPolyjms.expand()

GeqPolyjms = collect(GeqPolyjms,ujmh)
GeqPolyjms = collect(GeqPolyjms,ujms)
GeqPolyjms = collect(GeqPolyjms,ujps)
GeqPolyjms = collect(GeqPolyjms,ujph)

GjmsujmhC = GeqPolyjms.subs(ujmh,1).subs(ujms,0).subs(ujps,0).subs(ujph,0)
GjmsujmsC = GeqPolyjms.subs(ujmh,0).subs(ujms,1).subs(ujps,0).subs(ujph,0)
GjmsujpsC = GeqPolyjms.subs(ujmh,0).subs(ujms,0).subs(ujps,1).subs(ujph,0)
GjmsujphC = GeqPolyjms.subs(ujmh,0).subs(ujms,0).subs(ujps,0).subs(ujph,1)

GeqPolyjps = GeqPoly.subs(x,dx/6)
GeqPolyjps = GeqPolyjps.expand()

GeqPolyjps = collect(GeqPolyjps,ujmh)
GeqPolyjps = collect(GeqPolyjps,ujms)
GeqPolyjps = collect(GeqPolyjps,ujps)
GeqPolyjps = collect(GeqPolyjps,ujph)

GjpsujmhC = GeqPolyjps.subs(ujmh,1).subs(ujms,0).subs(ujps,0).subs(ujph,0)
GjpsujmsC = GeqPolyjps.subs(ujmh,0).subs(ujms,1).subs(ujps,0).subs(ujph,0)
GjpsujpsC = GeqPolyjps.subs(ujmh,0).subs(ujms,0).subs(ujps,1).subs(ujph,0)
GjpsujphC = GeqPolyjps.subs(ujmh,0).subs(ujms,0).subs(ujps,0).subs(ujph,1)


GeqPolyjph = GeqPoly.subs(x,dx/2)
GeqPolyjph = GeqPolyjph.expand()

GeqPolyjph = collect(GeqPolyjph,ujmh)
GeqPolyjph = collect(GeqPolyjph,ujms)
GeqPolyjph = collect(GeqPolyjph,ujps)
GeqPolyjph = collect(GeqPolyjph,ujph)

GjphujmhC = GeqPolyjph.subs(ujmh,1).subs(ujms,0).subs(ujps,0).subs(ujph,0)
GjphujmsC = GeqPolyjph.subs(ujmh,0).subs(ujms,1).subs(ujps,0).subs(ujph,0)
GjphujpsC = GeqPolyjph.subs(ujmh,0).subs(ujms,0).subs(ujps,1).subs(ujph,0)
GjphujphC = GeqPolyjph.subs(ujmh,0).subs(ujms,0).subs(ujps,0).subs(ujph,1)