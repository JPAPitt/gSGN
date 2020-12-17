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
Gjmh, Gjph,Gjms, Gjps = symbols('Gjmh Gjph Gjms Gjps')
a, b = symbols('a b')

Geq = u*h - b1/2*diff(h**3*diff(u,x), x)



Pua3,Pua2,Pua1,Pua0 = Polynomial(ujmh,ujms,ujps,ujph,dx)
Pha3,Pha2,Pha1,Pha0 = Polynomial(hjmh,hjms,hjps,hjph,dx)
PGa3,PGa2,PGa1,PGa0 = Polynomial(Gjmh,Gjms,Gjps,Gjph,dx)

PuPoly = Pua3*x**3 + Pua2*x**2 + Pua1*x + Pua0 
PhPoly = Pha3*x**3 + Pha2*x**2 + Pha1*x + Pha0 
PGPoly = PGa3*x**3 + PGa2*x**2 + PGa1*x + PGa0 


GeqPoly = Geq.subs(u,PuPoly).subs(h,PhPoly).doit()

# UHA = integrate(u*h,(x,a,b))
# UHAPoly = UHA.subs(u,PuPoly).subs(h,PhPoly).doit()
# Divh3ux = h**3*diff(u,x) 
# Divh3uxUHAPoly = Divh3ux.subs(u,PuPoly).subs(h,PhPoly).doit()

# GeqIntjmhjms = UHAPoly.subs(a,-dx/2).subs(b,-dx/6) - b1/2*( Divh3uxUHAPoly.subs(x,-dx/6) - Divh3uxUHAPoly.subs(x,-dx/2))
# GeqIntjmhjms = GeqIntjmhjms.expand()
# GeqIntjmhjms = collect(GeqIntjmhjms,ujmh)
# GeqIntjmhjms = collect(GeqIntjmhjms,ujms)
# GeqIntjmhjms = collect(GeqIntjmhjms,ujps)
# GeqIntjmhjms = collect(GeqIntjmhjms,ujph)

# #b vector
GIntjmhjms  = PGPoly.integrate((x,-dx/2,-dx/6))
GIntjmhjms = GIntjmhjms.simplify()
GIntjmsjps  = PGPoly.integrate((x,-dx/6,dx/6))
GIntjmsjps = GIntjmsjps.simplify()
GIntjpsjph  = PGPoly.integrate((x,dx/6,dx/2))
GIntjpsjph = GIntjpsjph.simplify()

# # A matrix
# #j -1/2 to j-1/6
GeqPolyIntjjmhjjms = GeqPoly.integrate((x,-dx/2,-dx/6))
GeqPolyIntjjmhjjms = GeqPolyIntjjmhjjms.expand()
GeqPolyIntjjmhjjms = collect(GeqPolyIntjjmhjjms,ujmh)
GeqPolyIntjjmhjjms = collect(GeqPolyIntjjmhjjms,ujms)
GeqPolyIntjjmhjjms = collect(GeqPolyIntjjmhjjms,ujps)
GeqPolyIntjjmhjjms = collect(GeqPolyIntjjmhjjms,ujph)

GeqPolyIntjjmhjjmsujmhC = GeqPolyIntjjmhjjms.subs(ujmh,1).subs(ujms,0).subs(ujps,0).subs(ujph,0)
GeqPolyIntjjmhjjmsujmsC = GeqPolyIntjjmhjjms.subs(ujmh,0).subs(ujms,1).subs(ujps,0).subs(ujph,0)
GeqPolyIntjjmhjjmsujpsC = GeqPolyIntjjmhjjms.subs(ujmh,0).subs(ujms,0).subs(ujps,1).subs(ujph,0)
GeqPolyIntjjmhjjmsujphC = GeqPolyIntjjmhjjms.subs(ujmh,0).subs(ujms,0).subs(ujps,0).subs(ujph,1)


# #j -1/6 to j+1/6
GeqPolyIntjjmsjjps = GeqPoly.integrate((x,-dx/6,dx/6))
GeqPolyIntjjmsjjps = GeqPolyIntjjmsjjps.expand()
GeqPolyIntjjmsjjps = collect(GeqPolyIntjjmsjjps,ujmh)
GeqPolyIntjjmsjjps = collect(GeqPolyIntjjmsjjps,ujms)
GeqPolyIntjjmsjjps = collect(GeqPolyIntjjmsjjps,ujps)
GeqPolyIntjjmsjjps = collect(GeqPolyIntjjmsjjps,ujph)

GeqPolyIntjjmsjjpsujmhC = GeqPolyIntjjmsjjps.subs(ujmh,1).subs(ujms,0).subs(ujps,0).subs(ujph,0)
GeqPolyIntjjmsjjpsujmsC = GeqPolyIntjjmsjjps.subs(ujmh,0).subs(ujms,1).subs(ujps,0).subs(ujph,0)
GeqPolyIntjjmsjjpsujpsC = GeqPolyIntjjmsjjps.subs(ujmh,0).subs(ujms,0).subs(ujps,1).subs(ujph,0)
GeqPolyIntjjmsjjpsujphC = GeqPolyIntjjmsjjps.subs(ujmh,0).subs(ujms,0).subs(ujps,0).subs(ujph,1)


# #j +1/6 to j+1/2
GeqPolyIntjjpsjjph = GeqPoly.integrate((x,dx/6,dx/2))
GeqPolyIntjjpsjjph = GeqPolyIntjjpsjjph.expand()
GeqPolyIntjjpsjjph = collect(GeqPolyIntjjpsjjph,ujmh)
GeqPolyIntjjpsjjph = collect(GeqPolyIntjjpsjjph,ujms)
GeqPolyIntjjpsjjph = collect(GeqPolyIntjjpsjjph,ujps)
GeqPolyIntjjpsjjph = collect(GeqPolyIntjjpsjjph,ujph)

GeqPolyIntjjpsjjphujmhC = GeqPolyIntjjpsjjph.subs(ujmh,1).subs(ujms,0).subs(ujps,0).subs(ujph,0)
GeqPolyIntjjpsjjphujmsC = GeqPolyIntjjpsjjph.subs(ujmh,0).subs(ujms,1).subs(ujps,0).subs(ujph,0)
GeqPolyIntjjpsjjphujpsC = GeqPolyIntjjpsjjph.subs(ujmh,0).subs(ujms,0).subs(ujps,1).subs(ujph,0)
GeqPolyIntjjpsjjphujphC = GeqPolyIntjjpsjjph.subs(ujmh,0).subs(ujms,0).subs(ujps,0).subs(ujph,1)