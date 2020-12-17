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
# PuPoly = Pua3*x**3 + Pua2*x**2 + Pua1*x + Pua0 
# PhPoly = Pha3*x**3 + Pha2*x**2 + Pha1*x + Pha0 
uhA = integrate(PqPoly,(x,-dx/2,dx/2))
uhA = uhA.simplify()
uhA = uhA.subs(qjmh,ujmh*hjmh).subs(qjph,ujph*hjph).subs(dqjmh,dujmh*hjmh + ujmh*dhjmh).subs(dqjph,dujph*hjph + ujph*dhjph)


GeqInt = uhA - b1/2*(DivTermjph - DivTermjmh)
GeqInt = GeqInt.expand()
GeqInt = collect(GeqInt,ujmh)
GeqInt = collect(GeqInt,ujph)
GeqInt = collect(GeqInt,dujmh)
GeqInt = collect(GeqInt,dujph)


Cjujph = GeqInt.subs(ujmh,0).subs(dujmh,0).subs(dujph,0).subs(ujph,1)
Cjujmh = GeqInt.subs(ujmh,1).subs(dujmh,0).subs(dujph,0).subs(ujph,0)
Cjdujph = GeqInt.subs(ujmh,0).subs(dujmh,0).subs(dujph,1).subs(ujph,0)
Cjdujmh = GeqInt.subs(ujmh,0).subs(dujmh,1).subs(dujph,0).subs(ujph,0)

jCoeffVec = Matrix([[Cjujmh,Cjdujmh,Cjujph,Cjdujph]])
juVec = Matrix([[ujmh,dujmh,ujph,dujph]]).T

jm1CoeffVec = jCoeffVec.subs(hjmh,hjm3h).subs(dhjmh,dhjm3h).subs(hjph,hjmh).subs(dhjph,dhjmh)
jm1uVec = juVec.subs(ujmh,ujm3h).subs(dujmh,dujm3h).subs(ujph,ujmh).subs(dujph,dujmh)

jp1CoeffVec = jCoeffVec.subs(hjph,hjp3h).subs(dhjph,dhjp3h).subs(hjmh,hjph).subs(dhjmh,dhjph)
jp1uVec = juVec.subs(ujph,ujp3h).subs(dujph,dujp3h).subs(ujmh,ujph).subs(dujmh,dujph)

# Cjm1ujph = Cjujph.subs(hjph,hjmh).subs(dhjph,dhjmh)
# Cjm1ujmh = Cjujmh.subs(hjmh,hjm3h).subs(dhjmh,dhjm3h)
# Cjm1dujph = Cjdujph.subs(hjph,hjmh).subs(dhjph,dhjmh)
# Cjm1dujmh = Cjdujmh.subs(hjmh,hjm3h).subs(dhjmh,dhjm3h)

# Cjp1ujph = Cjujph.subs(hjph,hjp3h).subs(dhjph,dhjp3h)
# Cjp1ujmh = Cjujmh.subs(hjmh,hjph).subs(dhjmh,dhjph)
# Cjp1dujph = Cjdujph.subs(hjph,hjp3h).subs(dhjph,dhjp3h)
# Cjp1dujmh = Cjdujmh.subs(hjmh,hjph).subs(dhjmh,dhjph)


# # Cjm1ujmh = uhAC.subs(ujmh,1).subs(dujmh,0).subs(dujph,0).subs(ujph,0)
# # Cjm1dujph = uhAC.subs(ujmh,0).subs(dujmh,0).subs(dujph,1).subs(ujph,0)
# # Cjm1dujmh = uhAC.subs(ujmh,0).subs(dujmh,1).subs(dujph,0).subs(ujph,0)

# GVec =  Matrix([[Gajm1,Gaj,Gajp1]])

# hMat = Matrix([[Cjm1ujph  + Cjujph,Cjm1dujph  + Cjdujph],[Cjp1ujph  + Cjujph,Cjp1dujph  + Cjdujph]])
# uVec = Matrix([[ujmh,dujmh,ujph,dujph]])