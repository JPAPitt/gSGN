"""
Python code to generate csv files that contain integrals of the basis functions
needed to generate the matrices of the Finite Element Method.

by Jordan Pitt 11/10/2018
"""

"""
############################# IMPORTS #########################################
"""


from sympy import * #Related website: https://www.sympy.org/en/index.html
from IPython.display import display

# symbols
x,t,dx,dt,ga,b1 = symbols('x,t,dx,dt,g,beta1', positive = True, nonzero = True,real=True)

h = Function('h')(x,t)
u = Function('u')(x,t)

# #Symbols
ha3,ha2,ha1,ha0 = symbols('a_h^3 a_h^2 a_h^1 a_h^0',real=True)
ua3,ua2,ua1,ua0 = symbols('a_u^3 a_u^2 a_u^1 a_u^0',real=True)
Ga3,Ga2,Ga1,Ga0 = symbols('a_G^3 a_G^2 a_G^1 a_G^0',real=True)
xmxj = symbols('xmxj', positive = True, nonzero = True,real=True)

Gaj = symbols('\\bar{G}_j',real=True)

Gjph,Gjps,Gjms,Gjmh,Gj,hjph,hjmh,hj  = symbols('G_{j+1/2} G_{j+1/6} G_{j-1/6} G_{j-1/2} G_j h_{j+1/2} h_{j-1/2} h_j',real=True)
dhjph,dhjmh,dhj  = symbols('h^d_{j+1/2} h^d_{j-1/2} h^d_{j}',real=True)
dujph,dujmh  = symbols('u^d_{j+1/2} u^d_{j-1/2}',real=True)

# #Terms

Geq = u*h - b1*diff(h**3*diff(u,x), x)
GeqInt = integrate(u*h,(x,-dx/2,dx/2)) - (b1*hjph**3*dujph - b1*hjmh**3*dujmh)
Geq =Geq.doit()

GPoly = Ga3*x**3+Ga2*x**2+Ga1*x + Ga0
hPoly = ha3*x**3+ha2*x**2+ha1*x + ha0
uPoly = ua3*x**3+ua2*x**2+ua1*x + ua0

GeqPolyLHS = Geq.subs(h,hPoly).subs(u,uPoly)
GeqPolyLHS = GeqPolyLHS.doit()
GeqPolyLHS = GeqPolyLHS.expand()
GeqPolyLHSPoly = collect(GeqPolyLHS,x)

GeqPolyLHSPolyInt = integrate(GeqPolyLHSPoly, (x,-dx/2,dx/2))
EQ1CA = GeqPolyLHSPolyInt - dx*Gaj

Geqjph = Geq.subs(diff(h,x),dhjph).subs(h,hjph).subs(u,uPoly).doit()
Geqjph = Geqjph.subs(x,dx/2)
Geqjmh = Geq.subs(diff(h,x),dhjmh).subs(h,hjmh).subs(u,uPoly).doit()
Geqjmh = Geqjmh.subs(x,-dx/2)

Geqj = Geq.subs(diff(h,x),dhj).subs(h,hj).subs(u,uPoly).doit()
Geqj = Geqj.subs(x,0)

GeqIntS = GeqInt.subs(integrate(u*h,(x,-dx/2,dx/2)) )

# uSolEdges = linsolve([Geqjph - Gjph, Geqj - Gj, GeqInt - Gaj, Geqjmh - Gjmh], (ua3,ua2,ua1,ua0 ))

# EQPjph = GeqPolyLHSPoly.subs(x,dx/2) - GPoly.subs(x,dx/2)
# EQPjmh = GeqPolyLHSPoly.subs(x,-dx/2) - GPoly.subs(x,-dx/2)


# dGPoly = diff(GPoly,x)
# dGeqPolyLHSPoly = diff(GeqPolyLHSPoly,x)
# EQdPjph = dGeqPolyLHSPoly.subs(x,dx/2) - dGPoly.subs(x,dx/2)
# EQdPjmh = dGeqPolyLHSPoly.subs(x,-dx/2) - dGPoly.subs(x,-dx/2)


# uSolEdges = linsolve([EQPjph, EQPjmh, EQdPjph, EQdPjmh], (ua3,ua2,ua1,ua0 ))

# GeqPolyLHS = Geq.subs(h,hPoly).subs(u,uPoly)