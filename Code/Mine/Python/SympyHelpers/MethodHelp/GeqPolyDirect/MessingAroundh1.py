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
x = symbols('x')
t,dx,dt,ga,b1 = symbols('t,dx,dt,g,beta1', positive = True, nonzero = True,real=True)

def PolyFromPoints(yjmh,yjms,yjps,yjph,dx):
    a3 = (-9*yjmh + 27*yjms - 27*yjps + 9*yjph)/ (2*dx**3)
    a2 = (9*yjmh - 9*yjms - 9*yjps + 9*yjph  )/ (4*dx**2)
    a1 = (yjmh  - 27*yjms  + 27*yjps  - yjph  )/ (8*dx)
    a0 = (-yjmh+ 9*yjms + 9*yjps  - yjph)/ 16
    return a0,a1,a2,a3
    
h = Function('h')(x,t)
u = Function('u')(x,t)

# #Symbols
Ga3,Ga2,Ga1,Ga0 = symbols('Ga3 Ga2 Ga1 Ga0',real=True)
ha3,ha2,ha1,ha0 = symbols('ha3 ha2 ha1 ha0',real=True)
ua3,ua2,ua1,ua0 = symbols('ua3 ua2 ua1 ua0',real=True)

#FEM gives a a way to solve, but you need the basis to agree at boundaries. Maybe...

hPoly = ha1*(x) + ha0
uPoly = ua3*(x)**3 + ua2*(x)**2 + ua1*(x) + ua0
GPoly = Ga3*(x)**3 + Ga2*(x)**2 + Ga1*(x) + Ga0

#Want the polynomial to be at least 2nd order accurate for derivative
#Want conservation to be respected.
GPolyEq = hPoly*uPoly - b1*diff(hPoly**3 *diff(uPoly,x), x)
GPolyEq = GPolyEq.expand()
GPolyEqx = collect(GPolyEq,x) 

#Argument - don't need the higher terms, since dx is small.
# GPolyEqxCutoff = GPolyEqx.subs(x**4,0).subs(x**5,0).subs(x**6,0).subs(x**7,0).subs(x**8,0).subs(x**9,0)
GPolyEqxCutoff = GPolyEqx

GPolyEqxCutoffx0 = GPolyEqxCutoff.subs(x**9,0).subs(x**8,0).subs(x**7,0).subs(x**6,0).subs(x**5,0).subs(x**4,0).subs(x**3,0).subs(x**2,0).subs(x**1,0)

GPolyEqxCutoffx1 = GPolyEqxCutoff - GPolyEqxCutoffx0
GPolyEqxCutoffx1 = GPolyEqxCutoffx1.subs(x**9,0).subs(x**8,0).subs(x**7,0).subs(x**6,0).subs(x**5,0).subs(x**4,0).subs(x**3,0).subs(x**2,0).subs(x**1,1)

GPolyEqxCutoffx2 = GPolyEqxCutoff - GPolyEqxCutoffx0
GPolyEqxCutoffx2 = GPolyEqxCutoffx2.subs(x**9,0).subs(x**8,0).subs(x**7,0).subs(x**6,0).subs(x**5,0).subs(x**4,0).subs(x**3,0).subs(x**2,1).subs(x**1,0)

GPolyEqxCutoffx3 = GPolyEqxCutoff - GPolyEqxCutoffx0
GPolyEqxCutoffx3 = GPolyEqxCutoffx3.subs(x**9,0).subs(x**8,0).subs(x**7,0).subs(x**6,0).subs(x**5,0).subs(x**4,0).subs(x**3,1).subs(x**2,0).subs(x**1,0)


# GPolyEqxCutoffx4 = GPolyEqxCutoff - GPolyEqxCutoffx0
# GPolyEqxCutoffx4 = GPolyEqxCutoffx3.subs(x**9,0).subs(x**8,0).subs(x**7,0).subs(x**6,0).subs(x**5,0).subs(x**4,1).subs(x**3,0).subs(x**2,0).subs(x**1,0)

# GPolyEqxCutoffx5 = GPolyEqxCutoff - GPolyEqxCutoffx0
# GPolyEqxCutoffx5 = GPolyEqxCutoffx5.subs(x**9,0).subs(x**8,0).subs(x**7,0).subs(x**6,0).subs(x**5,1).subs(x**4,0).subs(x**3,0).subs(x**2,0).subs(x**1,0)

# GPolyEqxCutoffx6 = GPolyEqxCutoff - GPolyEqxCutoffx0
# GPolyEqxCutoffx6 = GPolyEqxCutoffx6.subs(x**9,0).subs(x**8,0).subs(x**7,0).subs(x**6,1).subs(x**5,0).subs(x**4,0).subs(x**3,0).subs(x**2,0).subs(x**1,0)

# GPolyEqxCutoffx7 = GPolyEqxCutoff - GPolyEqxCutoffx0
# GPolyEqxCutoffx7 = GPolyEqxCutoffx5.subs(x**9,0).subs(x**8,0).subs(x**7,1).subs(x**6,0).subs(x**5,0).subs(x**4,0).subs(x**3,0).subs(x**2,0).subs(x**1,0)

# GPolyEqxCutoffx8 = GPolyEqxCutoff - GPolyEqxCutoffx0
# GPolyEqxCutoffx8 = GPolyEqxCutoffx8.subs(x**9,0).subs(x**8,1).subs(x**7,0).subs(x**6,0).subs(x**5,0).subs(x**4,0).subs(x**3,0).subs(x**2,0).subs(x**1,0)

# GPolyEqxCutoffx9 = GPolyEqxCutoff - GPolyEqxCutoffx0
# GPolyEqxCutoffx9 = GPolyEqxCutoffx9.subs(x**9,1).subs(x**8,0).subs(x**7,0).subs(x**6,0).subs(x**5,0).subs(x**4,0).subs(x**3,0).subs(x**2,0).subs(x**1,0)

uCoeffs = solve([GPolyEqxCutoffx0 -Ga0 ,GPolyEqxCutoffx1 -Ga1,GPolyEqxCutoffx2 -Ga2, GPolyEqxCutoffx3 -Ga3], [ua3,ua2,ua1,ua0])
ua3Sol = uCoeffs[ua3]
ua2Sol = uCoeffs[ua2]
ua1Sol = uCoeffs[ua1]
ua0Sol = uCoeffs[ua0]