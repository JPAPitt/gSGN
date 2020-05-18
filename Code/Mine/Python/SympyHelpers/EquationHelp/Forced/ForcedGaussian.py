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
x,t,xbeg,xend = symbols('x t xbeg xend',real=True)

ga,a0,a1,a2,a3,a4,a5,k,c = symbols('ga,a0,a1,a2,a3,a4,a5,k,c', positive = True, nonzero = True)
#k = sqrt(3*a1)/ (2*sqrt(a0*(a0 + a1)))
#c = sqrt(ga*(a0 + a1))
#h = a0 + a1*sech(k*(x - c*t))**2
#u = c*(1 - a0 / h)

phi = x - a2*t
exp1 = exp((phi - a3)**2 / (2*a4))
h = a0 + a1*exp1
u = a5*exp1

G = u*h - diff(h**3 *diff(u,x)**2/3,x)

## LHS's of equations
dhdt = diff(h,t).doit()

simpledhdt = dhdt.subs(exp1,'EXPPHI1').subs(phi,'PHI')



dGdt = diff(G,t).doit()

simpledGdt= dGdt.subs(h,'h(i)').subs(exp1,'EXPPHI1').subs(phi,'PHI').subs(phi,'PHI')


#RHS of equations

dFluxhdx = diff(u*h,x).doit()
simpledFluxhdx = dFluxhdx.subs(exp1,'EXPPHI1').subs(phi,'PHI')


dFluxGdx= diff(u*G + ga*h**2 - 2*h**3 /3 * diff(u,x),x).doit()

simpledFluxGdx = dFluxGdx.subs(h,'h(i)').subs(exp1,'EXPPHI1').subs(phi,'PHI').subs(phi,'PHI')