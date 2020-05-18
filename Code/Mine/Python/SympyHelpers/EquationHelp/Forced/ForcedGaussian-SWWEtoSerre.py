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

ga,a0,a1,a2,a3,a4,eta = symbols('ga,a0,a1,a2,a3,a4,eta', positive = True, nonzero = True)
#k = sqrt(3*a1)/ (2*sqrt(a0*(a0 + a1)))
#c = sqrt(ga*(a0 + a1))
#h = a0 + a1*sech(k*(x - c*t))**2
#u = c*(1 - a0 / h)

phi = x - a2*t
exp1 = exp(-(phi)**2 / (2*a3))
h = a0 + a1*exp1
u = a4*exp1

G = u*h - eta/3*diff(h**3 *diff(u,x),x)

simpleG = simplify(G.subs(u,'u(i)').subs(h,'h(i)').subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')

## LHS's of equations
dhdt = diff(h,t).doit()

simpledhdt = simplify(dhdt.subs(u,'ui').subs(h,'hi').subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')



dGdt = diff(G,t).doit()

simpledGdt= simplify(dGdt.subs(u,'ui').subs(h,'hi').subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')


#RHS of equations

dFluxhdx = diff(u*h,x).doit()
simpledFluxhdx = simplify(dFluxhdx.subs(u,'ui').subs(h,'hi').subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')


dFluxGdx= diff(u*G + (ga*h**2)/2 - (2*eta /3) * (h**3) * (diff(u,x)**2) ,x).doit()

simpledFluxGdx = simplify(dFluxGdx.subs(u,'ui').subs(h,'hi').subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')