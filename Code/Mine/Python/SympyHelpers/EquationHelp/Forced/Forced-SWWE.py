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

ga,a0,a1,a2,a3,a4,eps = symbols('ga,a0,a1,a2,a3,a4,eps', positive = True, nonzero = True)

#can come up with any smooth solutions which are functions in x,t

phi = x - a2*t
exp1 = exp(-(phi)**2 / (2*a3))
h = a0 + a1*exp1
u = a4*exp1

#momentum
p = u*h

#simplify and substitute, to use some obvious temporary variables
#ui = u*(x,t)
#hi = h*(x,t)
simpleG = simplify(p.subs(u,'u(i)').subs(h,'h(i)'))
print('p(i) = ' + str(expand(simpleG)))

## LHS's of equations

#Analytic value of dh*/dt
dhdt = diff(h,t).doit()

#simplify and substitute, to use some obvious temporary variables
#ui = u*(x,t)
#hi = h*(x,t)
#EXPPHI1 = exp(-(phi)**2 / (2*a3))
#-PHI = a2*t - x
simpledhdt = simplify(dhdt.subs(u,'ui').subs(h,'hi').subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')

print()
print('dhdt = ' + str(expand(simpledhdt)))

#Analytic value of d(p)*/dt
#simplify and substitute, to use some obvious temporary variables
dpdt = diff(p,t).doit()
simpledGdt= simplify(dpdt.subs(u,'ui').subs(h,'hi').subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')

print()
print('dGdt = ' + str(expand(simpledGdt)))

## LHS's of equations


#Analytic value of df(h)*/dx - the derivative of the flux function for mass equation
dFluxhdx = diff(u*h,x).doit()
simpledFluxhdx = simplify(dFluxhdx.subs(u,'ui').subs(h,'hi').subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')
print()
print('dfluxhdx = ' + str(expand(simpledFluxhdx)))

#Analytic value of df(uh)*/dx - the derivative of the flux function for momentum equation
dFluxGdx= diff(u*p + (ga*h**2)/2 ,x).doit()

simpledFluxGdx = simplify(dFluxGdx.subs(u,'ui').subs(h,'hi').subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')
print()
print('dfluxGdx = ' + str(expand(simpledFluxGdx)))

