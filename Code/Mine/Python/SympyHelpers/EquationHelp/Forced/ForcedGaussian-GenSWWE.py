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

phi = x - a2*t
exp1 = exp(-(phi)**2 / (2*a3))
h = a0 + a1*exp1
u = a4*exp1

G = u*h - eps*diff(h**3 *diff(u,x),x)

simpleG = simplify(G.subs(u,'u(i)').subs(h,'h(i)').subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')
print('G(i) = ' + str(expand(simpleG)))

## LHS's of equations
dhdt = diff(h,t).doit()

simpledhdt = simplify(dhdt.subs(u,'ui').subs(h,'hi').subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')

print()
print('dhdt = ' + str(expand(simpledhdt)))

dGdt = diff(G,t).doit()

simpledGdt= simplify(dGdt.subs(u,'ui').subs(h,'hi').subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')

print()
print('dGdt = ' + str(expand(simpledGdt)))

#RHS of equations

dFluxhdx = diff(u*h,x).doit()
simpledFluxhdx = simplify(dFluxhdx.subs(u,'ui').subs(h,'hi').subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')
print()
print('dfluxhdx = ' + str(expand(simpledFluxhdx)))

dFluxGdx= diff(u*G + (ga*h**2)/2 - (eps*h**2)*(2*h*diff(u,x)**2 \
        + ga*h*diff(h,(x,2)) + (ga/2)*diff(h,x)**2 ) ,x).doit()

simpledFluxGdx = simplify(dFluxGdx.subs(u,'ui').subs(h,'hi').subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')
print()
print('dfluxGdx = ' + str(expand(simpledFluxGdx)))

