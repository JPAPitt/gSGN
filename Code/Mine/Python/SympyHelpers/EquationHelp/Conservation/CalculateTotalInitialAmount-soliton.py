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

ga,a0,a1,k,c = symbols('ga,a0,a1,k,c', positive = True, nonzero = True)
#k = sqrt(3*a1)/ (2*sqrt(a0*(a0 + a1)))
#c = sqrt(ga*(a0 + a1))
#h = a0 + a1*sech(k*(x - c*t))**2
#u = c*(1 - a0 / h)


h = Function("h")(x)
u = c*(1- a0 / h)

#Conservation mass
#Integral(sech(k*(-c*t + x))**2, (x, xbeg, xend)) = tanh(k*(x - c*t)) /k
#Sympy doesnt do this integral itself yet
hp1tot  = integrate(a0, (x,xbeg,xend), conds='none')
hp2tot = a1*tanh(k*(x - c*t)) / k
hp2tot = hp2tot.subs(x,xend) - hp2tot.subs(x,xbeg) 

htot = hp1tot + hp2tot

#Conservation G
G = u*h - diff(h**3*diff(u,x)/3,x)



