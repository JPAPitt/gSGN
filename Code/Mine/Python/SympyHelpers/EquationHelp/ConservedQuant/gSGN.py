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
x,t,g,beta1,beta2 = symbols('x t g \\beta_1 \\beta_2 ',real=True)

#Serre equations
h = Function('h')(x,t)
u = Function('u')(x,t)



dhdt = diff(h,t)

fluxh = u*h
dfluxhdx = diff(fluxh,x)

LHSh = dhdt + dfluxhdx  

duhdt = diff(u*h,t)
Phi = diff(u,x)**2 - diff(diff(u,x),t) - u*diff(u,(x,2))
Regh = h*diff(h,(x,2)) + (diff(h,x)**2)/2
Gamma = (1 + 3/2*beta1)*h*Phi - 3/2*beta2*g*Regh
dfluxuhdx = diff(u**2 *h  + g*h**2/2 + h**2*Gamma/3,x)

LHSuh = (duhdt + dfluxuhdx ).doit()


G = u*h - ((1 + 3/2*beta1))/3*diff(h**3*diff(u,x),x)
dGdt = diff(G,t)
dfluxGdx = diff(u*G + g*h**2/2 - 2/3*(1 + 3/2*beta1)*h**3*diff(u,x)**2 - 1/2*beta2*h**2*g*Regh,x)

LHSG = (dGdt + dfluxGdx).doit()

LHSdiffuhG = (LHSuh - LHSG).subs(diff(h,t), diff(-u*h,x)).doit()