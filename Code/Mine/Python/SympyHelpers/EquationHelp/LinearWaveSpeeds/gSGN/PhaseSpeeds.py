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
g,beta1,beta2,k,w,u0,h0 = symbols(' g \\beta_1 \\beta_2 k \omega u_0 h_0',real=True)

#Serre equations
wp = u0*k + k*sqrt(g*h0)*sqrt((beta2*h0**2*k**2 + 2)/ ((S('2/3') + beta1)*h0**2*k**2 + 2 ))

wm = u0*k - k*sqrt(g*h0)*sqrt((beta2*h0**2*k**2 + 2)/ ((S('2/3') + beta1)*h0**2*k**2 + 2 ))

#phase speed
pvp = (wp /k)
pvm = (wm /k)

#SWWE phase speed

SWWpvp = pvp.subs(beta2,S('2/3') + beta1)
SWWpvm = pvm.subs(beta2,S('2/3') + beta1)

#Serre phase speed

Serrepvp = pvp.subs(beta1,0).subs(beta2,0)
Serrepvm = pvm.subs(beta1,0).subs(beta2,0)

#group speed
gvp = diff(wp,k)
gvm = diff(wm,k)

#SWWE group speed

SWWgvp = gvp.subs(beta2,S('2/3') + beta1)
SWWgvm = gvm.subs(beta2,S('2/3') + beta1)

#Serre group speed

Serregvp = gvp.subs(beta1,0).subs(beta2,0)
Serregvm = gvm.subs(beta1,0).subs(beta2,0)