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
x,t,xi,ga,c, a0, a1, A, C,k, tau = symbols('x t xi g c a_0 a_1 A C k \\tau',real=True)



#Serre equations
h = Function('h')(x,t)
u = Function('u')(x,t)
n = Function('n')(x,t)

#Travelling Wave Speed equations
np = a1*exp(-sqrt(3)*abs(x - c*t))
heq = 1 + np
ueq = c*n/(1 + n)

G = u*h - diff(h**3 / 3 * diff(u,x),x)

print('h')
display(heq)
print('u')
display(ueq)

print('dh/dx')
display(simplify(diff(heq,x)))

print('dh/dx')
display(simplify(diff(heq,(x,2))))

print('du/dx')
display(simplify(diff(ueq,x)))

print('d^2u/dx^2')
display(simplify(diff(ueq,(x,2))))



