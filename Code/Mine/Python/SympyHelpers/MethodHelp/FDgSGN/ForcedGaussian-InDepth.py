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
x,t,ga,b1, b2 = symbols('x t ga b1a7 b2a7',real=True)

#Serre equations
h = Function('h')(x,t)
u = Function('u')(x,t)


#Mass conservation
dhdt = diff(h,t).doit()
Fluxh = u*h
dFluxhdx = diff(u*h,x).doit()

MassEq = dhdt + dFluxhdx 

#Momentum conservation
duhdt = diff(u*h,t).doit()
Fluxuh = u**2*h + ga*h**2/2 + h**2/2*(b1*h*( diff(u,x)**2 - diff(diff(u,x),t) - u*diff(u,(x,2)) ) - b2*ga*(h*diff(h,(x,2)) + diff(h,x)**2/2))
dFluxuhdx = diff(Fluxuh,x).doit()

MomeEq = duhdt + dFluxuhdx 
MomeEq = MomeEq.subs(diff(h,t),-dFluxhdx ).simplify()
uEq = MomeEq/h


a0,a1,a2,a3,a4, PPHI = symbols('a0,a1,a2,a3,a4, PPHI', positive = True, nonzero = True)

phi = x - a2*t
exp1 = exp(-(phi)**2 / (2*a3))
hforce = a0 + a1*exp1

uforce = a4*exp1

#have to be a bit careful because of large values
uEqForced = uEq.subs(h,hforce).subs(u,uforce)
uEqForced  = uEqForced.doit()
uEqForced1  =uEqForced.subs(exp1,'EXPPHI').subs(phi,PPHI).subs(exp(-PPHI**2/(2*a3)),'EXPPHI')
# uEqForced  = uEqForced.expand()
# uEqForced1  =uEqForced.subs(exp1,'EXPPHI1').subs(phi,'PHI')



