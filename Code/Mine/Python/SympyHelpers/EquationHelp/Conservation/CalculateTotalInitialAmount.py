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
x,a,b,g= symbols('x xbeg xend ga')
h = Function("h")(x)
u = Function("u")(x)

#derived functions
IGp2 = - diff(  h**3*diff(u,x)/3,x)


Hp1 = u*u*h
Hp2 = g*h*h
Hp3 = h**3 *diff(u,x)**2 

# H = (Hp1 + Hp2 + Hp3/3) / 2


#integrals
Totalh = integrate(h, (x,a,b))

#TotalGp1 = integrate(Gp1, x) - same as uh
TotalGp2 =  IGp2 #integrate(Gp2, (x,a,b))

Totalhu = integrate(u*h,(x,a,b))

#TotalH = integrate(H , x)


#Substitute functions
a0,a1,a2,a3,a4= symbols('a0 a1 a2 a3 a4', positive=True, finite = True)
hsub = a0 + a1*sin(a2*x)
usub = a3*sin(a4*x)

Totalh = Totalh.subs(h, hsub).doit()

TotalGP2Anti = IGp2.subs(h, hsub).subs(u, usub).doit()
TotalGp2 = TotalGP2Anti.subs(x,b) - TotalGP2Anti.subs(x,a)

##TotalHp1Anti = integrate(Hp1.subs(h, hsub).subs(u, usub),(x,a,b), conds='none').doit()
#print(1)
#TotalHp2Anti = integrate(Hp2.subs(h, hsub).subs(u, usub),(x,a,b), conds='none').doit()
#print(2)
TotalHp3Anti = Hp3.subs(h, hsub).subs(u, usub).doit()

TotalHp3Anti = integrate(TotalHp3Anti,(x,a,b), conds='none').doit()