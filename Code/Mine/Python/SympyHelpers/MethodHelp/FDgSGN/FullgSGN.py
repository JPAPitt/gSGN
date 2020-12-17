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
x,t,dx, dt,ga,b1,b2 = symbols('x,t,dx,dt,ga,\\beta_1,\\beta_2', positive = True, nonzero = True,real=True)

h = Function('h')(x,t)
u = Function('u')(x,t)

# #Terms

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

#Symbols
hinm1,him1nm1,him2nm1,hip1nm1,hip2nm1 =symbols('h_i^{n-1},h_{i-1}^{n-1},h_{i-2}^{n-1},h_{i+1}^{n-1},h_{i+2}^{n-1}',real=True)
hin,him1n,him2n,hip1n,hip2n =symbols('h_i^n,h_{i-1}^n,h_{i-2}^n,h_{i+1}^n,h_{i+2}^n',real=True)
hinp1,him1np1,him2np1,hip1np1,hip2np1 =symbols('h_i^{n+1},h_{i-1}^{n+1},h_{i-2}^{n+1},h_{i+1}^{n+1},h_{i+2}^{n+1}',real=True)

uinm1,uim1nm1,uim2nm1,uip1nm1,uip2nm1 =symbols('u_i^{n-1},u_{i-1}^{n-1},u_{i-2}^{n-1},u_{i+1}^{n-1},u_{i+2}^{n-1}',real=True)
uin,uim1n,uim2n,uip1n,uip2n =symbols('u_i^n,u_{i-1}^n,u_{i-2}^n,u_{i+1}^n,u_{i+2}^n',real=True)
uinp1,uim1np1,uim2np1,uip1np1,uip2np1 =symbols('u_i^{n+1},u_{i-1}^{n+1},u_{i-2}^{n+1},u_{i+1}^{n+1},u_{i+2}^{n+1}',real=True)

#Replace with approximations - Mass
MassEqFD = MassEq.subs(diff(h,t), (hinp1-hinm1) /(2*dt) ).subs(diff(h,x), (hip1n-him1n) /(2*dx) ).subs(diff(u,x), (uip1n-uim1n) /(2*dx) ).subs(u, uin).subs(h, hin)
hEvolve = solve(MassEqFD, hinp1)[0]
hEvolve = hEvolve.simplify()



# #Replace with approximations - Momentum
uEqFD = uEq.subs(diff(u,(x,2),t), ((uip1np1 - 2*uinp1  + uim1np1) - (uip1nm1 - 2*uinm1 + uim1nm1)) /(2*dx*dx*dt)).subs(diff(u,(x,3)), (-uim2n/2 +uim1n - uip1n + uip2n/2 ) /(dx**3)) \
.subs(diff(h,(x,3)), (-him2n/2 +him1n - hip1n + hip2n/2 ) /(dx**3)).subs(diff(u,(x,2)), (uim1n - 2*uin + uip1n ) /(dx**2)).subs(diff(h,(x,2)), (him1n - 2*hin + hip1n ) /(dx**2))\
.subs(diff(u,x,t),((uip1np1 - uim1np1) - (uip1nm1 - uim1nm1)) /(4*dx*dt)).subs(diff(u,x),((uip1n - uim1n)/(2*dx))).subs(diff(u,t),((uinp1 - uinm1)/(2*dt))).subs(diff(h,x),((hip1n - him1n)/(2*dx))).subs(u,uin).subs(h,hin)


uEqFDatN = uEqFD.subs(uinp1,0).subs(uip1np1,0).subs(uip2np1,0).subs(uim1np1,0).subs(uim2np1,0).subs(uim1nm1,0).subs(uinm1,0).subs(uip1nm1,0)

uEqFDNp1AndNm1 = uEqFD - uEqFDatN
uEqFDNp1AndNm1 = uEqFDNp1AndNm1.simplify()

uEqFDNm1 = uEqFDNp1AndNm1.subs(uinp1,0).subs(uip1np1,0).subs(uip2np1,0).subs(uim1np1,0).subs(uim2np1,0)

uEqFDNp1 = uEqFDNp1AndNm1.subs(uinm1,0).subs(uip1nm1,0).subs(uim1nm1,0)


uEqFDCheck = uEqFD - uEqFDatN - uEqFDNm1 - uEqFDNp1
uEqFDCheck = uEqFDCheck.simplify()