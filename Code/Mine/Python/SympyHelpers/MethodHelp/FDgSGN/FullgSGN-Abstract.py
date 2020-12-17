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
x,t,dx, dt,ga,b1,b2 = symbols('x,t,dx,dt,g,beta1,beta2', positive = True, nonzero = True,real=True)

h = Function('h')(x,t)
u = Function('u')(x,t)
p = Function('p')(x,t)

# #Symbols
hinm1,him1nm1,him2nm1,hip1nm1,hip2nm1 =symbols('h_i^{n-1},h_{i-1}^{n-1},h_{i-2}^{n-1},h_{i+1}^{n-1},h_{i+2}^{n-1}',real=True)
hin,him1n,him2n,hip1n,hip2n =symbols('h_i^n,h_{i-1}^n,h_{i-2}^n,h_{i+1}^n,h_{i+2}^n',real=True)
hinp1,him1np1,him2np1,hip1np1,hip2np1 =symbols('h_i^{n+1},h_{i-1}^{n+1},h_{i-2}^{n+1},h_{i+1}^{n+1},h_{i+2}^{n+1}',real=True)

uinm1,uim1nm1,uim2nm1,uip1nm1,uip2nm1 =symbols('u_i^{n-1},u_{i-1}^{n-1},u_{i-2}^{n-1},u_{i+1}^{n-1},u_{i+2}^{n-1}',real=True)
uin,uim1n,uim2n,uip1n,uip2n =symbols('u_i^n,u_{i-1}^n,u_{i-2}^n,u_{i+1}^n,u_{i+2}^n',real=True)
uinp1,uim1np1,uim2np1,uip1np1,uip2np1 =symbols('u_i^{n+1},u_{i-1}^{n+1},u_{i-2}^{n+1},u_{i+1}^{n+1},u_{i+2}^{n+1}',real=True)

chx,cux,chxx,cuxx,chxxx,cuxxx,pux, puxx, puxxx =symbols('chx cux chxx cuxx chxxx cuxxx pux puxx puxxxx',real=True)
hinf, uinf , uinm1f  =symbols('hbc(i) ubc(i) pubc(i)',real=True)


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

uEqLHS = uEq.subs(diff(u,t),0)
uEqLHS = uEqLHS.simplify()

uEqRHS = uEq - uEqLHS

uEqLHS = uEqLHS.subs(diff(h,(x,3)), chxxx).subs(diff(u,(x,3)), cuxxx) \
    .subs(diff(h,(x,2)), chxx).subs(diff(u,(x,2)), cuxx) \
    .subs(diff(h,(x,1)), chx).subs(diff(u,(x,1)), cux) \
    .subs(u,uinf).subs(h,hinf)


uEqRHS = uEqRHS.simplify()

uEqRHSnm1 = uEqRHS.subs(diff(u,t),-p/(2*dt) )
uEqRHSnm1 = uEqRHSnm1.doit()
uEqRHSnm1F =  uEqRHSnm1.subs(diff(h,(x,3)), chxxx).subs(diff(p,(x,3)), puxxx) \
    .subs(diff(h,(x,2)), chxx).subs(diff(p,(x,2)),puxx) \
    .subs(diff(h,(x,1)), chx).subs(diff(p,(x,1)), pux) \
    .subs(p,uinm1f).subs(h,hinf)

uEqRHSnp1 = uEqRHS.subs(diff(u,t),p/(2*dt) )
uEqRHSnp1  = uEqRHSnp1.doit()
uEqRHSnp1 = uEqRHSnp1.simplify()

uEqRHSnp1F =  uEqRHSnp1.subs(diff(h,(x,3)), chxxx).subs(diff(h,(x,2)), chxx) \
    .subs(diff(h,(x,1)), chx).subs(h,hinf)

uEqRHSnp1F = uEqRHSnp1F.subs(diff(p,(x,2)), (uim1np1 - 2*uinp1 + uip1np1 ) /(dx**2)).subs(diff(p,x), (uip1np1-uim1np1) /(2*dx) ).subs(p,uinp1)


uEqRHSnp1FDiag = uEqRHSnp1F.subs(uim1np1,0).subs(uip1np1,0).subs(uinp1,1)
uEqRHSnp1FDiag = uEqRHSnp1FDiag.simplify()

uEqRHSnp1FSubDiag = uEqRHSnp1F.subs(uim1np1,1).subs(uip1np1,0).subs(uinp1,0)
uEqRHSnp1FSubDiag = uEqRHSnp1FSubDiag.simplify()

uEqRHSnp1FSupDiag = uEqRHSnp1F.subs(uim1np1,0).subs(uip1np1,1).subs(uinp1,0)
uEqRHSnp1FSupDiag = uEqRHSnp1FSupDiag.simplify()
