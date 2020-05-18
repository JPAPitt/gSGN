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
x,t,g,beta1,beta2 = symbols('x t g \\beta_1 \\beta_2',real=True)

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
Gamma = (1 + 3*beta1/2)*h*Phi - 3*beta2*g*Regh/2
dfluxuhdx = diff(u**2 *h  + g*h**2/2 +h**2*Gamma/3,x)

LHSuh = duhdt + dfluxuhdx

h0,u0,delta = symbols('h_0 u_0 \delta', real = True)

n = Function('n')(x,t)
v = Function('v')(x,t)

hsub = h0 + delta*n + O(delta**2)
usub = u0 + delta*v + O(delta**2)

#get linearised equations
LHShlin = LHSh.subs(h, hsub).subs(u, usub).doit().expand()
LHShlin = LHShlin.subs(delta**2,0).removeO()

LHSuhlin = LHSuh.subs(h, hsub).subs(u, usub).doit().expand()
LHSuhlin = LHSuhlin.subs(delta**2,0).removeO()
LHSuhlin = LHSuhlin.subs(delta*diff(n,t), - diff(h0*delta*v + u0*delta*n,x)).doit().expand()

print('Linearised Equations')
print('h:')
display(LHShlin)
print('\n uh:')
display(LHSuhlin)
print('\n\n')

# wave structure 
H,U, k,w = symbols('H, U, k, \omega')
nsub = H*exp(I*(k*x - w*t))
vsub = U*exp(I*(k*x - w*t))

LHShlinwave = simplify(LHShlin.subs(delta,1).subs(n,nsub).subs(v,vsub).doit())
LHSuhlinwave = simplify(LHSuhlin.subs(delta,1).subs(n,nsub).subs(v,vsub).doit())

LHShlinewavecon = LHShlinwave.args[1]
LHSuhlinewavecon =LHSuhlinwave.args[-2]


ACoeffh = LHShlinewavecon.subs(H,1).subs(U,0)
BCoeffh = LHShlinewavecon.subs(U,1).subs(H,0)

ACoeffuh = LHSuhlinewavecon.subs(H,1).subs(U,0)
BCoeffuh = LHSuhlinewavecon.subs(U,1).subs(H,0)


CoeffMatrix = Matrix([[ACoeffh,BCoeffh],[ACoeffuh,BCoeffuh]])

DetCoeffMatrix = CoeffMatrix.det()
DetCoeffMatC = DetCoeffMatrix.subs(w,0).expand()
DetCoeffMatB = collect((DetCoeffMatrix - DetCoeffMatC).expand(),w).subs(w**2,0).subs(w,1).expand()
DetCoeffMatA = collect((DetCoeffMatrix - DetCoeffMatB*w - DetCoeffMatC ).expand(),w).subs(w**2,1)

QuadSolp = (-DetCoeffMatB +sqrt(DetCoeffMatB**2 - 4*DetCoeffMatA*DetCoeffMatC))/ (2*DetCoeffMatA)
QuadSolm = (-DetCoeffMatB -sqrt(DetCoeffMatB**2 - 4*DetCoeffMatA*DetCoeffMatC))/ (2*DetCoeffMatA)

QuadSolp = QuadSolp.expand()
QuadSolm = QuadSolm.expand()

wavespeed = solve(DetCoeffMatrix,w)

print('Coefficient Matrix')
display(CoeffMatrix)
print('Determinant')
display(collect(DetCoeffMatrix.expand(),w))
print('Quadratic solution')
display(QuadSolp)
display(QuadSolm)
print('wavespeeds')
display(wavespeed[0].simplify())
display(wavespeed[1].simplify())
