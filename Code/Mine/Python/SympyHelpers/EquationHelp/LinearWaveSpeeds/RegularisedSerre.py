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
x,t,g,eps = symbols('x t g \epsilon',real=True)

#Serre equations
h = Function('h')(x,t)
u = Function('u')(x,t)



dhdt = diff(h,t)

fluxh = u*h
dfluxhdx = diff(fluxh,x)

LHSh = dhdt + dfluxhdx  

duhdt = diff(u*h,t)
Reg = h*(diff(u,x)**2 - diff(diff(u,x),t) -u*diff(u,(x,2))) -g*(h*diff(h,(x,2)) + (diff(h,x)**2)/2) 
dfluxuhdx = diff(u**2 *h  + g*h**2/2 + eps*h**2*Reg,x)

LHSuh = duhdt + dfluxuhdx 

H,U,delta = symbols('H U \delta', real = True)

n = Function('n')(x,t)
v = Function('v')(x,t)

hsub = H + delta*n + O(delta**2)
usub = U + delta*v + O(delta**2)

#get linearised equations
LHShlin = LHSh.subs(h, hsub).subs(u, usub).doit().expand()
LHShlin = LHShlin.subs(delta**2,0).removeO()

LHSuhlin = LHSuh.subs(h, hsub).subs(u, usub).doit().expand()
LHSuhlin = LHSuhlin.subs(delta**2,0).removeO()
LHSuhlin = LHSuhlin.subs(delta*diff(n,t), - diff(H*delta*v + U*delta*n,x)).doit().expand()

print('Linearised Equations')
print('h:')
display(LHShlin)
print('\n uh:')
display(LHSuhlin)
print('\n\n')

# wave structure 
A, B, k,w = symbols('A, B, \kappa \omega')
nsub = A*exp(I*(k*x - w*t))
vsub = B*exp(I*(k*x - w*t))

LHShlinwave = simplify(LHShlin.subs(delta,1).subs(n,nsub).subs(v,vsub).doit())
LHSuhlinwave = simplify(LHSuhlin.subs(delta,1).subs(n,nsub).subs(v,vsub).doit())

LHShlinewavecon = LHShlinwave.args[1]
LHSuhlinewavecon =LHSuhlinwave.args[2]

ACoeffh = LHShlinewavecon.subs(A,1).subs(B,0)
BCoeffh = LHShlinewavecon.subs(B,1).subs(A,0)

ACoeffuh = LHSuhlinewavecon.subs(A,1).subs(B,0)
BCoeffuh = LHSuhlinewavecon.subs(B,1).subs(A,0)

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
display(wavespeed[0])
display(wavespeed[1])