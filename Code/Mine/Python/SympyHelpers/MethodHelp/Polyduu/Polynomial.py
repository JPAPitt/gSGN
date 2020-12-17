from sympy import *
from sympy.solvers.solveset import linsolve

ujmh, ujph,dujmh, dujph = symbols('u_{j-1/2} u_{j+1/2} u^d_{j-1/2} u^d_{j+1/2}')
hjmh, hjph,dhjmh, dhjph = symbols('h_{j-1/2} h_{j+1/2} h^d_{j-1/2} h^d_{j+1/2}')
Gjmh, Gjph,dGjmh, dGjph = symbols('G_{j-1/2} G_{j+1/2} G^d_{j-1/2} G^d_{j+1/2}')

qjmh, qjph,dqjmh, dqjph = symbols('q_{j-1/2} q_{j+1/2} q^d_{j-1/2} q^d_{j+1/2}')
a3, a2,a1,a0 = symbols('a_3 a_2 a_1 a_0')
x, dx= symbols('x dx')

qPoly = a3*x**3 + a2*x**2 + a1*x + a0


dqPoly = diff(qPoly,x)

qPolyjmh = qPoly.subs(x,-dx/2)
dqPolyjmh = dqPoly.subs(x,-dx/2)

qPolyjph = qPoly.subs(x,dx/2)
dqPolyjph = dqPoly.subs(x,dx/2)

Poly_Coeff = linsolve([qPolyjmh - qjmh,qPolyjph - qjph,dqPolyjmh - dqjmh,dqPolyjph - dqjph], (a3, a2,a1,a0))

Poly_Coeffa3 = Poly_Coeff.args[0][0]
Poly_Coeffa2 = Poly_Coeff.args[0][1]
Poly_Coeffa1 = Poly_Coeff.args[0][2]
Poly_Coeffa0 = Poly_Coeff.args[0][3]

qPolySol = Poly_Coeffa3*x**3 + Poly_Coeffa2*x**2 + Poly_Coeffa1*x + Poly_Coeffa0

qPolySolInt = integrate(qPolySol,(x,-dx/2,dx/2))
# Poly_Coeffa2 = Poly_Coeff[a2]
# Poly_Coeffa1 = Poly_Coeff[a1]
# Poly_Coeffa0 = Poly_Coeff[a0]