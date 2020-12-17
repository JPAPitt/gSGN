from sympy import *
from sympy.solvers.solveset import linsolve

k = 4

qaj,qajm1,qajm2,qajm3,qajp1,qajp2,qajp3,dx = symbols('qA[j] qA[j-1] qA[j-2] qA[j-3] qA[j+1] qA[j+2] qA[j+3] dx')

pa,pb,pc,pd = symbols('p_a p_b p_c p_d')


x,xj = symbols('x x_j')

IntCell = integrate(pa*(x)**3 + pb*(x)**2 + pc*(x) + pd  ,x)


IntCelljm2 =  (IntCell.subs(x, -3*dx/2) -  IntCell.subs(x, -5*dx/2)) / dx
IntCelljm1 =  (IntCell.subs(x, -dx/2) -  IntCell.subs(x, -3*dx/2)) / dx
IntCellj =  (IntCell.subs(x, dx/2) -  IntCell.subs(x, -dx/2)) / dx
IntCelljp1 =  (IntCell.subs(x, 3*dx/2) -  IntCell.subs(x, dx/2)) / dx
IntCelljp2 =  (IntCell.subs(x, 5*dx/2) -  IntCell.subs(x, 3*dx/2)) / dx



# solve for polynomial coefficients

CubicSoljm2tojp2_Coeff = linsolve([IntCelljm2 -qajm2, IntCelljm1 -qajm1, IntCelljp1 -qajp1, IntCelljp2 -qajp2], (pa,pb,pc,pd ))

# #extract polynomials Coefficients
CubicSoljm2tojp2_Coeff = CubicSoljm2tojp2_Coeff.args[0]



