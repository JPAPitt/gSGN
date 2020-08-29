from sympy import *
from sympy.solvers.solveset import linsolve

k = 4

qaj,qajm1,qajm2,qajm3,qajm4,qajp1,qajp2,qajp3,qajp4,dx = symbols('qA[j] qA[j-1] qA[j-2] qA[j-3] qA[j-4] qA[j+1] qA[j+2] qA[j+3] qA[j+4] dx')

pa,pb,pc,pd = symbols('p_a p_b p_c p_d')

x,xj = symbols('x x_j')

IntCell = integrate(pa*(x)**3 + pb*(x)**2 + pc*(x) + pd  ,x)

IntCelljm4 =  (IntCell.subs(x, -7*dx/2) -  IntCell.subs(x, -9*dx/2)) / dx
IntCelljm3 =  (IntCell.subs(x, -5*dx/2) -  IntCell.subs(x, -7*dx/2)) / dx
IntCelljm2 =  (IntCell.subs(x, -3*dx/2) -  IntCell.subs(x, -5*dx/2)) / dx

IntCelljm1 =  (IntCell.subs(x, -dx/2) -  IntCell.subs(x, -3*dx/2)) / dx
IntCellj =  (IntCell.subs(x, dx/2) -  IntCell.subs(x, -dx/2)) / dx
IntCelljp1 =  (IntCell.subs(x, 3*dx/2) -  IntCell.subs(x, dx/2)) / dx
IntCelljp2 =  (IntCell.subs(x, 5*dx/2) -  IntCell.subs(x, 3*dx/2)) / dx

IntCelljp3 =  (IntCell.subs(x, 7*dx/2) -  IntCell.subs(x, 5*dx/2)) / dx
IntCelljp4 =  (IntCell.subs(x, 9*dx/2) -  IntCell.subs(x, 7*dx/2)) / dx


# solve for polynomial coefficients

CubicSoljm4tojm1_Coeff = linsolve([IntCelljm4 -qajm4, IntCelljm3 -qajm3, IntCelljm2 -qajm2, IntCelljm1 -qajm1], (pa,pb,pc,pd ))
CubicSoljp1tojp4_Coeff = linsolve([IntCelljp1 -qajp1, IntCelljp2 -qajp2, IntCelljp3 -qajp3, IntCelljp4 -qajp4], (pa,pb,pc,pd ))

#extract polynomials Coefficients
CubicSoljm4tojm1_Coeff = CubicSoljm4tojm1_Coeff.args[0]
CubicSoljp1tojp4_Coeff =CubicSoljp1tojp4_Coeff.args[0]

CubicSoljm4tojm1 = CubicSoljm4tojm1_Coeff[0]*(x)**3 + CubicSoljm4tojm1_Coeff[1]*(x)**2  + CubicSoljm4tojm1_Coeff[2]*(x) + CubicSoljm4tojm1_Coeff[3]
CubicSoljp1tojp4 = CubicSoljp1tojp4_Coeff[0]*(x)**3 + CubicSoljp1tojp4_Coeff[1]*(x)**2  + CubicSoljp1tojp4_Coeff[2]*(x) + CubicSoljp1tojp4_Coeff[3]

Intersections = solve(CubicSoljm4tojm1 - CubicSoljp1tojp4,x)