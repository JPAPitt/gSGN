from sympy import *
from sympy.solvers.solveset import linsolve

qaj,qajm1,qajm2,qajm3,qajp1,qajp2,qajp3,dx = symbols('qa[j] qa[j-1] qa[j-2] qa[j-3] qa[j+1] qa[j+2] qa[j+3] dx')

pa,pb = symbols('p_a p_b')

x,xj = symbols('x x_j')

IntCell = integrate(pa*(x) + pb  ,x)



IntCelljm1 =  (IntCell.subs(x, -dx/2) -  IntCell.subs(x, -3*dx/2)) / dx
IntCellj =  (IntCell.subs(x, dx/2) -  IntCell.subs(x, -dx/2)) / dx
IntCelljp1 =  (IntCell.subs(x, 3*dx/2) -  IntCell.subs(x, dx/2)) / dx



LineSolvejm1 = linsolve([IntCelljm1 -qajm1, IntCellj -qaj], (pa,pb))
LineSolvejp1 = linsolve([IntCellj -qaj, IntCelljp1 -qajp1], (pa,pb))


# Polyjm2toj = CubicSoljm2toj_Coeff[0]*(x)**2 + CubicSoljm2toj_Coeff[1]*(x) + CubicSoljm2toj_Coeff[2]
# Polyjm1tojp1 = CubicSoljm1tojp1_Coeff[0]*(x)**2 + CubicSoljm1tojp1_Coeff[1]*(x) + CubicSoljm1tojp1_Coeff[2]
# Polyjtojp2 = CubicSoljtojp2_Coeff[0]*(x)**2 + CubicSoljtojp2_Coeff[1]*(x) + CubicSoljtojp2_Coeff[2]
# CubicSoljtojp3 = linsolve([IntCellj -qaj, IntCelljp1 -qajp1, IntCelljp2 -qajp2, IntCelljp3 -qajp3], (pa,pb,pc,pd ))