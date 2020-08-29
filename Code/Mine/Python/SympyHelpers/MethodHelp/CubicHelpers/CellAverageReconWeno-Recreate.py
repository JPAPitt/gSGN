from sympy import *
from sympy.solvers.solveset import linsolve

k = 3
dx = symbols('dx',positive=True)

qaj,qajm1,qajm2,qajp1,qajp2 = symbols('q_j q_{j-1} q_{j-2} q_{j+1} q_{j+2}')

pa,pb,pc= symbols('p_a p_b p_c')

x,xj = symbols('x x_j')

IntCell = integrate(pa*(x)**2 + pb*(x) + pc ,x)



IntCelljm2 =  (IntCell.subs(x, -3*dx/2) -  IntCell.subs(x, -5*dx/2)) / dx
IntCelljm1 =  (IntCell.subs(x, -dx/2) -  IntCell.subs(x, -3*dx/2)) / dx
IntCellj =  (IntCell.subs(x, dx/2) -  IntCell.subs(x, -dx/2)) / dx
IntCelljp1 =  (IntCell.subs(x, 3*dx/2) -  IntCell.subs(x, dx/2)) / dx
IntCelljp2 =  (IntCell.subs(x, 5*dx/2) -  IntCell.subs(x, 3*dx/2)) / dx


# solve for polynomial coefficients
CubicSoljm2toj_Coeff = linsolve([IntCelljm2 -qajm2, IntCelljm1 -qajm1, IntCellj -qaj], (pa,pb,pc))
CubicSoljm1tojp1_Coeff = linsolve([IntCelljm1 -qajm1, IntCellj -qaj, IntCelljp1 -qajp1], (pa,pb,pc))
CubicSoljtojp2_Coeff = linsolve([IntCellj -qaj, IntCelljp1 -qajp1, IntCelljp2 -qajp2], (pa,pb,pc))


#extract polynomials Coefficients
CubicSoljm2toj_Coeff = CubicSoljm2toj_Coeff.args[0]
CubicSoljm1tojp1_Coeff = CubicSoljm1tojp1_Coeff.args[0]
CubicSoljtojp2_Coeff = CubicSoljtojp2_Coeff.args[0]


Polyjm2toj = CubicSoljm2toj_Coeff[0]*(x)**2 + CubicSoljm2toj_Coeff[1]*(x) + CubicSoljm2toj_Coeff[2]
Polyjm1tojp1 = CubicSoljm1tojp1_Coeff[0]*(x)**2 + CubicSoljm1tojp1_Coeff[1]*(x) + CubicSoljm1tojp1_Coeff[2]
Polyjtojp2 = CubicSoljtojp2_Coeff[0]*(x)**2 + CubicSoljtojp2_Coeff[1]*(x) + CubicSoljtojp2_Coeff[2]



#l term
BJ = 0
for l in range(k):
    Bjlterm = dx**(2*l - 1) *integrate( diff(Polyjm2toj,(x,l))**2, (x,-dx/2,dx/2))
    BJ = BJ + Bjlterm
    print(Bjlterm)
    
BJ = BJ.simplify()