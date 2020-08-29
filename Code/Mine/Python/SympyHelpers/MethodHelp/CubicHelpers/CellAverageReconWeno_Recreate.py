from sympy import *
from sympy.solvers.solveset import linsolve

k = 3

qaj,qajm1,qajm2,qajm3,qajp1,qajp2,qajp3,dx = symbols('qA[j] qA[j-1] qA[j-2] qA[j-3] qA[j+1] qA[j+2] qA[j+3] dx')

pa,pb,pc,pd = symbols('p_a p_b p_c p_d')

qa,qb,qc,qd, qe, qf, qg = symbols('q_a q_b q_c q_d q_e q_f q_g')

g1,g2,g3,g4= symbols('g_1 g_2 g_3 g_4')

x,xj = symbols('x x_j')

IntCell = integrate( pb*(x)**2 + pc*(x) + pd  ,x)


IntCelljm2 =  (IntCell.subs(x, -3*dx/2) -  IntCell.subs(x, -5*dx/2)) / dx

IntCelljm1 =  (IntCell.subs(x, -dx/2) -  IntCell.subs(x, -3*dx/2)) / dx
IntCellj =  (IntCell.subs(x, dx/2) -  IntCell.subs(x, -dx/2)) / dx
IntCelljp1 =  (IntCell.subs(x, 3*dx/2) -  IntCell.subs(x, dx/2)) / dx
IntCelljp2 =  (IntCell.subs(x, 5*dx/2) -  IntCell.subs(x, 3*dx/2)) / dx



# solve for polynomial coefficients

CubicSoljm2toj_Coeff = linsolve([IntCelljm2 -qajm2, IntCelljm1 -qajm1, IntCellj -qaj], (pb,pc,pd ))
CubicSoljm1tojp1_Coeff = linsolve([IntCelljm1 -qajm1, IntCellj -qaj, IntCelljp1 -qajp1], (pb,pc,pd ))
CubicSoljtojp2_Coeff = linsolve([IntCellj -qaj, IntCelljp1 -qajp1, IntCelljp2 -qajp2], (pb,pc,pd ))

# #extract polynomials Coefficients
CubicSoljm2toj_Coeff = CubicSoljm2toj_Coeff.args[0]
CubicSoljm1tojp1_Coeff = CubicSoljm1tojp1_Coeff.args[0]
CubicSoljtojp2_Coeff = CubicSoljtojp2_Coeff.args[0]

Polyjm2toj = CubicSoljm2toj_Coeff[0]*(x)**2 + CubicSoljm2toj_Coeff[1]*(x) + CubicSoljm2toj_Coeff[2]
Polyjm1tojp1 = CubicSoljm1tojp1_Coeff[0]*(x)**2 + CubicSoljm1tojp1_Coeff[1]*(x) + CubicSoljm1tojp1_Coeff[2]
Polyjtojp2 = CubicSoljtojp2_Coeff[0]*(x)**2 + CubicSoljtojp2_Coeff[1]*(x)+ CubicSoljtojp2_Coeff[2]

# # #Smoothness indicators
Bjm2toj = 0
Bjm1tojp1 = 0
Bjtojp2 = 0
for l in range(1,k):
    
    Bjm2tojterm = integrate( dx**(2*l - 1) *diff(Polyjm2toj,(x,l))**2, (x,-dx/2,dx/2))
    Bjm1tojp1term = integrate( dx**(2*l - 1) *diff(Polyjm1tojp1,(x,l))**2, (x,-dx/2,dx/2))
    Bjtojp2term = integrate( dx**(2*l - 1) *diff(Polyjtojp2,(x,l))**2, (x,-dx/2,dx/2))
    display(Bjm2tojterm.simplify())
    
    Bjm2toj = Bjm2toj +  Bjm2tojterm
    Bjm1tojp1 = Bjm1tojp1 +  Bjm1tojp1term
    Bjtojp2 = Bjtojp2 +  Bjtojp2term

    
Bjm2toj = Bjm2toj.simplify()
Bjm1tojp1 =Bjm1tojp1.simplify()
Bjtojp2 =Bjtojp2.simplify()

Bjm2toj_Notes = (13*(qajm2 - 2*qajm1 +qaj)**2) / 12 + (qajm2 - 4*qajm1 + 3*qaj)**2 /4

DiffB = Bjm2toj -  Bjm2toj_Notes

display( DiffB.simplify())
 