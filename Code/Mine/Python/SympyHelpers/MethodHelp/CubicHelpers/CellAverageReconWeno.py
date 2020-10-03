from sympy import *
from sympy.solvers.solveset import linsolve

k = 3

qaj,qajm1,qajm2,qajm3,qajp1,qajp2,qajp3,dx = symbols('cellavgbc_q(i) cellavgbc_q(i-1) cellavgbc_q(i-2) cellavgbc_q(i-3) cellavgbc_q(i+1) cellavgbc_q(i+2) cellavgbc_q(i+3) dx')

pa,pb,pc,pd = symbols('p_a p_b p_c p_d')

qa,qb,qc,qd, qe, qf, qg = symbols('q_a q_b q_c q_d q_e q_f q_g')

g1,g2,g3,g4= symbols('g_1 g_2 g_3 g_4')

x,xj = symbols('x x_j')

IntCell = integrate(pa*(x)**3 + pb*(x)**2 + pc*(x) + pd  ,x)


IntCelljm3 =  (IntCell.subs(x, -5*dx/2) -  IntCell.subs(x, -7*dx/2)) / dx
IntCelljm2 =  (IntCell.subs(x, -3*dx/2) -  IntCell.subs(x, -5*dx/2)) / dx

IntCelljm1 =  (IntCell.subs(x, -dx/2) -  IntCell.subs(x, -3*dx/2)) / dx
IntCellj =  (IntCell.subs(x, dx/2) -  IntCell.subs(x, -dx/2)) / dx
IntCelljp1 =  (IntCell.subs(x, 3*dx/2) -  IntCell.subs(x, dx/2)) / dx
IntCelljp2 =  (IntCell.subs(x, 5*dx/2) -  IntCell.subs(x, 3*dx/2)) / dx

IntCelljp3 =  (IntCell.subs(x, 7*dx/2) -  IntCell.subs(x, 5*dx/2)) / dx

IntqCell = integrate(qa*(x)**6 +qb*(x)**5 + qc*(x)**4 + qd*(x)**3 + qe*(x)**2 + qf*(x) + qg  ,x)

IntqCelljm3 =  (IntqCell.subs(x, -5*dx/2) -  IntqCell.subs(x, -7*dx/2)) / dx
IntqCelljm2 =  (IntqCell.subs(x, -3*dx/2) -  IntqCell.subs(x, -5*dx/2)) / dx
IntqCelljm1 =  (IntqCell.subs(x, -dx/2) -  IntqCell.subs(x, -3*dx/2)) / dx
IntqCellj =  (IntqCell.subs(x, dx/2) -  IntqCell.subs(x, -dx/2)) / dx
IntqCelljp1 =  (IntqCell.subs(x, 3*dx/2) -  IntqCell.subs(x, dx/2)) / dx
IntqCelljp2 =  (IntqCell.subs(x, 5*dx/2) -  IntqCell.subs(x, 3*dx/2)) / dx
IntqCelljp3 =  (IntqCell.subs(x, 7*dx/2) -  IntqCell.subs(x, 5*dx/2)) / dx



# solve for polynomial coefficients

CubicSoljm3toj_Coeff = linsolve([IntCelljm3 -qajm3, IntCelljm2 -qajm2, IntCelljm1 -qajm1, IntCellj -qaj], (pa,pb,pc,pd ))
CubicSoljm2tojp1_Coeff = linsolve([IntCelljm2 -qajm2, IntCelljm1 -qajm1, IntCellj -qaj, IntCelljp1 -qajp1], (pa,pb,pc,pd ))
CubicSoljm1tojp2_Coeff = linsolve([IntCelljm1 -qajm1, IntCellj -qaj, IntCelljp1 -qajp1, IntCelljp2 -qajp2], (pa,pb,pc,pd ))
CubicSoljtojp3_Coeff = linsolve([IntCellj -qaj, IntCelljp1 -qajp1, IntCelljp2 -qajp2, IntCelljp3 -qajp3], (pa,pb,pc,pd ))

# #extract polynomials Coefficients
CubicSoljm3toj_Coeff = CubicSoljm3toj_Coeff.args[0]
CubicSoljm2tojp1_Coeff = CubicSoljm2tojp1_Coeff.args[0]
CubicSoljm1tojp2_Coeff = CubicSoljm1tojp2_Coeff.args[0]
CubicSoljtojp3_Coeff = CubicSoljtojp3_Coeff.args[0]

Polyjm3toj = CubicSoljm3toj_Coeff[0]*(x)**3 + CubicSoljm3toj_Coeff[1]*(x)**2 + CubicSoljm3toj_Coeff[2]*(x) + CubicSoljm3toj_Coeff[3]
Polyjm2tojp1 = CubicSoljm2tojp1_Coeff[0]*(x)**3 + CubicSoljm2tojp1_Coeff[1]*(x)**2 + CubicSoljm2tojp1_Coeff[2]*(x) + CubicSoljm2tojp1_Coeff[3]
Polyjm1tojp2 = CubicSoljm1tojp2_Coeff[0]*(x)**3 + CubicSoljm1tojp2_Coeff[1]*(x)**2 + CubicSoljm1tojp2_Coeff[2]*(x) + CubicSoljm1tojp2_Coeff[3]
Polyjtojp3 = CubicSoljtojp3_Coeff[0]*(x)**3 + CubicSoljtojp3_Coeff[1]*(x)**2 + CubicSoljtojp3_Coeff[2]*(x) + CubicSoljtojp3_Coeff[3]

# # #Smoothness indicators
Bjm3toj = 0
Bjm2tojp1 = 0
Bjm1tojp2 = 0
Bjtojp3 = 0
for l in range(1,k):
    Bjm3tojlterm = dx**(2*l - 1) *integrate( diff(Polyjm3toj,(x,l))**2, (x,-dx/2,dx/2))
    Bjm2tojp1lterm = dx**(2*l - 1) *integrate( diff(Polyjm2tojp1,(x,l))**2, (x,-dx/2,dx/2))
    Bjm1tojp2lterm = dx**(2*l - 1) *integrate( diff(Polyjm1tojp2,(x,l))**2, (x,-dx/2,dx/2))
    Bjtojp3lterm = dx**(2*l - 1) *integrate( diff(Polyjtojp3,(x,l))**2, (x,-dx/2,dx/2))
    
    Bjm3toj = Bjm3toj +  Bjm3tojlterm
    Bjm2tojp1 = Bjm2tojp1 +  Bjm2tojp1lterm
    Bjm1tojp2 = Bjm1tojp2 +  Bjm1tojp2lterm
    Bjtojp3 = Bjtojp3 +  Bjtojp3lterm
    
Bjm3toj = Bjm3toj.simplify()
Bjm2tojp1 =Bjm2tojp1.simplify()
Bjm1tojp2 =Bjm1tojp2.simplify()
Bjtojp3 =Bjtojp3.simplify()

print('pjm3toj')
print(CubicSoljm3toj_Coeff[0])
print(CubicSoljm3toj_Coeff[1])
print(CubicSoljm3toj_Coeff[2])
print(CubicSoljm3toj_Coeff[3])

print('pjm2tojp1')
print(CubicSoljm2tojp1_Coeff[0])
print(CubicSoljm2tojp1_Coeff[1])
print(CubicSoljm2tojp1_Coeff[2])
print(CubicSoljm2tojp1_Coeff[3])

print('pjm1tojp2')
print(CubicSoljm1tojp2_Coeff[0])
print(CubicSoljm1tojp2_Coeff[1])
print(CubicSoljm1tojp2_Coeff[2])
print(CubicSoljm1tojp2_Coeff[3])

print('pjtojp3')
print(CubicSoljtojp3_Coeff[0])
print(CubicSoljtojp3_Coeff[1])
print(CubicSoljtojp3_Coeff[2])
print(CubicSoljtojp3_Coeff[3])

print('Bs')
print(Bjm3toj)
print(Bjm2tojp1)
print(Bjm1tojp2)
print(Bjtojp3)

 