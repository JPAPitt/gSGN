from sympy import *
from sympy.solvers.solveset import linsolve

k = 3

qaj,qajm1,qajm2,qajm3,qajp1,qajp2,qajp3,dx = symbols('qA[j] qA[j-1] qA[j-2] qA[j-3] qA[j+1] qA[j+2] qA[j+3] dx')

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

#extract polynomials Coefficients
CubicSoljm3toj_Coeff = CubicSoljm3toj_Coeff.args[0]
CubicSoljm2tojp1_Coeff = CubicSoljm2tojp1_Coeff.args[0]
CubicSoljm1tojp2_Coeff = CubicSoljm1tojp2_Coeff.args[0]
CubicSoljtojp3_Coeff = CubicSoljtojp3_Coeff.args[0]

Polyjm3toj = CubicSoljm3toj_Coeff[0]*(x)**3 + CubicSoljm3toj_Coeff[1]*(x)**2 + CubicSoljm3toj_Coeff[2]*(x) + CubicSoljm3toj_Coeff[3]
Polyjm2tojp1 = CubicSoljm2tojp1_Coeff[0]*(x)**3 + CubicSoljm2tojp1_Coeff[1]*(x)**2 + CubicSoljm2tojp1_Coeff[2]*(x) + CubicSoljm2tojp1_Coeff[3]
Polyjm1tojp2 = CubicSoljm1tojp2_Coeff[0]*(x)**3 + CubicSoljm1tojp2_Coeff[1]*(x)**2 + CubicSoljm1tojp2_Coeff[2]*(x) + CubicSoljm1tojp2_Coeff[3]
Polyjtojp3 = CubicSoljtojp3_Coeff[0]*(x)**3 + CubicSoljtojp3_Coeff[1]*(x)**2 + CubicSoljtojp3_Coeff[2]*(x) + CubicSoljtojp3_Coeff[3]


# optimal weightings - just pick 1/4 for now
Sixjm3jp3_Coeff = linsolve([IntqCelljm3 -qajm3,IntqCelljm2 -qajm2, IntqCelljm1 -qajm1, IntqCellj  -qaj,IntqCelljp1 -qajp1, IntqCelljp2 -qajp2, IntqCelljp3 -qajp3], (qa,qb,qc,qd,qe,qf,qg ))
Sixjm3jp3_Coeff = Sixjm3jp3_Coeff.args[0]

Polyjm3jp3 =  Sixjm3jp3_Coeff[0]*x**6 + Sixjm3jp3_Coeff[1]*x**5 + Sixjm3jp3_Coeff[2]*x**4 + Sixjm3jp3_Coeff[3]*x**3 + Sixjm3jp3_Coeff[4]*x**2 + Sixjm3jp3_Coeff[5]*x + Sixjm3jp3_Coeff[6]


vjph_big = Polyjm3jp3.subs(x,dx/2)
vjph_jm3toj = Polyjm3toj.subs(x,dx/2)
vjph_jm2tojp1 = Polyjm2tojp1.subs(x,dx/2)
vjph_jm1tojp2 = Polyjm1tojp2.subs(x,dx/2)
vjph_jtojp3 =Polyjtojp3.subs(x,dx/2)

#Ideal Weights at cell edge
#Only  vjph_jm3toj contains jm3 information, so
cjm3toj =S('1/35')
cjtojp3 =S('4/35')

#Eliminate both
vjph_big_elimedge = vjph_big - cjm3toj*vjph_jm3toj  - cjtojp3*vjph_jtojp3

#Only  vjph_jm2toj contains jm2 information, so
cjm2tojp1 =S('12/35')
cjm1tojp2 =S('18/35')

vjph_big_elimedge_a =  vjph_big_elimedge - cjm2tojp1*vjph_jm2tojp1 - cjm1tojp2*vjph_jm1tojp2


# IdealLinearWeights = solve( [vjph_big  - (g1*vjph_jm3toj + g2*vjph_jm2tojp1 + g3*vjph_jm1tojp2 + g4*vjph_jtojp3 ), g1 + g2 + g3 + g4 - 1, g3>0, g4 > 0 ], (g1,g2,g3,g4) )