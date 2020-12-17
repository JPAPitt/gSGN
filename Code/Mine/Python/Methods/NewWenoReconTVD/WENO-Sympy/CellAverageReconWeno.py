from sympy import *
from sympy.solvers.solveset import linsolve

k = 4

uaj,uajm1,uajm2,uajp1,uajp2,dx = symbols('\\bar{u}_i \\bar{u}_{i-1} \\bar{u}_{i-2} \\bar{u}_{i+1} \\bar{u}_{i+2} dx')

pa,pb,pc = symbols('p_a p_b p_c')
pia,pib,pic = symbols('p^*_a p^*_b p^*_c')

x,xmxj,xaj,xj,dxj,xmxjmdxj = symbols('x xmxj x^*_j x_j dx_j xmxjmdxj')

#what cells are we interpolating
#j-2,j-1,j
#want pc = bar{u}_j
#want pb = (bar{u}_{j} - bar{u}_{j-1})/dx
IntCell = integrate(pa*(xmxj)**2 + pb*(xmxj) + pc,xmxj)

#Cell integrals
IntCelljm2 =  (IntCell.subs(xmxj, -3*dx/2) -  IntCell.subs(xmxj, -5*dx/2)) / dx
IntCelljm1 =  (IntCell.subs(xmxj, -dx/2) -  IntCell.subs(xmxj, -3*dx/2)) / dx
IntCellj =  (IntCell.subs(xmxj, dx/2) -  IntCell.subs(xmxj, -dx/2)) / dx
IntCelljp1 =  (IntCell.subs(xmxj, 3*dx/2) -  IntCell.subs(xmxj, dx/2)) / dx
IntCelljp2 =  (IntCell.subs(xmxj, 5*dx/2) -  IntCell.subs(xmxj, 3*dx/2)) / dx


# # solve for polynomial coefficients
Polyjm2toj_Coeff = linsolve([IntCelljm2 -uajm2, IntCelljm1 -uajm1, IntCellj -uaj], (pa,pb,pc))
Polyjm1tojp1_Coeff = linsolve([IntCelljm1 -uajm1, IntCellj -uaj, IntCelljp1 -uajp1], (pa,pb,pc))
Polyjtojp2_Coeff = linsolve([IntCellj -uaj, IntCelljp1 -uajp1, IntCelljp2 -uajp2], (pa,pb,pc))

Polyjm2toj_Coeffs = Polyjm2toj_Coeff.args[0]  
Polyjm1tojp1_Coeffs = Polyjm1tojp1_Coeff.args[0]  
Polyjtojp2_Coeffs = Polyjtojp2_Coeff.args[0]  


#Re-arrange find solution - polynomial passes through cell average at FCAVjm2toj points (inside x_j pm dx/2)
FCAVjjm2toj = solve([Polyjm2toj_Coeffs[0]*(xmxj)**2 + Polyjm2toj_Coeffs[1]*(xmxj) + Polyjm2toj_Coeffs[2] - uaj ],(xmxj))
FCAVjjm2tojCheck0 = Polyjm2toj_Coeffs[0]*(FCAVjjm2toj[0][0])**2 + Polyjm2toj_Coeffs[1]*(FCAVjjm2toj[0][0]) + Polyjm2toj_Coeffs[2]
FCAVjjm2tojCheck0 = FCAVjjm2tojCheck0.simplify()
FCAVjjm2tojCheck1 = Polyjm2toj_Coeffs[0]*(FCAVjjm2toj[1][0])**2 + Polyjm2toj_Coeffs[1]*(FCAVjjm2toj[1][0]) + Polyjm2toj_Coeffs[2]
FCAVjjm2tojCheck1 = FCAVjjm2tojCheck1.simplify()

#Rewrite polynomial in right way
PolyF1 = pa*(x - xj)**2 + pb*(x - xj) + pc
PolyF1 = PolyF1.expand()

PolyF2 = pia*(x - xj - dxj)**2 + pib*(x - xj - dxj) + pic
PolyF2 = PolyF2.expand()

Polyn = solve([pa - pia, (pb) -(pib - 2*pia*dxj)],(pia,pib))

#Rewrite this poly:
#olyxmxjjm2toj = Polyjm2toj_Coeffs[0]*(xmxj)**2 + Polyjm2toj_Coeffs[1]*(xmxj) + Polyjm2toj_Coeffs[2]
# x^*_j - x_j =  FCAVjjm2toj

PCAjm2tojA = Polyjm2toj_Coeffs[0]
PCAjm2tojB = 2*PCAjm2tojA*(FCAVjjm2toj[0][0]) + Polyjm2toj_Coeffs[1]
PCAjm2tojB = PCAjm2tojB.simplify()
PCAjm2tojC = uaj

IntCelljm2toj = integrate(PCAjm2tojA*(xmxjmdxj)**2 + PCAjm2tojB*(xmxjmdxj) + PCAjm2tojC,xmxjmdxj)
IntCelljm2tojCjm2 =  (IntCelljm2toj.subs(xmxjmdxj,-3*dx/2-FCAVjjm2toj[0][0]) -  IntCelljm2toj.subs(xmxjmdxj,-5*dx/2-FCAVjjm2toj[0][0])) / dx
IntCelljm1tojp1Cjm2 =  (IntCelljm2toj.subs(xmxjmdxj,-dx/2-FCAVjjm2toj[0][0]) -  IntCelljm2toj.subs(xmxjmdxj,-3*dx/2-FCAVjjm2toj[0][0])) / dx
IntCelljtojp2Cjm2 =  (IntCelljm2toj.subs(xmxjmdxj,dx/2-FCAVjjm2toj[0][0]) -  IntCelljm2toj.subs(xmxjmdxj,-dx/2-FCAVjjm2toj[0][0])) / dx
