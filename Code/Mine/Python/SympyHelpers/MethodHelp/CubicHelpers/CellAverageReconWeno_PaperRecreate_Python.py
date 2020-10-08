from sympy import *
from sympy.solvers.solveset import linsolve

k = 3

vi,vim1,vim2,vip1,vip2,dx = symbols('vi vim1 vim2 vip1 vip2 dx')

pa,pb,pc = symbols('p_a p_b p_c')


x,xj = symbols('x x_j')

IntCell = integrate(pa*(x)**2 + pb*(x) + pc  ,x)

IntCelljm2 =  (IntCell.subs(x, -3*dx/2) -  IntCell.subs(x, -5*dx/2)) / dx
IntCelljm1 =  (IntCell.subs(x, -dx/2) -  IntCell.subs(x, -3*dx/2)) / dx
IntCellj =  (IntCell.subs(x, dx/2) -  IntCell.subs(x, -dx/2)) / dx
IntCelljp1 =  (IntCell.subs(x, 3*dx/2) -  IntCell.subs(x, dx/2)) / dx
IntCelljp2 =  (IntCell.subs(x, 5*dx/2) -  IntCell.subs(x, 3*dx/2)) / dx



# solve for polynomial coefficients

CubicSoljm2toj_Coeff = linsolve([ IntCelljm2 -vim2, IntCelljm1 -vim1, IntCellj -vi], (pa,pb,pc ))
CubicSoljm1tojp1_Coeff = linsolve([ IntCelljm1 -vim1, IntCellj -vi, IntCelljp1 -vip1], (pa,pb,pc ))
CubicSoljtojp2_Coeff = linsolve([ IntCellj -vi, IntCelljp1 -vip1, IntCelljp2 -vip2], (pa,pb,pc ))


# #extract polynomials Coefficients
CubicSoljm2toj_Coeff = CubicSoljm2toj_Coeff.args[0]
CubicSoljm1tojp1_Coeff  = CubicSoljm1tojp1_Coeff.args[0]
CubicSoljtojp2_Coeff = CubicSoljtojp2_Coeff.args[0]

Polyjm2toj = CubicSoljm2toj_Coeff[0]*(x)**2 + CubicSoljm2toj_Coeff[1]*(x) + CubicSoljm2toj_Coeff[2]
Polyjm1tojp1 = CubicSoljm1tojp1_Coeff[0]*(x)**2 + CubicSoljm1tojp1_Coeff[1]*(x) + CubicSoljm1tojp1_Coeff[2]
Polyjtojp2 = CubicSoljtojp2_Coeff[0]*(x)**2 + CubicSoljtojp2_Coeff[1]*(x) + CubicSoljtojp2_Coeff[2]


# # #Smoothness indicators
Bjm2toj = 0
Bjm1tojp1 = 0
Bjtojp2 = 0
for l in range(1,k):
    Bjm2tojlterm = dx**(2*l - 1) *integrate( diff(Polyjm2toj,(x,l))**2, (x,-dx/2,dx/2))
    Bjm1tojp1lterm = dx**(2*l - 1) *integrate( diff(Polyjm1tojp1,(x,l))**2, (x,-dx/2,dx/2))
    Bjtojp2lterm = dx**(2*l - 1) *integrate( diff(Polyjtojp2,(x,l))**2, (x,-dx/2,dx/2))

    
    Bjm2toj = Bjm2toj +  Bjm2tojlterm
    Bjm1tojp1 = Bjm1tojp1 +  Bjm1tojp1lterm
    Bjtojp2 = Bjtojp2 +  Bjtojp2lterm
    
Bjm2toj = Bjm2toj.simplify()
Bjm1tojp1 =Bjm1tojp1.simplify()
Bjtojp2 =Bjtojp2.simplify()


#Print same

# Polyjm2toj.subs(x,dx/2)