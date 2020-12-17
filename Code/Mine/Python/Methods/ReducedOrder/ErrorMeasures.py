from sympy import *
from sympy.solvers.solveset import linsolve

#
#If we redo as primitive, we can define just one polynomial, that can expand and contract
#As we go along, thus derive one expression, and then add the needed weights. 

uaj,uajm1,uajm2,uajp1,uajp2,dx = symbols('\\bar{u}_i \\bar{u}_{i-1} \\bar{u}_{i-2} \\bar{u}_{i+1} \\bar{u}_{i+2} dx')

pa,pb,pc,pd = symbols('p_a p_b p_c p_d')

x,xmxj,dx = symbols('x xmxj dx')

Poly1 = (uajp1 - uajm1)/(2*dx)*(xmxj) + uaj
Poly1Int = integrate(Poly1,xmxj)
Poly1Intjm2 = (Poly1Int.subs(xmxj,-3*dx/2) - Poly1Int.subs(xmxj,-5*dx/2)) / dx
Poly1Intjm2 = Poly1Intjm2.simplify()

# Poly1Intj = (Poly1Int.subs(xmxj,dx/2) - Poly1Int.subs(xmxj,-dx/2)) / dx
# Poly1Intj = Poly1Intj.simplify()

Poly1Intjp2 = (Poly1Int.subs(xmxj,5*dx/2) - Poly1Int.subs(xmxj,3*dx/2)) / dx
Poly1Intjp2 = Poly1Intjp2.simplify()


Poly2 = (-2*uaj + uajm1 + uajp1)/(2*dx**2)*(xmxj)**2 + (-uajm1 + uajp1)/(2*dx)*(xmxj) + 13*uaj/12 - uajm1/24 - uajp1/24
Poly2Int = integrate(Poly2,xmxj)
Poly2Intjp2 = (Poly2Int.subs(xmxj,5*dx/2) - Poly2Int.subs(xmxj,3*dx/2)) / dx
Poly2Intjp2 = Poly2Intjp2.simplify()

Poly2Intjm2 = (Poly2Int.subs(xmxj,-3*dx/2) - Poly2Int.subs(xmxj,-5*dx/2)) / dx
Poly2Intjm2 = Poly2Intjm2.simplify()



    


