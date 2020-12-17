from sympy import *
from sympy.solvers.solveset import linsolve

#
#If we redo as primitive, we can define just one polynomial, that can expand and contract
#As we go along, thus derive one expression, and then add the needed weights. 

k = 4

uaj,uajm1,uajm2,uajp1,uajp2,dx = symbols('\\bar{u}_i \\bar{u}_{i-1} \\bar{u}_{i-2} \\bar{u}_{i+1} \\bar{u}_{i+2} dx')

pa,pb,pc = symbols('p_a p_b p_c')

x,xmxj,xaj,xj,dx = symbols('x xmxj x^*_j x_j dx')

xm5jh,xm3jh,xmjh,xpjh,xp3jh,xp5jh = symbols('x_{j-5/2} x_{j-3/2} x_{j-1/2} x_{j+1/2} x_{j+3/2} x_{j+5/2}')

ual = [uajm2,uajm1,uajp1,uajp2]
xl = [xm5jh,xm3jh,xpjh,xp3jh,xp5jh]

def GenSum(k,m,r,ci,x,xl,dx):
    Sum = 0
    for l in range(k+1):
        if l != m:
            Prod = 1
            for q in range(k+1):
                if q!=l and q!=m:
                    Prod = Prod*(x - xl[ci -r + q])
            Sum = Sum + Prod
    Norm = 1
    for l in range(k+1):
        if l!=m:
             Norm =  Norm*(m-l)*dx
            
    return Sum/Norm
    

def GeneratePoly(k,r,ci,ual,cix,x,xl,dx):
    Poly = 0
    for m in range(k+1):
        for j in range(m):
            Poly = Poly +ual[ci+j-r]*dx*GenSum(k,m,r,cix,x,xl,dx)

    return Poly
            
k = 4
ci = 2
cix = 2

Polyo4r0 = GeneratePoly(k,2,ci,ual,cix,x,xl,dx)
Polyo4r0N =  Polyo4r0.subs(xm5jh,xj - 5*dx/2).subs(xm3jh,xj - 3*dx/2).subs(xp3jh,xj + 3*dx/2) \
        .subs(xpjh,xj + dx/2).subs(xp5jh,xj + 5*dx/2)
Polyo4r0N = Polyo4r0N.subs(x - xj, xmxj)
Polyo4r0N = Polyo4r0N.simplify()
Polyo4r0N =  collect(Polyo4r0N,xmxj)


Polyo4r0NCoeff0  = Polyo4r0N.subs(xmxj,0)
Polyo4r0NCoeff0 = Polyo4r0NCoeff0.simplify()

Polyo4r0NCoeff1 = Polyo4r0N - Polyo4r0NCoeff0
Polyo4r0NCoeff1 = Polyo4r0NCoeff1.subs(xmxj**3,0).subs(xmxj**2,0)
Polyo4r0NCoeff1 = Polyo4r0NCoeff1.simplify()

Polyo4r0NCoeff2 = Polyo4r0N - Polyo4r0NCoeff0 - Polyo4r0NCoeff1
Polyo4r0NCoeff2 = Polyo4r0NCoeff2.subs(xmxj**3,0)
Polyo4r0NCoeff2 = Polyo4r0NCoeff2.simplify()

Polyo4r0NCoeff3 = Polyo4r0N - Polyo4r0NCoeff0 - Polyo4r0NCoeff1 - Polyo4r0NCoeff2
Polyo4r0NCoeff3 = Polyo4r0NCoeff3.simplify()

Polyo4full = Polyo4r0NCoeff3 + Polyo4r0NCoeff2 + Polyo4r0NCoeff1 + Polyo4r0NCoeff0
Polyo4fullInt = integrate(Polyo4full,xmxj)

Polyo4fullIntm2 = (Polyo4fullInt.subs(xmxj,-3*dx/2) - Polyo4fullInt.subs(xmxj,-5*dx/2)) / dx
Polyo4fullIntm2 = Polyo4fullIntm2.simplify()
Polyo4fullIntm1 = (Polyo4fullInt.subs(xmxj,-dx/2) - Polyo4fullInt.subs(xmxj,-3*dx/2)) / dx
Polyo4fullIntm1 = Polyo4fullIntm1.simplify()
Polyo4fullInt = (Polyo4fullInt.subs(xmxj,dx/2) - Polyo4fullInt.subs(xmxj,-dx/2)) / dx
Polyo4fullInt = Polyo4fullInt.simplify()
Polyo4fullIntp1 = (Polyo4fullInt.subs(xmxj,3*dx/2) - Polyo4fullInt.subs(xmxj,dx/2)) / dx
Polyo4fullIntp1 = Polyo4fullIntp1.simplify()
Polyo4fullIntp2 = (Polyo4fullInt.subs(xmxj,5*dx/2) - Polyo4fullInt.subs(xmxj,3*dx/2)) / dx
Polyo4fullIntp2 = Polyo4fullIntp2.simplify()

# Polyo4r0NCoeff3  = Polyo4r0N.coeff(xmxj,3)
# Polyo4r0NCoeff2  = Polyo4r0N.coeff(xmxj,2)
# Polyo4r0NCoeff1  = Polyo4r0N.coeff(xmxj,1)
# Polyo4r0NCoeff0  = Polyo4r0N.coeff(xmxj,0)
# Polyo3r1 = GeneratePoly(k,1,ci,ual,cix,x,xl,dx)
# Polyo3r2 = GeneratePoly(k,2,ci,ual,cix,x,xl,dx)

# Polyo2r0 = GeneratePoly(2,0,ci,ual,cix,x,xl,dx)
# Polyo2r1 = GeneratePoly(2,1,ci,ual,cix,x,xl,dx)

# Polyo1r0 = GeneratePoly(1,0,ci,ual,cix,x,xl,dx)
      
    


