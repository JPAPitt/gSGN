"""
Python code to generate csv files that contain integrals of the basis functions
needed to generate the matrices of the Finite Element Method.

by Jordan Pitt 11/10/2018
"""

"""
############################# IMPORTS #########################################
"""


from sympy import * #Related website: https://www.sympy.org/en/index.html
from IPython.display import display

# symbols
x,t,dx,dt,ga,b1 = symbols('x,t,dx,dt,g,beta1', positive = True, nonzero = True,real=True)

def PolyFromPoints(yjmh,yjms,yjps,yjph,dx):
    a3 = (-9*yjmh + 27*yjms - 27*yjps + 9*yjph)/ (2*dx**3)
    a2 = (9*yjmh - 9*yjms - 9*yjps + 9*yjph  )/ (4*dx**2)
    a1 = (yjmh  - 27*yjms  + 27*yjps  - yjph  )/ (8*dx)
    a0 = (-yjmh+ 9*yjms + 9*yjps  - yjph)/ 16
    return a0,a1,a2,a3
    
def SolveForuEdges(beta1,hjmh,hjms,hjps,hjph,Gjmh,Gjms,Gjps,Gjph,Gaj,dx):
    
    if (abs(beta1) > 10.0**(-10)):
        Gj = -Gjmh/16 + 9*Gjms/16 - Gjph/16 + 9*Gjps/16
        hj = -hjmh/16 + 9*hjms/16 - hjph/16 + 9*hjps/16
        dhjmh = (-11*hjmh + 18*hjms + 2*hjph - 9*hjps)/(2*dx)
        dhj = (hjmh - 27*hjms - hjph + 27*hjps)/(8*dx)
        dhjph = (-2*hjmh + 9*hjms + 11*hjph - 18*hjps)/(2*dx)
        
        print(dhjmh,dhj,dhjph)
        
        a0 = -11*beta1*hjmh**3/(2*dx**2) + beta1*hjph**3/dx**2 + hjmh/8
        a1 = 9*beta1*hjmh**3/dx**2 - 9*beta1*hjph**3/(2*dx**2) + 3*hjms/8
        a2 = -9*beta1*hjmh**3/(2*dx**2) + 9*beta1*hjph**3/dx**2 + 3*hjps/8
        a3 = beta1*hjmh**3/dx**2 - 11*beta1*hjph**3/(2*dx**2) + hjph/8
        
        
        b0 = 33*beta1*dhjmh*hjmh**2/(2*dx) - 18*beta1*hjmh**3/dx**2 + hjmh
        b1 = -27*beta1*dhjmh*hjmh**2/dx + 45*beta1*hjmh**3/dx**2
        b2 = 27*beta1*dhjmh*hjmh**2/(2*dx) - 36*beta1*hjmh**3/dx**2
        b3 = -3*beta1*dhjmh*hjmh**2/dx + 9*beta1*hjmh**3/dx**2
        
        c0 = 3*beta1*dhjph*hjph**2/dx + 9*beta1*hjph**3/dx**2
        c1 = -27*beta1*dhjph*hjph**2/(2*dx) - 36*beta1*hjph**3/dx**2
        c2 = 27*beta1*dhjph*hjph**2/dx + 45*beta1*hjph**3/dx**2
        c3 = -33*beta1*dhjph*hjph**2/(2*dx) - 18*beta1*hjph**3/dx**2 + hjph
    
    
        d0 = -3*beta1*dhj*hj**2/(8*dx) - 9*beta1*hj**3/(2*dx**2) - hj/16
        d1 = 81*beta1*dhj*hj**2/(8*dx) + 9*beta1*hj**3/(2*dx**2) + 9*hj/16
        d2 = -81*beta1*dhj*hj**2/(8*dx) + 9*beta1*hj**3/(2*dx**2) + 9*hj/16
        d3 = 3*beta1*dhj*hj**2/(8*dx) - 9*beta1*hj**3/(2*dx**2) - hj/16
    
        print(a0,a1,a2,a3)
        print(b0,b1,b2,b3)
        print(c0,c1,c2,c3)
        print(d0,d1,d2,d3)

        Div = (a0*b1*c2*d3 - a0*b1*c3*d2 - a0*b2*c1*d3 + a0*b2*c3*d1 \
               + a0*b3*c1*d2 - a0*b3*c2*d1 - a1*b0*c2*d3 + a1*b0*c3*d2 + a1*b2*c0*d3 - a1*b2*c3*d0 \
               - a1*b3*c0*d2 + a1*b3*c2*d0 + a2*b0*c1*d3 - a2*b0*c3*d1 - a2*b1*c0*d3 + a2*b1*c3*d0 \
               + a2*b3*c0*d1 - a2*b3*c1*d0 - a3*b0*c1*d2 + a3*b0*c2*d1 + a3*b1*c0*d2 - a3*b1*c2*d0 \
               - a3*b2*c0*d1 + a3*b2*c1*d0)
        print(Div)
        if abs(Div) < 10.0**(-15):
            ujph = Gjph/hjph 
        else:
            ujph = (-Gaj*b0*c1*d2 + Gaj*b0*c2*d1 + Gaj*b1*c0*d2 - Gaj*b1*c2*d0 - Gaj*b2*c0*d1 + Gaj*b2*c1*d0 \
               + Gj*a0*b1*c2 - Gj*a0*b2*c1 - Gj*a1*b0*c2 + Gj*a1*b2*c0 + Gj*a2*b0*c1 - Gj*a2*b1*c0 \
               + Gjmh*a0*c1*d2 - Gjmh*a0*c2*d1 - Gjmh*a1*c0*d2 + Gjmh*a1*c2*d0 + Gjmh*a2*c0*d1 \
               - Gjmh*a2*c1*d0 - Gjph*a0*b1*d2 + Gjph*a0*b2*d1 + Gjph*a1*b0*d2 - Gjph*a1*b2*d0 \
               - Gjph*a2*b0*d1 + Gjph*a2*b1*d0)/Div
                                                 
        Div = (a0*b1*c2*d3 \
               - a0*b1*c3*d2 - a0*b2*c1*d3 + a0*b2*c3*d1 + a0*b3*c1*d2 - a0*b3*c2*d1 \
               - a1*b0*c2*d3 + a1*b0*c3*d2 + a1*b2*c0*d3 - a1*b2*c3*d0 - a1*b3*c0*d2 \
               + a1*b3*c2*d0 + a2*b0*c1*d3 - a2*b0*c3*d1 - a2*b1*c0*d3 + a2*b1*c3*d0 \
               + a2*b3*c0*d1 - a2*b3*c1*d0 - a3*b0*c1*d2 + a3*b0*c2*d1 + a3*b1*c0*d2 \
               - a3*b1*c2*d0 - a3*b2*c0*d1 + a3*b2*c1*d0)     
        print(Div)
        if abs(Div) < 10.0**(-15):
            ujps = Gjps/hjps
        else:                                     
            ujps = (Gaj*b0*c1*d3 - Gaj*b0*c3*d1 - Gaj*b1*c0*d3 + Gaj*b1*c3*d0 + Gaj*b3*c0*d1 \
                   - Gaj*b3*c1*d0 - Gj*a0*b1*c3 + Gj*a0*b3*c1 + Gj*a1*b0*c3 - Gj*a1*b3*c0 \
                   - Gj*a3*b0*c1 + Gj*a3*b1*c0 - Gjmh*a0*c1*d3 + Gjmh*a0*c3*d1 + Gjmh*a1*c0*d3 \
                   - Gjmh*a1*c3*d0 - Gjmh*a3*c0*d1 + Gjmh*a3*c1*d0 + Gjph*a0*b1*d3 - Gjph*a0*b3*d1 \
                   - Gjph*a1*b0*d3 + Gjph*a1*b3*d0 + Gjph*a3*b0*d1 - Gjph*a3*b1*d0)/Div


        Div = (a0*b1*c2*d3 \
                - a0*b1*c3*d2 - a0*b2*c1*d3 + a0*b2*c3*d1 + a0*b3*c1*d2 \
                - a0*b3*c2*d1 - a1*b0*c2*d3 + a1*b0*c3*d2 + a1*b2*c0*d3 \
                - a1*b2*c3*d0 - a1*b3*c0*d2 + a1*b3*c2*d0 + a2*b0*c1*d3 \
                - a2*b0*c3*d1 - a2*b1*c0*d3 + a2*b1*c3*d0 + a2*b3*c0*d1 \
                - a2*b3*c1*d0 - a3*b0*c1*d2 + a3*b0*c2*d1 + a3*b1*c0*d2 \
                - a3*b1*c2*d0 - a3*b2*c0*d1 + a3*b2*c1*d0)   
        print(Div)
        if abs(Div) < 10.0**(-15):
            ujms = Gjms/hjms
        else:                                     
            ujms = (-Gaj*b0*c2*d3 + Gaj*b0*c3*d2 + Gaj*b2*c0*d3 - Gaj*b2*c3*d0 \
                - Gaj*b3*c0*d2 + Gaj*b3*c2*d0 + Gj*a0*b2*c3 - Gj*a0*b3*c2 \
                - Gj*a2*b0*c3 + Gj*a2*b3*c0 + Gj*a3*b0*c2 - Gj*a3*b2*c0 + Gjmh*a0*c2*d3 \
                - Gjmh*a0*c3*d2 - Gjmh*a2*c0*d3 + Gjmh*a2*c3*d0 + Gjmh*a3*c0*d2 \
                - Gjmh*a3*c2*d0 - Gjph*a0*b2*d3 + Gjph*a0*b3*d2 + Gjph*a2*b0*d3 \
                - Gjph*a2*b3*d0 - Gjph*a3*b0*d2 + Gjph*a3*b2*d0) / Div
       
        Div = (a0*b1*c2*d3 - a0*b1*c3*d2 - a0*b2*c1*d3 \
                + a0*b2*c3*d1 + a0*b3*c1*d2 - a0*b3*c2*d1 - a1*b0*c2*d3 + a1*b0*c3*d2 \
                + a1*b2*c0*d3 - a1*b2*c3*d0 - a1*b3*c0*d2 + a1*b3*c2*d0 + a2*b0*c1*d3 \
                - a2*b0*c3*d1 - a2*b1*c0*d3 + a2*b1*c3*d0 + a2*b3*c0*d1 - a2*b3*c1*d0 \
                - a3*b0*c1*d2 + a3*b0*c2*d1 + a3*b1*c0*d2 - a3*b1*c2*d0 - a3*b2*c0*d1 + a3*b2*c1*d0)    
        print(Div)
        if abs(Div) < 10.0**(-15):
            ujmh = Gjmh/hjmh
        else:                 
            ujmh = (Gaj*b1*c2*d3 - Gaj*b1*c3*d2 - Gaj*b2*c1*d3 + Gaj*b2*c3*d1 \
                + Gaj*b3*c1*d2 - Gaj*b3*c2*d1 - Gj*a1*b2*c3 + Gj*a1*b3*c2 + Gj*a2*b1*c3 \
                - Gj*a2*b3*c1 - Gj*a3*b1*c2 + Gj*a3*b2*c1 - Gjmh*a1*c2*d3 + Gjmh*a1*c3*d2\
                + Gjmh*a2*c1*d3 - Gjmh*a2*c3*d1 - Gjmh*a3*c1*d2 + Gjmh*a3*c2*d1 \
                + Gjph*a1*b2*d3 - Gjph*a1*b3*d2 - Gjph*a2*b1*d3 + Gjph*a2*b3*d1 \
                + Gjph*a3*b1*d2 - Gjph*a3*b2*d1)/Div
    else:
        ujph = Gjph/hjph
        ujps = Gjps/hjps
        ujms = Gjms/hjms
        ujmh = Gjmh/hjmh
                                                  
    return ujmh,ujms,ujps, ujph


def SolveForuNodes(beta1,hjmh,hjms,hjps,hjph,Gjmh,Gjms,Gjps,Gjph,Gaj,dx):
    
    if (abs(beta1) > 10.0**(-10)):
        Gj = -Gjmh/16 + 9*Gjms/16 - Gjph/16 + 9*Gjps/16
        hj = -hjmh/16 + 9*hjms/16 - hjph/16 + 9*hjps/16
        dhjmh = (-11*hjmh + 18*hjms + 2*hjph - 9*hjps)/(2*dx)
        dhj = (hjmh - 27*hjms - hjph + 27*hjps)/(8*dx)
        dhjph = (-2*hjmh + 9*hjms + 11*hjph - 18*hjps)/(2*dx)
        
        a0 = -11*beta1*hjmh**3/(2*dx) + beta1*hjph**3/dx + dx*hjmh/4
        a1 = 9*beta1*hjmh**3/dx - 9*beta1*hjph**3/(2*dx) + 3*dx*hjms/8
        a2 = -9*beta1*hjmh**3/(2*dx) + 9*beta1*hjph**3/dx + 3*dx*hjps/8
        a3 = beta1*hjmh**3/dx - 11*beta1*hjph**3/(2*dx)
        
        b0 = 33*beta1*dhjmh*hjmh**2/(2*dx) - 18*beta1*hjmh**3/dx**2 + hjmh
        b1 = -27*beta1*dhjmh*hjmh**2/dx + 45*beta1*hjmh**3/dx**2
        b2 = 27*beta1*dhjmh*hjmh**2/(2*dx) - 36*beta1*hjmh**3/dx**2
        b3 = -3*beta1*dhjmh*hjmh**2/dx + 9*beta1*hjmh**3/dx**2
        
        c0 = 3*beta1*dhjph*hjph**2/dx + 9*beta1*hjph**3/dx**2
        c1 = -27*beta1*dhjph*hjph**2/(2*dx) - 36*beta1*hjph**3/dx**2
        c2 = 27*beta1*dhjph*hjph**2/dx + 45*beta1*hjph**3/dx**2
        c3 = -33*beta1*dhjph*hjph**2/(2*dx) - 18*beta1*hjph**3/dx**2 + hjph
    
    
        d0 = -3*beta1*dhj*hj**2/(8*dx) - 9*beta1*hj**3/(2*dx**2) - hj/16
        d1 = 81*beta1*dhj*hj**2/(8*dx) + 9*beta1*hj**3/(2*dx**2) + 9*hj/16
        d2 = -81*beta1*dhj*hj**2/(8*dx) + 9*beta1*hj**3/(2*dx**2) + 9*hj/16
        d3 = 3*beta1*dhj*hj**2/(8*dx) - 9*beta1*hj**3/(2*dx**2) - hj/16
    
        print(a0,a1,a2,a3)
        print(b0,b1,b2,b3)
        print(c0,c1,c2,c3)
        print(d0,d1,d2,d3)
                
        ujph = (-Gjmh*b0*c1*d2 + Gjmh*b0*c2*d1 + Gjmh*b1*c0*d2 - Gjmh*b1*c2*d0 - Gjmh*b2*c0*d1 + Gjmh*b2*c1*d0 + Gjms*a0*c1*d2 - Gjms*a0*c2*d1 - Gjms*a1*c0*d2 + Gjms*a1*c2*d0 + Gjms*a2*c0*d1 - Gjms*a2*c1*d0 + Gjph*a0*b1*c2 - Gjph*a0*b2*c1 - Gjph*a1*b0*c2 + Gjph*a1*b2*c0 + Gjph*a2*b0*c1 - Gjph*a2*b1*c0 - Gjps*a0*b1*d2 + Gjps*a0*b2*d1 + Gjps*a1*b0*d2 - Gjps*a1*b2*d0 - Gjps*a2*b0*d1 + Gjps*a2*b1*d0)/(a0*b1*c2*d3 - a0*b1*c3*d2 - a0*b2*c1*d3 + a0*b2*c3*d1 + a0*b3*c1*d2 - a0*b3*c2*d1 - a1*b0*c2*d3 + a1*b0*c3*d2 + a1*b2*c0*d3 - a1*b2*c3*d0 - a1*b3*c0*d2 + a1*b3*c2*d0 + a2*b0*c1*d3 - a2*b0*c3*d1 - a2*b1*c0*d3 + a2*b1*c3*d0 + a2*b3*c0*d1 - a2*b3*c1*d0 - a3*b0*c1*d2 + a3*b0*c2*d1 + a3*b1*c0*d2 - a3*b1*c2*d0 - a3*b2*c0*d1 + a3*b2*c1*d0)
        ujps = (Gjmh*b0*c1*d3 - Gjmh*b0*c3*d1 - Gjmh*b1*c0*d3 + Gjmh*b1*c3*d0 + Gjmh*b3*c0*d1 - Gjmh*b3*c1*d0 - Gjms*a0*c1*d3 + Gjms*a0*c3*d1 + Gjms*a1*c0*d3 - Gjms*a1*c3*d0 - Gjms*a3*c0*d1 + Gjms*a3*c1*d0 - Gjph*a0*b1*c3 + Gjph*a0*b3*c1 + Gjph*a1*b0*c3 - Gjph*a1*b3*c0 - Gjph*a3*b0*c1 + Gjph*a3*b1*c0 + Gjps*a0*b1*d3 - Gjps*a0*b3*d1 - Gjps*a1*b0*d3 + Gjps*a1*b3*d0 + Gjps*a3*b0*d1 - Gjps*a3*b1*d0)/(a0*b1*c2*d3 - a0*b1*c3*d2 - a0*b2*c1*d3 + a0*b2*c3*d1 + a0*b3*c1*d2 - a0*b3*c2*d1 - a1*b0*c2*d3 + a1*b0*c3*d2 + a1*b2*c0*d3 - a1*b2*c3*d0 - a1*b3*c0*d2 + a1*b3*c2*d0 + a2*b0*c1*d3 - a2*b0*c3*d1 - a2*b1*c0*d3 + a2*b1*c3*d0 + a2*b3*c0*d1 - a2*b3*c1*d0 - a3*b0*c1*d2 + a3*b0*c2*d1 + a3*b1*c0*d2 - a3*b1*c2*d0 - a3*b2*c0*d1 + a3*b2*c1*d0)
        ujms = (-Gjmh*b0*c2*d3 + Gjmh*b0*c3*d2 + Gjmh*b2*c0*d3 - Gjmh*b2*c3*d0 - Gjmh*b3*c0*d2 + Gjmh*b3*c2*d0 + Gjms*a0*c2*d3 - Gjms*a0*c3*d2 - Gjms*a2*c0*d3 + Gjms*a2*c3*d0 + Gjms*a3*c0*d2 - Gjms*a3*c2*d0 + Gjph*a0*b2*c3 - Gjph*a0*b3*c2 - Gjph*a2*b0*c3 + Gjph*a2*b3*c0 + Gjph*a3*b0*c2 - Gjph*a3*b2*c0 - Gjps*a0*b2*d3 + Gjps*a0*b3*d2 + Gjps*a2*b0*d3 - Gjps*a2*b3*d0 - Gjps*a3*b0*d2 + Gjps*a3*b2*d0)/(a0*b1*c2*d3 - a0*b1*c3*d2 - a0*b2*c1*d3 + a0*b2*c3*d1 + a0*b3*c1*d2 - a0*b3*c2*d1 - a1*b0*c2*d3 + a1*b0*c3*d2 + a1*b2*c0*d3 - a1*b2*c3*d0 - a1*b3*c0*d2 + a1*b3*c2*d0 + a2*b0*c1*d3 - a2*b0*c3*d1 - a2*b1*c0*d3 + a2*b1*c3*d0 + a2*b3*c0*d1 - a2*b3*c1*d0 - a3*b0*c1*d2 + a3*b0*c2*d1 + a3*b1*c0*d2 - a3*b1*c2*d0 - a3*b2*c0*d1 + a3*b2*c1*d0)
        ujmh = (Gjmh*b1*c2*d3 - Gjmh*b1*c3*d2 - Gjmh*b2*c1*d3 + Gjmh*b2*c3*d1 + Gjmh*b3*c1*d2 - Gjmh*b3*c2*d1 - Gjms*a1*c2*d3 + Gjms*a1*c3*d2 + Gjms*a2*c1*d3 - Gjms*a2*c3*d1 - Gjms*a3*c1*d2 + Gjms*a3*c2*d1 - Gjph*a1*b2*c3 + Gjph*a1*b3*c2 + Gjph*a2*b1*c3 - Gjph*a2*b3*c1 - Gjph*a3*b1*c2 + Gjph*a3*b2*c1 + Gjps*a1*b2*d3 - Gjps*a1*b3*d2 - Gjps*a2*b1*d3 + Gjps*a2*b3*d1 + Gjps*a3*b1*d2 - Gjps*a3*b2*d1)/(a0*b1*c2*d3 - a0*b1*c3*d2 - a0*b2*c1*d3 + a0*b2*c3*d1 + a0*b3*c1*d2 - a0*b3*c2*d1 - a1*b0*c2*d3 + a1*b0*c3*d2 + a1*b2*c0*d3 - a1*b2*c3*d0 - a1*b3*c0*d2 + a1*b3*c2*d0 + a2*b0*c1*d3 - a2*b0*c3*d1 - a2*b1*c0*d3 + a2*b1*c3*d0 + a2*b3*c0*d1 - a2*b3*c1*d0 - a3*b0*c1*d2 + a3*b0*c2*d1 + a3*b1*c0*d2 - a3*b1*c2*d0 - a3*b2*c0*d1 + a3*b2*c1*d0)
    else:
        ujph = Gjph/hjph
        ujps = Gjps/hjps
        ujms = Gjms/hjms
        ujmh = Gjmh/hjmh
                                                  
    return ujmh,ujms,ujps, ujph

dx = 0.00001
h0 = 2
G0 = 1.5

print('A')
A = SolveForuEdges(2.0/3.0,h0,h0,h0,h0,G0,G0,G0,G0,G0,dx)   
print('B')
B = SolveForuEdges(2.0/3.0,h0+0.1,h0-0.1,h0+0.1,h0-0.1,G0+0.1,G0-0.1,G0+0.1,G0-0.1,G0,dx)   
print('C')
C = SolveForuEdges(2.0/3.0,h0+0.1,h0-0.1,h0+0.1,h0-0.1,G0,G0,G0,G0,G0,dx)   
uCPoly = PolyFromPoints(C[0],C[1],C[2],C[3],dx)
print('D')    
D = SolveForuEdges(2.0/3.0,h0,h0,h0,h0,G0+0.1,G0-0.1,G0+0.1,G0-0.1,G0,dx)   

#jpLets Check                                                     