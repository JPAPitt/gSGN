
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from matplotlib.pyplot import plot

def Recon(qA,x,dx,nG):
    n = len(qA)
    qcubic = zeros(4*n)
    
    # j= 0
    j = 0
    qcubic[4*j: 4*(j+nG)] = qA[j]*ones(4*nG)
    
    for j in range(nG,n-nG):

        P3r0a33,P3r0a32,P3r0a31,P3r0a30 = P3jtojp3(qA[j],qA[j+1],qA[j+2],qA[j+3],dx)
        P3r1a33,P3r1a32,P3r1a31,P3r1a30 = P3jm1tojp2(qA[j],qA[j-1],qA[j+1],qA[j+2],dx)
        P3r2a33,P3r2a32,P3r2a31,P3r2a30 = P3jm2tojp1(qA[j],qA[j-1],qA[j-2],qA[j+1],dx)
        P3r3a33,P3r3a32,P3r3a31,P3r3a30 =  P3jm3toj(qA[j],qA[j-1],qA[j-2],qA[j-3],dx)
        
        P2r0a22,P2r0a21,P2r0a20=   P2jtojp2(qA[j],qA[j+1],qA[j+2],dx)
        P2r1a22,P2r1a21,P2r1a20  =  P2jm1tojp1(qA[j-1],qA[j],qA[j+1],dx)
        P2r2a22,P2r2a21,P2r2a20  =  P2jm2toj(qA[j-2],qA[j-1],qA[j],dx)
        
        P1r0a11,P1r0a10  =  P1jtojp1(qA[j],qA[j+1],dx)
        P1r1a11,P1r1a10  =  P1jm1toj(qA[j],qA[j-1],dx)
        
        P1Wa11,P1Wa10  = P1jWeird(qA[j+1],qA[j],qA[j-1],dx)

        # w1 = minmodL([P1r0a11,P1r1a11,P1Wa11])  / P1Wa11
        # P1a11 = P1Wa11
        P1a11w1 = minmodL([P1r0a11,P1r1a11,P1Wa11])
        
        #w2 = minmodL([P2r0a22,P2r1a2,P2r2a22])  / P2r1a22
        #P2a22 = P2r1a22
        P2a22w2 = minmodL([P2r0a22,P2r1a22,P2r2a22]) 
        
        #w3 = minmodL([P3r0a33,P3r1a33,P3r2a33,P3r3a33])  / P3r2a33
        #P3a33 = P3r2a33
        P3a33w3 = minmodL([P3r0a33,P3r1a33,P3r2a33,P3r3a33])
        
        #Weight on second coefficient
        if ( abs(P3r2a33) > 10.0**(-12.0)):
            w3 = minmodL([P3r0a33,P3r1a33,P3r2a33,P3r3a33])  / P3r2a33
        else:
            w3 = 0
        #P3a32w3 = w3*P3r2a32
        P3a32w3Corr = w3*(P3r2a32 - P2a22w2) + P2a22w2
        
        
        #Weight on first coefficient, want it to reduce to Reconed P2a22
        if (abs(P2r1a22)  > 10.0**(-12.0)):
            w2 = minmodL([P2r0a22,P2r1a22,P2r2a22])  / P2r1a22
        else:
            w2 = 0
        P2a21w2 = w2 *P2r1a21
        
        P2a21w2Corr = w2*(P2r1a21 - P1a11w1) + P1a11w1
        
        P3a31w3Corr = w3*(P3r2a31 - P2a21w2Corr) + P2a21w2Corr
        
        na33 = P3a33w3
        na32 = P3a32w3Corr
 
        na31 = P3a31w3Corr
        na30 = qA[j] - na32/3*(dx/2)**2 
        

        qcubic[4*j] = na33*(-dx/2)**3 + na32*(-dx/2)**2 + na31*(-dx/2) + na30
        qcubic[4*j + 1] = na33*(-dx/6)**3 + na32*(-dx/6)**2 + na31*(-dx/6) + na30
        qcubic[4*j + 2] = na33*(dx/6)**3 + na32*(dx/6)**2 + na31*(dx/6) + na30
        qcubic[4*j + 3] = na33*(dx/2)**3 + na32*(dx/2)**2 + na31*(dx/2) + na30


    # j= n-3
    j = n-nG
    qcubic[4*j: 4*(j+nG)] = qA[j]*ones(4*nG)
        
    return qcubic


def PointsFromPolyFromPoints(a0,a1,a2,a3,x):
    return a0 + a1*x + a2*x**2 + a3*x**3


def PolyFromPoints(yjmh,yjms,yjps,yjph,dx):
    a3 = (-9*yjmh + 27*yjms - 27*yjps + 9*yjph)/ (2*dx**3)
    a2 = (9*yjmh - 9*yjms - 9*yjps + 9*yjph  )/ (4*dx**2)
    a1 = (yjmh  - 27*yjms  + 27*yjps  - yjph  )/ (8*dx)
    a0 = (-yjmh+ 9*yjms + 9*yjps  - yjph)/ 16
    return a0,a1,a2,a3

def SolveForuEdges(b11,hjmh,hjms,hjps,hjph,Gjmh,Gjms,Gjps,Gjph,dx):
    # print('Solve',beta1,ha0,ha1,ha2,ha3,Ga0,Ga1,Ga2,Ga3 ,dx)
    beta1 = b11/2.0
        
    if (abs(beta1) > 10.0**(-10)):
        Gaj = dx*(Gjmh + 3*Gjms + Gjph + 3*Gjps)/8
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
    
        ua0,ua1,ua2,ua3=PolyFromPoints(ujmh,ujms,ujps,ujph,dx)        
    else:
        hjmh = PointsFromPolyFromPoints(ha0,ha1,ha2,ha3,-dx/2)
        hjms = PointsFromPolyFromPoints(ha0,ha1,ha2,ha3,-dx/6)
        hjps = PointsFromPolyFromPoints(ha0,ha1,ha2,ha3,dx/6)
        hjph = PointsFromPolyFromPoints(ha0,ha1,ha2,ha3,dx/2)
        
        Gjmh = PointsFromPolyFromPoints(Ga0,Ga1,Ga2,Ga3,-dx/2)
        Gjms = PointsFromPolyFromPoints(Ga0,Ga1,Ga2,Ga3,-dx/6)
        Gjps = PointsFromPolyFromPoints(Ga0,Ga1,Ga2,Ga3,dx/6)
        Gjph = PointsFromPolyFromPoints(Ga0,Ga1,Ga2,Ga3,dx/2)
               
        ujmh  =  Gjmh/ hjmh  
        ujms =  Gjms/ hjms  
        ujps =  Gjps/ hjps  
        ujph   =  Gjph/ hjph      

        ua0,ua1,ua2,ua3=PolyFromPoints(ujmh,ujms,ujps,ujph,dx)
                                   
    return ua0,ua1,ua2,ua3


def CalculateG(beta1,ha0,ha1,ha2,ha3,ua0,ua1,ua2,ua3 ,dx):

        hjmh = PointsFromPolyFromPoints(ha0,ha1,ha2,ha3,-dx/2)
        hjms = PointsFromPolyFromPoints(ha0,ha1,ha2,ha3,-dx/6)
        hjps = PointsFromPolyFromPoints(ha0,ha1,ha2,ha3,dx/6)
        hjph = PointsFromPolyFromPoints(ha0,ha1,ha2,ha3,dx/2)

        dhjmh = PointsFromPolyFromPoints(ha1,2*ha2,3*ha3,0,-dx/2)
        dhjms = PointsFromPolyFromPoints(ha1,2*ha2,3*ha3,0,-dx/6)
        dhjps = PointsFromPolyFromPoints(ha1,2*ha2,3*ha3,0,dx/6)
        dhjph = PointsFromPolyFromPoints(ha1,2*ha2,3*ha3,0,dx/2)
        
        ujmh = PointsFromPolyFromPoints(ua0,ua1,ua2,ua3,-dx/2)
        ujms = PointsFromPolyFromPoints(ua0,ua1,ua2,ua3,-dx/6)
        ujps = PointsFromPolyFromPoints(ua0,ua1,ua2,ua3,dx/6)
        ujph = PointsFromPolyFromPoints(ua0,ua1,ua2,ua3,dx/2)

        dujmh = PointsFromPolyFromPoints(ua1,2*ua2,3*ua3,0,-dx/2)
        dujms = PointsFromPolyFromPoints(ua1,2*ua2,3*ua3,0,-dx/6)
        dujps = PointsFromPolyFromPoints(ua1,2*ua2,3*ua3,0,dx/6)
        dujph = PointsFromPolyFromPoints(ua1,2*ua2,3*ua3,0,dx/2)

        ddujmh = PointsFromPolyFromPoints(2*ua2,6*ua3,0,0,-dx/2)
        ddujms = PointsFromPolyFromPoints(2*ua2,6*ua3,0,0,-dx/6)
        ddujps = PointsFromPolyFromPoints(2*ua2,6*ua3,0,0,dx/6)
        ddujph = PointsFromPolyFromPoints(2*ua2,6*ua3,0,0,dx/2)
               
        Gjmh = ujmh*hjmh - b1/2.0*(3*hjmh**2*dhjmh*dujmh + hjmh**3*ddujmh )  
        Gjms = ujms*hjms - b1/2.0*(3*hjms**2*dhjms*dujms + hjms**3*ddujms ) 
        Gjps = ujps*hjps - b1/2.0*(3*hjps**2*dhjps*dujps + hjps**3*ddujps ) 
        Gjph = ujph*hjph - b1/2.0*(3*hjph**2*dhjph*dujph + hjph**3*ddujph ) 

                                   
        return  Gjmh, Gjms, Gjps,  Gjph

def InitialSoliton(x,dx,t,ga,a0,a1,b1):
    n = len(x)
    hp = []
    Gp = []
    up = []
    xp = []
    
    k = sqrt(3*a1) / (2*a0*sqrt(a0 + a1))
    c = sqrt(ga*(a0 + a1))
    for i in range(n):
        
        
        xjmh = x[i] - dx/2
        xjms = x[i] - dx/6
        xjps = x[i] + dx/6
        xjph = x[i] + dx/2
        
        #xjmh
        phi = xjmh - c*t
        sechkphi = 1.0 / cosh(k*phi)
        hjmh = a0 + a1*sechkphi**2
        ujmh = c*(1 - a0 /  hjmh)
        Gjmh = ujmh*hjmh + b1*a0*a1*c*(k**2)*sechkphi**4*hjmh - 2*b1*a0*a1**2*c*k**2*sechkphi**4*tanh(k*phi)**2  - 2*b1*a0*a1*c*k**2*sechkphi**2*hjmh*tanh(k*phi)**2 

        phi = xjms - c*t
        sechkphi = 1.0 / cosh(k*phi)
        hjms = a0 + a1*sechkphi**2
        ujms = c*(1 - a0 /  hjms)
        Gjms = ujms*hjms + b1*a0*a1*c*(k**2)*sechkphi**4*hjms - 2*b1*a0*a1**2*c*k**2*sechkphi**4*tanh(k*phi)**2  - 2*b1*a0*a1*c*k**2*sechkphi**2*hjms*tanh(k*phi)**2       
        

        phi = xjps - c*t
        sechkphi = 1.0 / cosh(k*phi)
        hjps = a0 + a1*sechkphi**2
        ujps = c*(1 - a0 /  hjps)
        Gjps = ujps*hjps + b1*a0*a1*c*(k**2)*sechkphi**4*hjps - 2*b1*a0*a1**2*c*k**2*sechkphi**4*tanh(k*phi)**2  - 2*b1*a0*a1*c*k**2*sechkphi**2*hjps*tanh(k*phi)**2  
 
        phi = xjph - c*t
        sechkphi = 1.0 / cosh(k*phi)
        hjph = a0 + a1*sechkphi**2
        ujph = c*(1 - a0 /  hjph)
        Gjph = ujph*hjph + b1*a0*a1*c*(k**2)*sechkphi**4*hjph - 2*b1*a0*a1**2*c*k**2*sechkphi**4*tanh(k*phi)**2  - 2*b1*a0*a1*c*k**2*sechkphi**2*hjph*tanh(k*phi)**2  
       

        hcell = [hjmh,hjms,hjps,hjph]
        ucell = [ujmh,ujms,ujps,ujph]
        Gcell = [Gjmh,Gjms,Gjps,Gjph]
        xhcell = [xjmh,xjms,xjps,xjph]
        hp = hp + hcell
        xp = xp + xhcell
        up = up + ucell
        Gp = Gp + Gcell
        
        
    return array(hp),array(Gp),array(up),array(xp)

def FEM(beta1,h,G,n,x,dx):
    up = []
    for i in range(n):
        # ha0,ha1,ha2,ha3 = PolyFromPoints(h[4*i],h[4*i+1],h[4*i+2],h[4*i+3],dx)
        # Ga0,Ga1,Ga2,Ga3 = PolyFromPoints(G[4*i],G[4*i+1],G[4*i+2],G[4*i+3],dx)
        # print('Initial h', ha0,ha1,ha2,ha3)
        # print('Initial G', Ga0,Ga1,Ga2,Ga3)
        ua0,ua1,ua2,ua3 = SolveForuEdges(beta1,h[4*i],h[4*i+1],h[4*i+2],h[4*i+3],G[4*i],G[4*i+1],G[4*i+2],G[4*i+3] ,dx)
        print('U Coeff', ua0,ua1,ua2,ua3)
        ujmh = PointsFromPolyFromPoints(ua0,ua1,ua2,ua3,-dx/2)
        ujms = PointsFromPolyFromPoints(ua0,ua1,ua2,ua3,-dx/6)
        ujps = PointsFromPolyFromPoints(ua0,ua1,ua2,ua3,dx/6)
        ujph = PointsFromPolyFromPoints(ua0,ua1,ua2,ua3,dx/2)
        print('U Value',ujmh,ujms,ujps,ujph)
        ucell = [ujmh,ujms,ujps,ujph]
        up = up + ucell
    
    return array(up)
    
def FEMforG(beta1,h,u,n,x,dx):
    
    Gp = []
    for i in range(n):
        ha0,ha1,ha2,ha3 = PolyFromPoints(h[4*i],h[4*i+1],h[4*i+2],h[4*i+3],dx)
        ua0,ua1,ua2,ua3 = PolyFromPoints(u[4*i],u[4*i+1],u[4*i+2],u[4*i+3],dx)
        Gjmh, Gjms, Gjps,  Gjph = CalculateG(b1,ha0,ha1,ha2,ha3,ua0,ua1,ua2,ua3 ,dx)
        Gcell = [Gjmh, Gjms, Gjps,  Gjph]
        
        Gp = Gp + Gcell 
    return Gp

sx = -1
ex = 1
n= 100
dx = (ex- sx)/n
xn = arange(sx,ex+ dx,dx)

ga = 10.0
a0 = 1.0
a1 = 0.7
t = 0
b1 = 2.0/3.0

h,G,u,x = InitialSoliton(xn,dx,t,ga,a0,a1,b1)

# Gn = FEMforG(b1,h,u,n+1,x,dx)

un = FEM(b1,h,G,len(xn),x,dx)

# i = n//2
# ha0,ha1,ha2,ha3 = PolyFromPoints(h[4*i],h[4*i+1],h[4*i+2],h[4*i+3],dx)
# ua0,ua1,ua2,ua3 = PolyFromPoints(u[4*i],u[4*i+1],u[4*i+2],u[4*i+3],dx)
# Gjmh, Gjms, Gjps,  Gjph = CalculateG(b1,ha0,ha1,ha2,ha3,ua0,ua1,ua2,ua3 ,dx)
# print(Gjmh, Gjms, Gjps,  Gjph)
# print(G[4*i:4*(i+1)])
# G = FEMforG(xn,h,u,dx)