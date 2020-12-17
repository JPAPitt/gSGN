
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

# def SolveForuEdges(beta1,ha0,ha1,ha2,ha3,Ga0,Ga1,Ga2,Ga3 ,dx):
#     print('Solve',beta1,ha0,ha1,ha2,ha3,Ga0,Ga1,Ga2,Ga3 ,dx)
#     if (abs(beta1) > 10.0**(-10)):
        
#         ua3 = 2.0*(54.0*Ga0*beta1**2*ha0**4*ha1*ha3**2 - 18.0*Ga0*beta1**2*ha0**4*ha2**2*ha3 + 288.0*Ga0*beta1**2*ha0**3*ha1**2*ha2*ha3 - 180.0*Ga0*beta1**2*ha0**3*ha1*ha2**3 - 108.0*Ga0*beta1**2*ha0**2*ha1**4*ha3 + 120.0*Ga0*beta1**2*ha0**2*ha1**3*ha2**2 + 12.0*Ga0*beta1**2*ha0*ha1**5*ha2 + 12.0*Ga0*beta1**2*ha1**7 - 21.0*Ga0*beta1*ha0**2*ha1**2*ha3 + 72.0*Ga0*beta1*ha0**2*ha1*ha2**2 - 22.0*Ga0*beta1*ha0*ha1**3*ha2 - 11.0*Ga0*beta1*ha1**5 + 2.0*Ga0*ha0**2*ha3 - 4.0*Ga0*ha0*ha1*ha2 + 2.0*Ga0*ha1**3 - 99.0*Ga1*beta1**2*ha0**5*ha3**2 - 213.0*Ga1*beta1**2*ha0**4*ha1*ha2*ha3 + 96.0*Ga1*beta1**2*ha0**4*ha2**3 + 120.0*Ga1*beta1**2*ha0**3*ha1**3*ha3 - 156.0*Ga1*beta1**2*ha0**3*ha1**2*ha2**2 - 24.0*Ga1*beta1**2*ha0**2*ha1**4*ha2 - 12.0*Ga1*beta1**2*ha0*ha1**6 + 10.0*Ga1*beta1*ha0**3*ha1*ha3 - 28.0*Ga1*beta1*ha0**3*ha2**2 + 33.0*Ga1*beta1*ha0**2*ha1**2*ha2 + 11.0*Ga1*beta1*ha0*ha1**4 + 2.0*Ga1*ha0**2*ha2 - 2.0*Ga1*ha0*ha1**2 + 66.0*Ga2*beta1**2*ha0**5*ha2*ha3 - 72.0*Ga2*beta1**2*ha0**4*ha1**2*ha3 + 84.0*Ga2*beta1**2*ha0**4*ha1*ha2**2 + 36.0*Ga2*beta1**2*ha0**3*ha1**3*ha2 + 12.0*Ga2*beta1**2*ha0**2*ha1**5 - 22.0*Ga2*beta1*ha0**4*ha3 - 44.0*Ga2*beta1*ha0**3*ha1*ha2 - 11.0*Ga2*beta1*ha0**2*ha1**3 + 2.0*Ga2*ha0**2*ha1 + 45.0*Ga3*beta1**2*ha0**5*ha1*ha3 - 48.0*Ga3*beta1**2*ha0**5*ha2**2 - 3.0*Ga3*beta1**2*ha0**4*ha1**2*ha2 - 12.0*Ga3*beta1**2*ha0**3*ha1**4 + 22.0*Ga3*beta1*ha0**4*ha2 + 11.0*Ga3*beta1*ha0**3*ha1**2 - 2.0*Ga3*ha0**3)/(ha0**4*(594.0*beta1**3*ha0**4*ha3**2 - 2124.0*beta1**3*ha0**3*ha1*ha2*ha3 + 1152.0*beta1**3*ha0**3*ha2**3 - 396.0*beta1**3*ha0**2*ha1**3*ha3 + 504.0*beta1**3*ha0**2*ha1**2*ha2**2 - 288.0*beta1**3*ha0*ha1**4*ha2 + 180.0*beta1**3*ha1**6 + 624.0*beta1**2*ha0**2*ha1*ha3 - 720.0*beta1**2*ha0**2*ha2**2 - 204.0*beta1**2*ha0*ha1**2*ha2 - 189.0*beta1**2*ha1**4 + 104.0*beta1*ha0*ha2 + 52.0*beta1*ha1**2 - 4.0))
#         ua2 = (54.0*Ga0*beta1**2*ha0**4*ha3**2 - 522.0*Ga0*beta1**2*ha0**3*ha1*ha2*ha3 + 144.0*Ga0*beta1**2*ha0**3*ha2**3 + 180.0*Ga0*beta1**2*ha0**2*ha1**3*ha3 + 36.0*Ga0*beta1**2*ha0**2*ha1**2*ha2**2 - 216.0*Ga0*beta1**2*ha0*ha1**4*ha2 - 108.0*Ga0*beta1**2*ha1**6 + 60.0*Ga0*beta1*ha0**2*ha1*ha3 - 72.0*Ga0*beta1*ha0**2*ha2**2 - 30.0*Ga0*beta1*ha0*ha1**2*ha2 + 78.0*Ga0*beta1*ha1**4 + 4.0*Ga0*ha0*ha2 - 4.0*Ga0*ha1**2 + 324.0*Ga1*beta1**2*ha0**4*ha2*ha3 - 243.0*Ga1*beta1**2*ha0**3*ha1**2*ha3 + 216.0*Ga1*beta1**2*ha0**3*ha1*ha2**2 + 324.0*Ga1*beta1**2*ha0**2*ha1**3*ha2 + 108.0*Ga1*beta1**2*ha0*ha1**5 - 18.0*Ga1*beta1*ha0**3*ha3 - 48.0*Ga1*beta1*ha0**2*ha1*ha2 - 78.0*Ga1*beta1*ha0*ha1**3 + 4.0*Ga1*ha0*ha1 + 126.0*Ga2*beta1**2*ha0**4*ha1*ha3 - 144.0*Ga2*beta1**2*ha0**4*ha2**2 - 252.0*Ga2*beta1**2*ha0**3*ha1**2*ha2 - 108.0*Ga2*beta1**2*ha0**2*ha1**4 + 72.0*Ga2*beta1*ha0**3*ha2 + 78.0*Ga2*beta1*ha0**2*ha1**2 - 4.0*Ga2*ha0**2 - 54.0*Ga3*beta1**2*ha0**5*ha3 + 72.0*Ga3*beta1**2*ha0**4*ha1*ha2 + 63.0*Ga3*beta1**2*ha0**3*ha1**3 - 42.0*Ga3*beta1*ha0**3*ha1)/(ha0**3*(594.0*beta1**3*ha0**4*ha3**2 - 2124.0*beta1**3*ha0**3*ha1*ha2*ha3 + 1152.0*beta1**3*ha0**3*ha2**3 - 396.0*beta1**3*ha0**2*ha1**3*ha3 + 504.0*beta1**3*ha0**2*ha1**2*ha2**2 - 288.0*beta1**3*ha0*ha1**4*ha2 + 180.0*beta1**3*ha1**6 + 624.0*beta1**2*ha0**2*ha1*ha3 - 720.0*beta1**2*ha0**2*ha2**2 - 204.0*beta1**2*ha0*ha1**2*ha2 - 189.0*beta1**2*ha1**4 + 104.0*beta1*ha0*ha2 + 52.0*beta1*ha1**2 - 4.0))
#         ua1 = 2.0*(18.0*Ga0*beta1**2*ha0**3*ha2*ha3 - 216.0*Ga0*beta1**2*ha0**2*ha1**2*ha3 + 252.0*Ga0*beta1**2*ha0**2*ha1*ha2**2 - 192.0*Ga0*beta1**2*ha0*ha1**3*ha2 + 216.0*Ga0*beta1**2*ha1**5 + 6.0*Ga0*beta1*ha0**2*ha3 - 48.0*Ga0*beta1*ha0*ha1*ha2 - 27.0*Ga0*beta1*ha1**3 + 2.0*Ga0*ha1 + 297.0*Ga1*beta1**2*ha0**3*ha1*ha3 - 288.0*Ga1*beta1**2*ha0**3*ha2**2 + 36.0*Ga1*beta1**2*ha0**2*ha1**2*ha2 - 216.0*Ga1*beta1**2*ha0*ha1**4 + 52.0*Ga1*beta1*ha0**2*ha2 + 27.0*Ga1*beta1*ha0*ha1**2 - 2.0*Ga1*ha0 - 66.0*Ga2*beta1**2*ha0**4*ha3 + 36.0*Ga2*beta1**2*ha0**3*ha1*ha2 + 156.0*Ga2*beta1**2*ha0**2*ha1**3 - 4.0*Ga2*beta1*ha0**2*ha1 + 48.0*Ga3*beta1**2*ha0**4*ha2 - 81.0*Ga3*beta1**2*ha0**3*ha1**2 - 6.0*Ga3*beta1*ha0**3)/(ha0**2*(594.0*beta1**3*ha0**4*ha3**2 - 2124.0*beta1**3*ha0**3*ha1*ha2*ha3 + 1152.0*beta1**3*ha0**3*ha2**3 - 396.0*beta1**3*ha0**2*ha1**3*ha3 + 504.0*beta1**3*ha0**2*ha1**2*ha2**2 - 288.0*beta1**3*ha0*ha1**4*ha2 + 180.0*beta1**3*ha1**6 + 624.0*beta1**2*ha0**2*ha1*ha3 - 720.0*beta1**2*ha0**2*ha2**2 - 204.0*beta1**2*ha0*ha1**2*ha2 - 189.0*beta1**2*ha1**4 + 104.0*beta1*ha0*ha2 + 52.0*beta1*ha1**2 - 4.0))
#         ua0 = (648.0*Ga0*beta1**3*ha0**4*ha3**2 - 2592.0*Ga0*beta1**3*ha0**3*ha1*ha2*ha3 + 1296.0*Ga0*beta1**3*ha0**3*ha2**3 - 864.0*Ga0*beta1**3*ha0**2*ha1**3*ha3 + 1296.0*Ga0*beta1**3*ha0**2*ha1**2*ha2**2 - 1080.0*Ga0*beta1**3*ha0*ha1**4*ha2 + 720.0*Ga0*beta1**3*ha1**6 + 702.0*Ga0*beta1**2*ha0**2*ha1*ha3 - 792.0*Ga0*beta1**2*ha0**2*ha2**2 - 378.0*Ga0*beta1**2*ha0*ha1**2*ha2 - 192.0*Ga0*beta1**2*ha1**4 + 108.0*Ga0*beta1*ha0*ha2 + 54.0*Ga0*beta1*ha1**2 - 4.0*Ga0 + 324.0*Ga1*beta1**3*ha0**4*ha2*ha3 + 648.0*Ga1*beta1**3*ha0**3*ha1**2*ha3 - 648.0*Ga1*beta1**3*ha0**3*ha1*ha2**2 + 432.0*Ga1*beta1**3*ha0**2*ha1**3*ha2 - 540.0*Ga1*beta1**3*ha0*ha1**5 - 18.0*Ga1*beta1**2*ha0**3*ha3 + 108.0*Ga1*beta1**2*ha0**2*ha1*ha2 + 3.0*Ga1*beta1**2*ha0*ha1**3 - 2.0*Ga1*beta1*ha0*ha1 - 72.0*Ga2*beta1**3*ha0**4*ha1*ha3 - 144.0*Ga2*beta1**3*ha0**4*ha2**2 - 144.0*Ga2*beta1**3*ha0**3*ha1**2*ha2 + 360.0*Ga2*beta1**3*ha0**2*ha1**4 + 72.0*Ga2*beta1**2*ha0**3*ha2 + 66.0*Ga2*beta1**2*ha0**2*ha1**2 - 4.0*Ga2*beta1*ha0**2 - 54.0*Ga3*beta1**3*ha0**5*ha3 + 216.0*Ga3*beta1**3*ha0**4*ha1*ha2 - 180.0*Ga3*beta1**3*ha0**3*ha1**3 - 60.0*Ga3*beta1**2*ha0**3*ha1)/(ha0*(594.0*beta1**3*ha0**4*ha3**2 - 2124.0*beta1**3*ha0**3*ha1*ha2*ha3 + 1152.0*beta1**3*ha0**3*ha2**3 - 396.0*beta1**3*ha0**2*ha1**3*ha3 + 504.0*beta1**3*ha0**2*ha1**2*ha2**2 - 288.0*beta1**3*ha0*ha1**4*ha2 + 180.0*beta1**3*ha1**6 + 624.0*beta1**2*ha0**2*ha1*ha3 - 720.0*beta1**2*ha0**2*ha2**2 - 204.0*beta1**2*ha0*ha1**2*ha2 - 189.0*beta1**2*ha1**4 + 104.0*beta1*ha0*ha2 + 52.0*beta1*ha1**2 - 4.0))
        
#         # if (abs(ha3) > 10.0**(-12)):
#         #     ua3 =(108*Ga0*beta1**2*ha0**4*ha1*ha3**2 - 36*Ga0*beta1**2*ha0**4*ha2**2*ha3 + 576*Ga0*beta1**2*ha0**3*ha1**2*ha2*ha3 - 360*Ga0*beta1**2*ha0**3*ha1*ha2**3 - 216*Ga0*beta1**2*ha0**2*ha1**4*ha3 + 240*Ga0*beta1**2*ha0**2*ha1**3*ha2**2 + 24*Ga0*beta1**2*ha0*ha1**5*ha2 + 24*Ga0*beta1**2*ha1**7 - 21*Ga0*beta1*ha0**2*ha1**2*ha3 + 72*Ga0*beta1*ha0**2*ha1*ha2**2 - 22*Ga0*beta1*ha0*ha1**3*ha2 - 11*Ga0*beta1*ha1**5 + Ga0*ha0**2*ha3 - 2*Ga0*ha0*ha1*ha2 + Ga0*ha1**3 - 198*Ga1*beta1**2*ha0**5*ha3**2 - 426*Ga1*beta1**2*ha0**4*ha1*ha2*ha3 + 192*Ga1*beta1**2*ha0**4*ha2**3 + 240*Ga1*beta1**2*ha0**3*ha1**3*ha3 - 312*Ga1*beta1**2*ha0**3*ha1**2*ha2**2 - 48*Ga1*beta1**2*ha0**2*ha1**4*ha2 - 24*Ga1*beta1**2*ha0*ha1**6 + 10*Ga1*beta1*ha0**3*ha1*ha3 - 28*Ga1*beta1*ha0**3*ha2**2 + 33*Ga1*beta1*ha0**2*ha1**2*ha2 + 11*Ga1*beta1*ha0*ha1**4 + Ga1*ha0**2*ha2 - Ga1*ha0*ha1**2 + 132*Ga2*beta1**2*ha0**5*ha2*ha3 - 144*Ga2*beta1**2*ha0**4*ha1**2*ha3 + 168*Ga2*beta1**2*ha0**4*ha1*ha2**2 + 72*Ga2*beta1**2*ha0**3*ha1**3*ha2 + 24*Ga2*beta1**2*ha0**2*ha1**5 - 22*Ga2*beta1*ha0**4*ha3 - 44*Ga2*beta1*ha0**3*ha1*ha2 - 11*Ga2*beta1*ha0**2*ha1**3 + Ga2*ha0**2*ha1 + 90*Ga3*beta1**2*ha0**5*ha1*ha3 - 96*Ga3*beta1**2*ha0**5*ha2**2 - 6*Ga3*beta1**2*ha0**4*ha1**2*ha2 - 24*Ga3*beta1**2*ha0**3*ha1**4 + 22*Ga3*beta1*ha0**4*ha2 + 11*Ga3*beta1*ha0**3*ha1**2 - Ga3*ha0**3)/(ha0**4*(1188*beta1**3*ha0**4*ha3**2 - 4248*beta1**3*ha0**3*ha1*ha2*ha3 + 2304*beta1**3*ha0**3*ha2**3 - 792*beta1**3*ha0**2*ha1**3*ha3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 + 624*beta1**2*ha0**2*ha1*ha3 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
#         #     ua2 = (54*Ga0*beta1**2*ha0**4*ha3**2 - 522*Ga0*beta1**2*ha0**3*ha1*ha2*ha3 + 144*Ga0*beta1**2*ha0**3*ha2**3 + 180*Ga0*beta1**2*ha0**2*ha1**3*ha3 + 36*Ga0*beta1**2*ha0**2*ha1**2*ha2**2 - 216*Ga0*beta1**2*ha0*ha1**4*ha2 - 108*Ga0*beta1**2*ha1**6 + 30*Ga0*beta1*ha0**2*ha1*ha3 - 36*Ga0*beta1*ha0**2*ha2**2 - 15*Ga0*beta1*ha0*ha1**2*ha2 + 39*Ga0*beta1*ha1**4 + Ga0*ha0*ha2 - Ga0*ha1**2 + 324*Ga1*beta1**2*ha0**4*ha2*ha3 - 243*Ga1*beta1**2*ha0**3*ha1**2*ha3 + 216*Ga1*beta1**2*ha0**3*ha1*ha2**2 + 324*Ga1*beta1**2*ha0**2*ha1**3*ha2 + 108*Ga1*beta1**2*ha0*ha1**5 - 9*Ga1*beta1*ha0**3*ha3 - 24*Ga1*beta1*ha0**2*ha1*ha2 - 39*Ga1*beta1*ha0*ha1**3 + Ga1*ha0*ha1 + 126*Ga2*beta1**2*ha0**4*ha1*ha3 - 144*Ga2*beta1**2*ha0**4*ha2**2 - 252*Ga2*beta1**2*ha0**3*ha1**2*ha2 - 108*Ga2*beta1**2*ha0**2*ha1**4 + 36*Ga2*beta1*ha0**3*ha2 + 39*Ga2*beta1*ha0**2*ha1**2 - Ga2*ha0**2 - 54*Ga3*beta1**2*ha0**5*ha3 + 72*Ga3*beta1**2*ha0**4*ha1*ha2 + 63*Ga3*beta1**2*ha0**3*ha1**3 - 21*Ga3*beta1*ha0**3*ha1)/(ha0**3*(1188*beta1**3*ha0**4*ha3**2 - 4248*beta1**3*ha0**3*ha1*ha2*ha3 + 2304*beta1**3*ha0**3*ha2**3 - 792*beta1**3*ha0**2*ha1**3*ha3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 + 624*beta1**2*ha0**2*ha1*ha3 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
#         #     ua1 = (36*Ga0*beta1**2*ha0**3*ha2*ha3 - 432*Ga0*beta1**2*ha0**2*ha1**2*ha3 + 504*Ga0*beta1**2*ha0**2*ha1*ha2**2 - 384*Ga0*beta1**2*ha0*ha1**3*ha2 + 432*Ga0*beta1**2*ha1**5 + 6*Ga0*beta1*ha0**2*ha3 - 48*Ga0*beta1*ha0*ha1*ha2 - 27*Ga0*beta1*ha1**3 + Ga0*ha1 + 594*Ga1*beta1**2*ha0**3*ha1*ha3 - 576*Ga1*beta1**2*ha0**3*ha2**2 + 72*Ga1*beta1**2*ha0**2*ha1**2*ha2 - 432*Ga1*beta1**2*ha0*ha1**4 + 52*Ga1*beta1*ha0**2*ha2 + 27*Ga1*beta1*ha0*ha1**2 - Ga1*ha0 - 132*Ga2*beta1**2*ha0**4*ha3 + 72*Ga2*beta1**2*ha0**3*ha1*ha2 + 312*Ga2*beta1**2*ha0**2*ha1**3 - 4*Ga2*beta1*ha0**2*ha1 + 96*Ga3*beta1**2*ha0**4*ha2 - 162*Ga3*beta1**2*ha0**3*ha1**2 - 6*Ga3*beta1*ha0**3)/(ha0**2*(1188*beta1**3*ha0**4*ha3**2 - 4248*beta1**3*ha0**3*ha1*ha2*ha3 + 2304*beta1**3*ha0**3*ha2**3 - 792*beta1**3*ha0**2*ha1**3*ha3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 + 624*beta1**2*ha0**2*ha1*ha3 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
#         #     ua0 = (1296*Ga0*beta1**3*ha0**4*ha3**2 - 5184*Ga0*beta1**3*ha0**3*ha1*ha2*ha3 + 2592*Ga0*beta1**3*ha0**3*ha2**3 - 1728*Ga0*beta1**3*ha0**2*ha1**3*ha3 + 2592*Ga0*beta1**3*ha0**2*ha1**2*ha2**2 - 2160*Ga0*beta1**3*ha0*ha1**4*ha2 + 1440*Ga0*beta1**3*ha1**6 + 702*Ga0*beta1**2*ha0**2*ha1*ha3 - 792*Ga0*beta1**2*ha0**2*ha2**2 - 378*Ga0*beta1**2*ha0*ha1**2*ha2 - 192*Ga0*beta1**2*ha1**4 + 54*Ga0*beta1*ha0*ha2 + 27*Ga0*beta1*ha1**2 - Ga0 + 648*Ga1*beta1**3*ha0**4*ha2*ha3 + 1296*Ga1*beta1**3*ha0**3*ha1**2*ha3 - 1296*Ga1*beta1**3*ha0**3*ha1*ha2**2 + 864*Ga1*beta1**3*ha0**2*ha1**3*ha2 - 1080*Ga1*beta1**3*ha0*ha1**5 - 18*Ga1*beta1**2*ha0**3*ha3 + 108*Ga1*beta1**2*ha0**2*ha1*ha2 + 3*Ga1*beta1**2*ha0*ha1**3 - Ga1*beta1*ha0*ha1 - 144*Ga2*beta1**3*ha0**4*ha1*ha3 - 288*Ga2*beta1**3*ha0**4*ha2**2 - 288*Ga2*beta1**3*ha0**3*ha1**2*ha2 + 720*Ga2*beta1**3*ha0**2*ha1**4 + 72*Ga2*beta1**2*ha0**3*ha2 + 66*Ga2*beta1**2*ha0**2*ha1**2 - 2*Ga2*beta1*ha0**2 - 108*Ga3*beta1**3*ha0**5*ha3 + 432*Ga3*beta1**3*ha0**4*ha1*ha2 - 360*Ga3*beta1**3*ha0**3*ha1**3 - 60*Ga3*beta1**2*ha0**3*ha1)/(ha0*(1188*beta1**3*ha0**4*ha3**2 - 4248*beta1**3*ha0**3*ha1*ha2*ha3 + 2304*beta1**3*ha0**3*ha2**3 - 792*beta1**3*ha0**2*ha1**3*ha3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 + 624*beta1**2*ha0**2*ha1*ha3 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
#         # elif (abs(ha2) > 10.0**(-12)):
#         #     ua3 = (-360*Ga0*beta1**2*ha0**3*ha1*ha2**3 + 240*Ga0*beta1**2*ha0**2*ha1**3*ha2**2 + 24*Ga0*beta1**2*ha0*ha1**5*ha2 + 24*Ga0*beta1**2*ha1**7 + 72*Ga0*beta1*ha0**2*ha1*ha2**2 - 22*Ga0*beta1*ha0*ha1**3*ha2 - 11*Ga0*beta1*ha1**5 - 2*Ga0*ha0*ha1*ha2 + Ga0*ha1**3 + 192*Ga1*beta1**2*ha0**4*ha2**3 - 312*Ga1*beta1**2*ha0**3*ha1**2*ha2**2 - 48*Ga1*beta1**2*ha0**2*ha1**4*ha2 - 24*Ga1*beta1**2*ha0*ha1**6 - 28*Ga1*beta1*ha0**3*ha2**2 + 33*Ga1*beta1*ha0**2*ha1**2*ha2 + 11*Ga1*beta1*ha0*ha1**4 + Ga1*ha0**2*ha2 - Ga1*ha0*ha1**2 + 168*Ga2*beta1**2*ha0**4*ha1*ha2**2 + 72*Ga2*beta1**2*ha0**3*ha1**3*ha2 + 24*Ga2*beta1**2*ha0**2*ha1**5 - 44*Ga2*beta1*ha0**3*ha1*ha2 - 11*Ga2*beta1*ha0**2*ha1**3 + Ga2*ha0**2*ha1 - 96*Ga3*beta1**2*ha0**5*ha2**2 - 6*Ga3*beta1**2*ha0**4*ha1**2*ha2 - 24*Ga3*beta1**2*ha0**3*ha1**4 + 22*Ga3*beta1*ha0**4*ha2 + 11*Ga3*beta1*ha0**3*ha1**2 - Ga3*ha0**3)/(ha0**4*(2304*beta1**3*ha0**3*ha2**3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
#         #     ua2 = (144*Ga0*beta1**2*ha0**3*ha2**3 + 36*Ga0*beta1**2*ha0**2*ha1**2*ha2**2 - 216*Ga0*beta1**2*ha0*ha1**4*ha2 - 108*Ga0*beta1**2*ha1**6 - 36*Ga0*beta1*ha0**2*ha2**2 - 15*Ga0*beta1*ha0*ha1**2*ha2 + 39*Ga0*beta1*ha1**4 + Ga0*ha0*ha2 - Ga0*ha1**2 + 216*Ga1*beta1**2*ha0**3*ha1*ha2**2 + 324*Ga1*beta1**2*ha0**2*ha1**3*ha2 + 108*Ga1*beta1**2*ha0*ha1**5 - 24*Ga1*beta1*ha0**2*ha1*ha2 - 39*Ga1*beta1*ha0*ha1**3 + Ga1*ha0*ha1 - 144*Ga2*beta1**2*ha0**4*ha2**2 - 252*Ga2*beta1**2*ha0**3*ha1**2*ha2 - 108*Ga2*beta1**2*ha0**2*ha1**4 + 36*Ga2*beta1*ha0**3*ha2 + 39*Ga2*beta1*ha0**2*ha1**2 - Ga2*ha0**2 + 72*Ga3*beta1**2*ha0**4*ha1*ha2 + 63*Ga3*beta1**2*ha0**3*ha1**3 - 21*Ga3*beta1*ha0**3*ha1)/(ha0**3*(2304*beta1**3*ha0**3*ha2**3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
#         #     ua1 = (504*Ga0*beta1**2*ha0**2*ha1*ha2**2 - 384*Ga0*beta1**2*ha0*ha1**3*ha2 + 432*Ga0*beta1**2*ha1**5 - 48*Ga0*beta1*ha0*ha1*ha2 - 27*Ga0*beta1*ha1**3 + Ga0*ha1 - 576*Ga1*beta1**2*ha0**3*ha2**2 + 72*Ga1*beta1**2*ha0**2*ha1**2*ha2 - 432*Ga1*beta1**2*ha0*ha1**4 + 52*Ga1*beta1*ha0**2*ha2 + 27*Ga1*beta1*ha0*ha1**2 - Ga1*ha0 + 72*Ga2*beta1**2*ha0**3*ha1*ha2 + 312*Ga2*beta1**2*ha0**2*ha1**3 - 4*Ga2*beta1*ha0**2*ha1 + 96*Ga3*beta1**2*ha0**4*ha2 - 162*Ga3*beta1**2*ha0**3*ha1**2 - 6*Ga3*beta1*ha0**3)/(ha0**2*(2304*beta1**3*ha0**3*ha2**3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
#         #     ua0 =(2592*Ga0*beta1**3*ha0**3*ha2**3 + 2592*Ga0*beta1**3*ha0**2*ha1**2*ha2**2 - 2160*Ga0*beta1**3*ha0*ha1**4*ha2 + 1440*Ga0*beta1**3*ha1**6 - 792*Ga0*beta1**2*ha0**2*ha2**2 - 378*Ga0*beta1**2*ha0*ha1**2*ha2 - 192*Ga0*beta1**2*ha1**4 + 54*Ga0*beta1*ha0*ha2 + 27*Ga0*beta1*ha1**2 - Ga0 - 1296*Ga1*beta1**3*ha0**3*ha1*ha2**2 + 864*Ga1*beta1**3*ha0**2*ha1**3*ha2 - 1080*Ga1*beta1**3*ha0*ha1**5 + 108*Ga1*beta1**2*ha0**2*ha1*ha2 + 3*Ga1*beta1**2*ha0*ha1**3 - Ga1*beta1*ha0*ha1 - 288*Ga2*beta1**3*ha0**4*ha2**2 - 288*Ga2*beta1**3*ha0**3*ha1**2*ha2 + 720*Ga2*beta1**3*ha0**2*ha1**4 + 72*Ga2*beta1**2*ha0**3*ha2 + 66*Ga2*beta1**2*ha0**2*ha1**2 - 2*Ga2*beta1*ha0**2 + 432*Ga3*beta1**3*ha0**4*ha1*ha2 - 360*Ga3*beta1**3*ha0**3*ha1**3 - 60*Ga3*beta1**2*ha0**3*ha1)/(ha0*(2304*beta1**3*ha0**3*ha2**3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
#         # elif (abs(ha1) > 10.0**(-12)):
#         #     ua3 = (Ga0*ha1**3 - Ga1*ha0*ha1**2 + Ga2*ha0**2*ha1 - Ga3*ha0**3)/(ha0**4*(15*beta1*ha1**2 - 1))
#         #     ua2 = (-36*Ga0*beta1*ha1**4 + Ga0*ha1**2 + 36*Ga1*beta1*ha0*ha1**3 - Ga1*ha0*ha1 - 36*Ga2*beta1*ha0**2*ha1**2 + Ga2*ha0**2 + 21*Ga3*beta1*ha0**3*ha1)/(ha0**3*(120*beta1**2*ha1**4 - 23*beta1*ha1**2 + 1))
#         #     ua1 = (432*Ga0*beta1**2*ha1**5 - 27*Ga0*beta1*ha1**3 + Ga0*ha1 - 432*Ga1*beta1**2*ha0*ha1**4 + 27*Ga1*beta1*ha0*ha1**2 - Ga1*ha0 + 312*Ga2*beta1**2*ha0**2*ha1**3 - 4*Ga2*beta1*ha0**2*ha1 - 162*Ga3*beta1**2*ha0**3*ha1**2 - 6*Ga3*beta1*ha0**3)/(ha0**2*(360*beta1**3*ha1**6 - 189*beta1**2*ha1**4 + 26*beta1*ha1**2 - 1))
#         #     ua0 = (1440*Ga0*beta1**3*ha1**6 - 192*Ga0*beta1**2*ha1**4 + 27*Ga0*beta1*ha1**2 - Ga0 - 1080*Ga1*beta1**3*ha0*ha1**5 + 3*Ga1*beta1**2*ha0*ha1**3 - Ga1*beta1*ha0*ha1 + 720*Ga2*beta1**3*ha0**2*ha1**4 + 66*Ga2*beta1**2*ha0**2*ha1**2 - 2*Ga2*beta1*ha0**2 - 360*Ga3*beta1**3*ha0**3*ha1**3 - 60*Ga3*beta1**2*ha0**3*ha1)/(ha0*(360*beta1**3*ha1**6 - 189*beta1**2*ha1**4 + 26*beta1*ha1**2 - 1))
#         # else:
#         #     ua3 =Ga3/ha0
#         #     ua2 =Ga2/ha0
#         #     ua1 =Ga1/ha0 + 6*Ga3*beta1*ha0
#         #     ua0 =Ga0/ha0 + 2*Ga2*beta1*ha0
            
#     else:
#         hjmh = PointsFromPolyFromPoints(ha0,ha1,ha2,ha3,-dx/2)
#         hjms = PointsFromPolyFromPoints(ha0,ha1,ha2,ha3,-dx/6)
#         hjps = PointsFromPolyFromPoints(ha0,ha1,ha2,ha3,dx/6)
#         hjph = PointsFromPolyFromPoints(ha0,ha1,ha2,ha3,dx/2)
        
#         Gjmh = PointsFromPolyFromPoints(Ga0,Ga1,Ga2,Ga3,-dx/2)
#         Gjms = PointsFromPolyFromPoints(Ga0,Ga1,Ga2,Ga3,-dx/6)
#         Gjps = PointsFromPolyFromPoints(Ga0,Ga1,Ga2,Ga3,dx/6)
#         Gjph = PointsFromPolyFromPoints(Ga0,Ga1,Ga2,Ga3,dx/2)
               
#         ujmh  =  Gjmh/ hjmh  
#         ujms =  Gjms/ hjms  
#         ujps =  Gjps/ hjps  
#         ujph   =  Gjph/ hjph      

#         ua0,ua1,ua2,ua3=PolyFromPoints(ujmh,ujms,ujps,ujph,dx)
                                   
#     return ua0,ua1,ua2,ua3


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

# def FEM(beta1,h,G,n,x,dx):
#     up = []
#     for i in range(n):
#         ha0,ha1,ha2,ha3 = PolyFromPoints(h[4*i],h[4*i+1],h[4*i+2],h[4*i+3],dx)
#         Ga0,Ga1,Ga2,Ga3 = PolyFromPoints(G[4*i],G[4*i+1],G[4*i+2],G[4*i+3],dx)
#         print('Initial h', ha0,ha1,ha2,ha3)
#         print('Initial G', Ga0,Ga1,Ga2,Ga3)
#         ua0,ua1,ua2,ua3 = SolveForuEdges(beta1,ha0,ha1,ha2,ha3,Ga0,Ga1,Ga2,Ga3 ,dx)
#         print('U Coeff', ua0,ua1,ua2,ua3)
#         ujmh = PointsFromPolyFromPoints(ua0,ua1,ua2,ua3,-dx/2)
#         ujms = PointsFromPolyFromPoints(ua0,ua1,ua2,ua3,-dx/6)
#         ujps = PointsFromPolyFromPoints(ua0,ua1,ua2,ua3,dx/6)
#         ujph = PointsFromPolyFromPoints(ua0,ua1,ua2,ua3,dx/2)
#         print('U Value',ujmh,ujms,ujps,ujph)
#         ucell = [ujmh,ujms,ujps,ujph]
#         up = up + ucell
    
#     return array(up)
    
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

Gn = FEMforG(b1,h,u,n+1,x,dx)

# un = FEM(b1,h,G,len(xn),x,dx)

# i = n//2
# ha0,ha1,ha2,ha3 = PolyFromPoints(h[4*i],h[4*i+1],h[4*i+2],h[4*i+3],dx)
# ua0,ua1,ua2,ua3 = PolyFromPoints(u[4*i],u[4*i+1],u[4*i+2],u[4*i+3],dx)
# Gjmh, Gjms, Gjps,  Gjph = CalculateG(b1,ha0,ha1,ha2,ha3,ua0,ua1,ua2,ua3 ,dx)
# print(Gjmh, Gjms, Gjps,  Gjph)
# print(G[4*i:4*(i+1)])
# G = FEMforG(xn,h,u,dx)