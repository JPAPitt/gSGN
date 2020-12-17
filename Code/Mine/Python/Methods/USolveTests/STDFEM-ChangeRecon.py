
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from matplotlib.pyplot import plot, loglog
from numpy.linalg import solve,lstsq,cond,norm

def minmodL(q):
    qmin = min(q)
    qmax = max(q)
    if qmin>0:
        return qmin
    if qmax<0 :
        return qmax
    else:
        return 0

#linear
def P1jtojp1(qaj,qajp1,dx):
    return (qajp1 - qaj)/ dx ,qaj

def P1jm1toj(qaj,qajm1,dx):
    return (qaj - qajm1)/ dx ,qaj

def P2jm1tojp1(qajm1,qaj,qajp1,dx):
    a = (qajp1 - 2*qaj + qajm1)/ (2*dx**2)
    b = (qajp1 - qajm1)/ (2*dx)
    c =  13*qaj/12 - qajp1/24 - qajm1/24
    return a,b,c

#Quadtraitc
def P2jtojp2(qaj,qajp1,qajp2,dx):
    a = (qajp2 - 2*qajp1 + qaj)/ (2*dx**2)
    b = (-qajp2 + 4*qajp1 - 3*qaj)/ (2*dx)
    c =  -qajp2/24 + qajp1/12 + 23*qaj/24
    return a,b,c

def P2jm2toj(qajm2,qajm1,qaj,dx):
    a = (qaj - 2*qajm1 + qajm2)/ (2*dx**2)
    b = (qajm2 - 4*qajm1 + 3*qaj)/ (2*dx)
    c =  -qajm2/24 + qajm1/12 + 23*qaj/24
    return a,b,c

def P4jm2tojp2(qajm2,qajm1,qaj,qajp1,qajp2,dx):
    a = -qajp1/(6*dx**4) + qajp2/(24*dx**4) - qajm1/(6*dx**4) + qajm2/(24*dx**4) + qaj/(4*dx**4)
    b = -qajp1/(6*dx**3) + qajp2/(12*dx**3) + qajm1/(6*dx**3) - qajm2/(12*dx**3)
    c = 3*qajp1/(4*dx**2) - qajp2/(16*dx**2) + 3*qajm1/(4*dx**2) - qajm2/(16*dx**2) - 11*qaj/(8*dx**2)
    d = 17*qajp1/(24*dx) - 5*qajp2/(48*dx) - 17*qajm1/(24*dx) + 5*qajm2/(48*dx)
    e = -29*qajp1/480 + 3*qajp2/640 - 29*qajm1/480 + 3*qajm2/640 + 1067*qaj/960
    return a,b,c,d,e

#Cubic
def P3jtojp3(qaj,qajp1,qajp2,qajp3,dx):
    a=(3*qajp1 - 3*qajp2 + qajp3 - qaj)/(6*dx**3)
    b=(-5*qajp1 + 4*qajp2 - qajp3 + 2*qaj)/(2*dx**2)
    c=(69*qajp1 - 33*qajp2 + 7*qajp3 - 43*qaj)/(24*dx)
    d=5*qajp1/24 - qajp2/6 + qajp3/24 + 11*qaj/12
    return a,b,c,d

def P3jm3toj(qaj,qajm1,qajm2,qajm3,dx):
    a=(-3*qajm1 + 3*qajm2 - qajm3 + qaj)/(6*dx**3)
    b=(-5*qajm1 + 4*qajm2 - qajm3 + 2*qaj)/(2*dx**2)
    c=(-69*qajm1 + 33*qajm2 - 7*qajm3 + 43*qaj)/(24*dx)
    d=5*qajm1/24 - qajm2/6 + qajm3/24 + 11*qaj/12
    return a,b,c,d

def P6jm3tojp3(qajm3,qajm2,qajm1,qaj,qajp1,qajp2,qajp3,dx):
    a = qajp1/(48*dx**6) - qajp2/(120*dx**6) + qajp3/(720*dx**6) + qajm1/(48*dx**6) - qajm2/(120*dx**6) + qajm3/(720*dx**6) - qaj/(36*dx**6)
    b = qajp1/(48*dx**5) - qajp2/(60*dx**5) + qajp3/(240*dx**5) - qajm1/(48*dx**5) + qajm2/(60*dx**5) - qajm3/(240*dx**5)
    c = -19*qajp1/(64*dx**4) + 3*qajp2/(32*dx**4) - 5*qajp3/(576*dx**4) - 19*qajm1/(64*dx**4) + 3*qajm2/(32*dx**4) - 5*qajm3/(576*dx**4) + 61*qaj/(144*dx**4)
    d = -83*qajp1/(288*dx**3) + 13*qajp2/(72*dx**3) - 7*qajp3/(288*dx**3) + 83*qajm1/(288*dx**3) - 13*qajm2/(72*dx**3) + 7*qajm3/(288*dx**3)
    e = 229*qajp1/(256*dx**2) - 77*qajp2/(640*dx**2) + 37*qajp3/(3840*dx**2) + 229*qajm1/(256*dx**2) - 77*qajm2/(640*dx**2) + 37*qajm3/(3840*dx**2) - 301*qaj/(192*dx**2)
    f = 1891*qajp1/(2304*dx) - 559*qajp2/(2880*dx) + 259*qajp3/(11520*dx) - 1891*qajm1/(2304*dx) + 559*qajm2/(2880*dx) - 259*qajm3/(11520*dx)
    g = -7621*qajp1/107520 + 159*qajp2/17920 - 5*qajp3/7168 - 7621*qajm1/107520 + 159*qajm2/17920 - 5*qajm3/7168 + 30251*qaj/26880

    return a,b,c,d,e,f,g

def Recon(qA,x,dx,nG):
    n = len(qA)
    qcubic = zeros(4*n)
    
    # j= 0
    j = 0
    qcubic[4*j: 4*(j+nG)] = qA[j]*ones(4*nG)
    
    for j in range(nG,n-nG):

        P3r0a33,P3r0a32,P3r0a31,P3r0a30 = P3jtojp3(qA[j],qA[j+1],qA[j+2],qA[j+3],dx)
        # P3r1a33,P3r1a32,P3r1a31,P3r1a30 = P3jm1tojp2(qA[j],qA[j-1],qA[j+1],qA[j+2],dx)
        # P3r2a33,P3r2a32,P3r2a31,P3r2a30 = P3jm2tojp1(qA[j],qA[j-1],qA[j-2],qA[j+1],dx)
        P3r3a33,P3r3a32,P3r3a31,P3r3a30 =  P3jm3toj(qA[j],qA[j-1],qA[j-2],qA[j-3],dx)
        P6ca66,P6ca65,P6ca64,P6ca63,P6ca62,P6ca61,P6ca60 = P6jm3tojp3(qA[j-3],qA[j-2],qA[j-1],qA[j],qA[j+1],qA[j+2],qA[j+3],dx)
        
        P2r0a22,P2r0a21,P2r0a20=   P2jtojp2(qA[j],qA[j+1],qA[j+2],dx)
        # P2r1a22,P2r1a21,P2r1a20  =  P2jm1tojp1(qA[j-1],qA[j],qA[j+1],dx)
        P2r2a22,P2r2a21,P2r2a20  =  P2jm2toj(qA[j-2],qA[j-1],qA[j],dx)
        P4ca44,P4ca43,P4ca42,P4ca41,P4ca40 = P4jm2tojp2(qA[j-2],qA[j-1],qA[j],qA[j+1],qA[j+2],dx)
        
        P1r0a11,P1r0a10  =  P1jtojp1(qA[j],qA[j+1],dx)
        P1r1a11,P1r1a10  =  P1jm1toj(qA[j],qA[j-1],dx)
        
        P2ca22,P2ca21,P2ca20  =  P2jm1tojp1(qA[j-1],qA[j],qA[j+1],dx)

        # w1 = minmodL([P1r0a11,P1r1a11,P1Wa11])  / P1Wa11
        # P1a11 = P1Wa11
        P1a11w1 = minmodL([P1r0a11,P2ca21,P1r1a11])
        
        #w2 = minmodL([P2r0a22,P2r1a2,P2r2a22])  / P2r1a22
        #P2a22 = P2r1a22
        P2a22w2 = minmodL([P2r0a22,P4ca42,P2r2a22]) 
        
        #w3 = minmodL([P3r0a33,P3r1a33,P3r2a33,P3r3a33])  / P3r2a33
        #P3a33 = P3r2a33
        P3a33w3 = minmodL([P3r0a33,P6ca63,P3r3a33])
        
        #Weight on second coefficient
        if ( abs(P6ca63) > 10.0**(-12.0)):
            w3 = P3a33w3  / P6ca63
        else:
            w3 = 0
        #P3a32w3 = w3*P3r2a32
        P3a32w3Corr = w3*(P6ca62 - P2a22w2) + P2a22w2
        
        
        #Weight on first coefficient, want it to reduce to Reconed P2a22
        if (abs(P4ca42)  > 10.0**(-12.0)):
            w2 = P2a22w2 / P4ca42
        else:
            w2 = 0
        P2a21w2 = w2 *P4ca41
        
        P2a21w2Corr = w2*(P4ca41 - P1a11w1) + P1a11w1
        
        P3a31w3Corr = w3*(P6ca61 - P2a21w2Corr) + P2a21w2Corr
        
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
    return a3,a2,a1,a0

def PolyFromRep(ujmh,ujph,dujmh,dujph,dx):
    a3 = (dx*dujph + dx*dujmh - 2*ujph + 2*ujmh)/dx**3
    a2 = (dujph- dujmh)/(2*dx)
    a1 = (-dx*dujph - dx*dujmh + 6*ujph - 6*ujmh)/(4*dx)
    a0 = -dx*dujph/8 + dx*dujmh/8 + ujph/2 + ujmh/2
    
    return a3,a2,a1,a0

def GetCellAverage(n,q,dx):
    
    qn = zeros(n)
    # cubic passes through y1, y2, y3, y4
    # pa = -9 Gjmh + 27 Gms - 27 Gps + 9 Gph / 2*dx3
    # pb = 9 Gjmh - 9 Gms - 9 Gps + 9 Gph / 4*dx2
    # pc = Gjmh - 27 Gms + 27 Gps - Gph / 8*dx
    # pd = -Gjmh + 9 Gms + 9 Gps - Gph / 16
    for i in range(n):
        qjmh = q[4*i]
        qjms = q[4*i + 1]
        qjps = q[4*i + 2]
        qjph = q[4*i + 3]
        #pa = (-9*qjmh + 27*qjms - 27*qjps  + 9*qjph )/ (2*dx**3)
        pb = (9*qjmh - 9*qjms - 9*qjps  + 9*qjph )/ (4*dx**2)
        #pc = (qjmh - 27*qjms + 27*qjps  - qjph )/ (8*dx)
        pd = (-qjmh + 9*qjms + 9*qjps  - qjph )/ 16.0

        #QTop = pa/4* (dx/2)**4 + pb/3* (dx/2)**3 + pc/2* (dx/2)**2 +  pd*(dx/2) 
        #QBot = pa/4* (dx/2)**4 - pb/3* (dx/2)**3 + pc/2* (dx/2)**2 -  pd*(dx/2) 
        
        qn[i] = 1.0/dx*(2*pb/3* (dx/2)**3 +  2*pd*(dx/2) )
    return qn


def InitialSolitonPointwise(x,dx,t,ga,a0,a1,b1):
    n = len(x)
    hp = []
    dhp = []
    ddhp = []
    
    Gp = []
    up = []
    dup = []
    ddup = []
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
        dhjmh = -2*a1*k*tanh(k*phi)*sechkphi**2
        ddhjmh = -2*a1*k**2*sechkphi**2*(sechkphi**2 -2*tanh(k*phi)**2 )
        ujmh = c*(1 - a0 /  hjmh)
        dujmh = a0*c*dhjmh/hjmh**2
        ddujmh = a0*c*(hjmh*ddhjmh-2*dhjmh**2)/hjmh**3
        Gjmh = ujmh*hjmh + b1*a0*a1*c*(k**2)*sechkphi**4*hjmh - 2*b1*a0*a1**2*c*k**2*sechkphi**4*tanh(k*phi)**2  - 2*b1*a0*a1*c*k**2*sechkphi**2*hjmh*tanh(k*phi)**2 

        phi = xjms - c*t
        sechkphi = 1.0 / cosh(k*phi)
        hjms = a0 + a1*sechkphi**2
        dhjms = -2*a1*k*tanh(k*phi)*sechkphi**2
        ddhjms = -2*a1*k**2*sechkphi**2*(sechkphi**2 -2*tanh(k*phi)**2 )
        ujms = c*(1 - a0 /  hjms)
        dujms = a0*c*dhjms/hjms**2
        ddujms = a0*c*(hjms*ddhjms-2*dhjms**2)/hjms**3
        Gjms = ujms*hjms + b1*a0*a1*c*(k**2)*sechkphi**4*hjms - 2*b1*a0*a1**2*c*k**2*sechkphi**4*tanh(k*phi)**2  - 2*b1*a0*a1*c*k**2*sechkphi**2*hjms*tanh(k*phi)**2       
        

        phi = xjps - c*t
        sechkphi = 1.0 / cosh(k*phi)
        hjps = a0 + a1*sechkphi**2
        dhjps = -2*a1*k*tanh(k*phi)*sechkphi**2
        ddhjps = -2*a1*k**2*sechkphi**2*(sechkphi**2 -2*tanh(k*phi)**2 )
        ujps = c*(1 - a0 /  hjps)
        dujps = a0*c*dhjps/hjps**2
        ddujps = a0*c*(hjps*ddhjps-2*dhjps**2)/hjps**3
        Gjps = ujps*hjps + b1*a0*a1*c*(k**2)*sechkphi**4*hjps - 2*b1*a0*a1**2*c*k**2*sechkphi**4*tanh(k*phi)**2  - 2*b1*a0*a1*c*k**2*sechkphi**2*hjps*tanh(k*phi)**2  
 
        phi = xjph - c*t
        sechkphi = 1.0 / cosh(k*phi)
        hjph = a0 + a1*sechkphi**2
        dhjph = -2*a1*k*tanh(k*phi)*sechkphi**2
        ddhjph = -2*a1*k**2*sechkphi**2*(sechkphi**2 -2*tanh(k*phi)**2 )
        ujph = c*(1 - a0 /  hjph)
        dujph = a0*c*dhjph/hjph**2
        ddujph = a0*c*(hjph*ddhjph-2*dhjph**2)/hjph**3
        Gjph = ujph*hjph + b1*a0*a1*c*(k**2)*sechkphi**4*hjph - 2*b1*a0*a1**2*c*k**2*sechkphi**4*tanh(k*phi)**2  - 2*b1*a0*a1*c*k**2*sechkphi**2*hjph*tanh(k*phi)**2  
       

        hcell = [hjmh,hjms,hjps,hjph]
        dhcell = [dhjmh,dhjms,dhjps,dhjph]
        ddhcell = [ddhjmh,ddhjms,ddhjps,ddhjph]
        ucell = [ujmh,ujms,ujps,ujph]
        ducell = [dujmh,dujms,dujps,dujph]
        dducell = [ddujmh,ddujms,ddujps,ddujph]
        Gcell = [Gjmh,Gjms,Gjps,Gjph]
        xhcell = [xjmh,xjms,xjps,xjph]
        hp = hp + hcell
        dhp = dhp + dhcell
        ddhp = ddhp + ddhcell
        xp = xp + xhcell
        up = up + ucell
        dup = dup + ducell
        ddup = ddup + dducell
        Gp = Gp + Gcell
        
        
    return array(hp),array(Gp),array(up),array(xp),array(dhp),array(ddhp),array(dup),array(ddup)



def FEMElem(hjmh,hjms,hjps,hjph,Gjmh,Gjms,Gjps,Gjph,beta1,dx):
    
    Ge = zeros((4,4))
    Ge[0,0] = dx/2*(16.0/105)
    Ge[0,1] = dx/2*(33.0/280)
    Ge[0,2] = dx/2*(-3.0/70)
    Ge[0,3] = dx/2*(19.0/840)
    Ge[1,0] = dx/2*(33.0/280)
    Ge[1,1] = dx/2*(27.0/35)
    Ge[1,2] = dx/2*(-27.0/280)
    Ge[1,3] = dx/2*(-3.0/70)
    Ge[2,0] = dx/2*(-3.0/70)
    Ge[2,1] = dx/2*(-27.0/280)
    Ge[2,2] = dx/2*(27.0/35)
    Ge[2,3] = dx/2*(33.0/280)
    Ge[3,0] = dx/2*(19.0/840)
    Ge[3,1] = dx/2*(-3.0/70)
    Ge[3,2] = dx/2*(33.0/280)
    Ge[3,3] = dx/2*(16.0/105)
    
    GMatE = dot(Ge,array([Gjmh,Gjms,Gjps,Gjph]))


    
    uhe = zeros((4,4))
    
    uhe[0,0] = dx/2*(17*hjmh/160 + 9*hjms/140 + hjph/168 - 27*hjps/1120)
    uhe[0,1] = dx/2*(9*hjmh/140 + 27*hjms/280 + 3*hjph/560 - 27*hjps/560)
    uhe[0,2] = dx/2*(-27*hjmh/1120 - 27*hjms/560 + 3*hjph/560 + 27*hjps/1120)
    uhe[0,3] = dx/2*(hjmh/168 + 3*hjms/560 + hjph/168 + 3*hjps/560)
    uhe[1,0] = dx/2*(9*hjmh/140 + 27*hjms/280 + 3*hjph/560 - 27*hjps/560)
    uhe[1,1] = dx/2*(27*hjmh/280 + 729*hjms/1120 + 27*hjph/1120)
    uhe[1,2] = dx/2*(-27*hjmh/560 - 27*hjph/560)
    uhe[1,3] = dx/2*(3*hjmh/560 + 27*hjms/1120 - 27*hjph/1120 - 27*hjps/560)
    uhe[2,0] = dx/2*(-27*hjmh/1120 - 27*hjms/560 + 3*hjph/560 + 27*hjps/1120)
    uhe[2,1] = dx/2*(-27*hjmh/560 - 27*hjph/560)
    uhe[2,2] = dx/2*(27*hjmh/1120 + 27*hjph/280 + 729*hjps/1120)
    uhe[2,3] = dx/2*(3*hjmh/560 - 27*hjms/560 + 9*hjph/140 + 27*hjps/280)
    uhe[3,0] = dx/2*(hjmh/168 + 3*hjms/560 + hjph/168 + 3*hjps/560)
    uhe[3,1] = dx/2*(3*hjmh/560 + 27*hjms/1120 - 27*hjph/1120 - 27*hjps/560)
    uhe[3,2] = dx/2*(3*hjmh/560 - 27*hjms/560 + 9*hjph/140 + 27*hjps/280)
    uhe[3,3] = dx/2*(hjmh/168 - 27*hjms/1120 + 17*hjph/160 + 9*hjps/140)
    

    
    
    h3uxe = zeros((4,4))
    h3uxe[0,0] = (beta1/2.0)*2/dx*(5377*hjmh**3/8960 + 490239*hjmh**2*hjms/640640 + 45223*hjmh**2*hjph/640640 - 7695*hjmh**2*hjps/23296 + 19035*hjmh*hjms**2/23296 + 2601*hjmh*hjms*hjph/20020 - 100521*hjmh*hjms*hjps/160160 + 1163*hjmh*hjph**2/183040 - 8109*hjmh*hjph*hjps/160160 + 14499*hjmh*hjps**2/116480 + 235467*hjms**3/366080 + 297837*hjms**2*hjph/2562560 - 728271*hjms**2*hjps/1281280 + 9963*hjms*hjph**2/2562560 - 1377*hjms*hjph*hjps/16016 + 124659*hjms*hjps**2/640640 + 953*hjph**3/73216 + 8397*hjph**2*hjps/1281280 + 13527*hjph*hjps**2/640640 + 729*hjps**3/1281280)
    h3uxe[0,1] = (beta1/2.0)*2/dx*(-1164861*hjmh**3/1281280 - 2657367*hjmh**2*hjms/2562560 - 251253*hjmh**2*hjph/2562560 + 584091*hjmh**2*hjps/1281280 - 2374353*hjmh*hjms**2/2562560 - 39609*hjmh*hjms*hjph/256256 + 963009*hjmh*hjms*hjps/1281280 - 28737*hjmh*hjph**2/2562560 + 75087*hjmh*hjph*hjps/1281280 - 197559*hjmh*hjps**2/1281280 - 518319*hjms**3/1281280 - 21141*hjms**2*hjph/183040 + 269001*hjms**2*hjps/512512 + 1377*hjms*hjph**2/116480 + 143613*hjms*hjph*hjps/1281280 - 518319*hjms*hjps**2/2562560 - 6939*hjph**3/116480 - 8343*hjph**2*hjps/197120 - 114453*hjph*hjps**2/2562560 - 150903*hjps**3/1281280)
    h3uxe[0,2] = (beta1/2.0)*2/dx*(503469*hjmh**3/1281280 + 108783*hjmh**2*hjms/320320 + 171*hjmh**2*hjph/4928 - 201771*hjmh**2*hjps/1281280 + 165483*hjmh*hjms**2/1281280 + 2187*hjmh*hjms*hjph/80080 - 729*hjmh*hjms*hjps/4928 + 15507*hjmh*hjph**2/1281280 - 243*hjmh*hjph*hjps/45760 + 43011*hjmh*hjps**2/1281280 - 59049*hjms**3/197120 + 729*hjms**2*hjph/232960 + 6561*hjms**2*hjps/116480 - 11097*hjms*hjph**2/232960 - 729*hjms*hjph*hjps/14560 + 6561*hjms*hjps**2/320320 + 47763*hjph**3/366080 + 12069*hjph**2*hjps/116480 + 13851*hjph*hjps**2/320320 + 6561*hjps**3/116480)
    h3uxe[0,3] = (beta1/2.0)*2/dx*(-107519*hjmh**3/1281280 - 173853*hjmh**2*hjms/2562560 - 18559*hjmh**2*hjph/2562560 + 8181*hjmh**2*hjps/256256 - 7209*hjmh*hjms**2/366080 - 3411*hjmh*hjms*hjph/1281280 + 30699*hjmh*hjms*hjps/1281280 - 18559*hjmh*hjph**2/2562560 - 3411*hjmh*hjph*hjps/1281280 - 4941*hjmh*hjps**2/1281280 + 78003*hjms**3/1281280 - 4941*hjms**2*hjph/1281280 - 6561*hjms**2*hjps/512512 + 8181*hjms*hjph**2/256256 + 30699*hjms*hjph*hjps/1281280 - 6561*hjms*hjps**2/512512 - 107519*hjph**3/1281280 - 173853*hjph**2*hjps/2562560 - 7209*hjph*hjps**2/366080 + 78003*hjps**3/1281280)
    h3uxe[1,0] = (beta1/2.0)*2/dx*(-1164861*hjmh**3/1281280 - 2657367*hjmh**2*hjms/2562560 - 251253*hjmh**2*hjph/2562560 + 584091*hjmh**2*hjps/1281280 - 2374353*hjmh*hjms**2/2562560 - 39609*hjmh*hjms*hjph/256256 + 963009*hjmh*hjms*hjps/1281280 - 28737*hjmh*hjph**2/2562560 + 75087*hjmh*hjph*hjps/1281280 - 197559*hjmh*hjps**2/1281280 - 518319*hjms**3/1281280 - 21141*hjms**2*hjph/183040 + 269001*hjms**2*hjps/512512 + 1377*hjms*hjph**2/116480 + 143613*hjms*hjph*hjps/1281280 - 518319*hjms*hjps**2/2562560 - 6939*hjph**3/116480 - 8343*hjph**2*hjps/197120 - 114453*hjph*hjps**2/2562560 - 150903*hjps**3/1281280)
    h3uxe[1,1] = (beta1/2.0)*2/dx*(323217*hjmh**3/232960 + 1858221*hjmh**2*hjms/1281280 + 13689*hjmh**2*hjph/98560 - 1630773*hjmh**2*hjps/2562560 + 19683*hjmh*hjms**2/18304 + 729*hjmh*hjms*hjph/3640 - 6561*hjmh*hjms*hjps/6160 + 18873*hjmh*hjph**2/640640 - 729*hjmh*hjph*hjps/16016 + 400221*hjmh*hjps**2/2562560 + 85293*hjms**3/98560 + 111537*hjms**2*hjph/1281280 - 19683*hjms**2*hjps/640640 - 21141*hjms*hjph**2/256256 - 6561*hjms*hjph*hjps/20020 + 846369*hjms*hjps**2/1281280 + 50409*hjph**3/183040 + 140697*hjph**2*hjps/640640 + 19683*hjph*hjps**2/183040 + 478953*hjps**3/512512)
    h3uxe[1,2] = (beta1/2.0)*2/dx*(-111429*hjmh**3/183040 - 1324593*hjmh**2*hjms/2562560 - 27135*hjmh**2*hjph/512512 + 292329*hjmh**2*hjps/1281280 - 98415*hjmh*hjms**2/512512 - 51759*hjmh*hjms*hjph/1281280 + 465831*hjmh*hjms*hjps/1281280 - 27135*hjmh*hjph**2/512512 - 51759*hjmh*hjph*hjps/1281280 - 6561*hjmh*hjps**2/1281280 - 662661*hjms**3/1281280 - 6561*hjms**2*hjph/1281280 - 1318761*hjms**2*hjps/2562560 + 292329*hjms*hjph**2/1281280 + 465831*hjms*hjph*hjps/1281280 - 1318761*hjms*hjps**2/2562560 - 111429*hjph**3/183040 - 1324593*hjph**2*hjps/2562560 - 98415*hjph*hjps**2/512512 - 662661*hjps**3/1281280)
    h3uxe[1,3] = (beta1/2.0)*2/dx*(47763*hjmh**3/366080 + 12069*hjmh**2*hjms/116480 + 15507*hjmh**2*hjph/1281280 - 11097*hjmh**2*hjps/232960 + 13851*hjmh*hjms**2/320320 - 243*hjmh*hjms*hjph/45760 - 729*hjmh*hjms*hjps/14560 + 171*hjmh*hjph**2/4928 + 2187*hjmh*hjph*hjps/80080 + 729*hjmh*hjps**2/232960 + 6561*hjms**3/116480 + 43011*hjms**2*hjph/1281280 + 6561*hjms**2*hjps/320320 - 201771*hjms*hjph**2/1281280 - 729*hjms*hjph*hjps/4928 + 6561*hjms*hjps**2/116480 + 503469*hjph**3/1281280 + 108783*hjph**2*hjps/320320 + 165483*hjph*hjps**2/1281280 - 59049*hjps**3/197120)
    h3uxe[2,0] = (beta1/2.0)*2/dx*(503469*hjmh**3/1281280 + 108783*hjmh**2*hjms/320320 + 171*hjmh**2*hjph/4928 - 201771*hjmh**2*hjps/1281280 + 165483*hjmh*hjms**2/1281280 + 2187*hjmh*hjms*hjph/80080 - 729*hjmh*hjms*hjps/4928 + 15507*hjmh*hjph**2/1281280 - 243*hjmh*hjph*hjps/45760 + 43011*hjmh*hjps**2/1281280 - 59049*hjms**3/197120 + 729*hjms**2*hjph/232960 + 6561*hjms**2*hjps/116480 - 11097*hjms*hjph**2/232960 - 729*hjms*hjph*hjps/14560 + 6561*hjms*hjps**2/320320 + 47763*hjph**3/366080 + 12069*hjph**2*hjps/116480 + 13851*hjph*hjps**2/320320 + 6561*hjps**3/116480)
    h3uxe[2,1] = (beta1/2.0)*2/dx*(-111429*hjmh**3/183040 - 1324593*hjmh**2*hjms/2562560 - 27135*hjmh**2*hjph/512512 + 292329*hjmh**2*hjps/1281280 - 98415*hjmh*hjms**2/512512 - 51759*hjmh*hjms*hjph/1281280 + 465831*hjmh*hjms*hjps/1281280 - 27135*hjmh*hjph**2/512512 - 51759*hjmh*hjph*hjps/1281280 - 6561*hjmh*hjps**2/1281280 - 662661*hjms**3/1281280 - 6561*hjms**2*hjph/1281280 - 1318761*hjms**2*hjps/2562560 + 292329*hjms*hjph**2/1281280 + 465831*hjms*hjph*hjps/1281280 - 1318761*hjms*hjps**2/2562560 - 111429*hjph**3/183040 - 1324593*hjph**2*hjps/2562560 - 98415*hjph*hjps**2/512512 - 662661*hjps**3/1281280)
    h3uxe[2,2] = (beta1/2.0)*2/dx*(50409*hjmh**3/183040 + 140697*hjmh**2*hjms/640640 + 18873*hjmh**2*hjph/640640 - 21141*hjmh**2*hjps/256256 + 19683*hjmh*hjms**2/183040 - 729*hjmh*hjms*hjph/16016 - 6561*hjmh*hjms*hjps/20020 + 13689*hjmh*hjph**2/98560 + 729*hjmh*hjph*hjps/3640 + 111537*hjmh*hjps**2/1281280 + 478953*hjms**3/512512 + 400221*hjms**2*hjph/2562560 + 846369*hjms**2*hjps/1281280 - 1630773*hjms*hjph**2/2562560 - 6561*hjms*hjph*hjps/6160 - 19683*hjms*hjps**2/640640 + 323217*hjph**3/232960 + 1858221*hjph**2*hjps/1281280 + 19683*hjph*hjps**2/18304 + 85293*hjps**3/98560)
    h3uxe[2,3] = (beta1/2.0)*2/dx*(-6939*hjmh**3/116480 - 8343*hjmh**2*hjms/197120 - 28737*hjmh**2*hjph/2562560 + 1377*hjmh**2*hjps/116480 - 114453*hjmh*hjms**2/2562560 + 75087*hjmh*hjms*hjph/1281280 + 143613*hjmh*hjms*hjps/1281280 - 251253*hjmh*hjph**2/2562560 - 39609*hjmh*hjph*hjps/256256 - 21141*hjmh*hjps**2/183040 - 150903*hjms**3/1281280 - 197559*hjms**2*hjph/1281280 - 518319*hjms**2*hjps/2562560 + 584091*hjms*hjph**2/1281280 + 963009*hjms*hjph*hjps/1281280 + 269001*hjms*hjps**2/512512 - 1164861*hjph**3/1281280 - 2657367*hjph**2*hjps/2562560 - 2374353*hjph*hjps**2/2562560 - 518319*hjps**3/1281280)
    h3uxe[3,0] = (beta1/2.0)*2/dx*(-107519*hjmh**3/1281280 - 173853*hjmh**2*hjms/2562560 - 18559*hjmh**2*hjph/2562560 + 8181*hjmh**2*hjps/256256 - 7209*hjmh*hjms**2/366080 - 3411*hjmh*hjms*hjph/1281280 + 30699*hjmh*hjms*hjps/1281280 - 18559*hjmh*hjph**2/2562560 - 3411*hjmh*hjph*hjps/1281280 - 4941*hjmh*hjps**2/1281280 + 78003*hjms**3/1281280 - 4941*hjms**2*hjph/1281280 - 6561*hjms**2*hjps/512512 + 8181*hjms*hjph**2/256256 + 30699*hjms*hjph*hjps/1281280 - 6561*hjms*hjps**2/512512 - 107519*hjph**3/1281280 - 173853*hjph**2*hjps/2562560 - 7209*hjph*hjps**2/366080 + 78003*hjps**3/1281280)
    h3uxe[3,1] = (beta1/2.0)*2/dx*(47763*hjmh**3/366080 + 12069*hjmh**2*hjms/116480 + 15507*hjmh**2*hjph/1281280 - 11097*hjmh**2*hjps/232960 + 13851*hjmh*hjms**2/320320 - 243*hjmh*hjms*hjph/45760 - 729*hjmh*hjms*hjps/14560 + 171*hjmh*hjph**2/4928 + 2187*hjmh*hjph*hjps/80080 + 729*hjmh*hjps**2/232960 + 6561*hjms**3/116480 + 43011*hjms**2*hjph/1281280 + 6561*hjms**2*hjps/320320 - 201771*hjms*hjph**2/1281280 - 729*hjms*hjph*hjps/4928 + 6561*hjms*hjps**2/116480 + 503469*hjph**3/1281280 + 108783*hjph**2*hjps/320320 + 165483*hjph*hjps**2/1281280 - 59049*hjps**3/197120)
    h3uxe[3,2] = (beta1/2.0)*2/dx*(-6939*hjmh**3/116480 - 8343*hjmh**2*hjms/197120 - 28737*hjmh**2*hjph/2562560 + 1377*hjmh**2*hjps/116480 - 114453*hjmh*hjms**2/2562560 + 75087*hjmh*hjms*hjph/1281280 + 143613*hjmh*hjms*hjps/1281280 - 251253*hjmh*hjph**2/2562560 - 39609*hjmh*hjph*hjps/256256 - 21141*hjmh*hjps**2/183040 - 150903*hjms**3/1281280 - 197559*hjms**2*hjph/1281280 - 518319*hjms**2*hjps/2562560 + 584091*hjms*hjph**2/1281280 + 963009*hjms*hjph*hjps/1281280 + 269001*hjms*hjps**2/512512 - 1164861*hjph**3/1281280 - 2657367*hjph**2*hjps/2562560 - 2374353*hjph*hjps**2/2562560 - 518319*hjps**3/1281280)
    h3uxe[3,3] = (beta1/2.0)*2/dx*(953*hjmh**3/73216 + 8397*hjmh**2*hjms/1281280 + 1163*hjmh**2*hjph/183040 + 9963*hjmh**2*hjps/2562560 + 13527*hjmh*hjms**2/640640 - 8109*hjmh*hjms*hjph/160160 - 1377*hjmh*hjms*hjps/16016 + 45223*hjmh*hjph**2/640640 + 2601*hjmh*hjph*hjps/20020 + 297837*hjmh*hjps**2/2562560 + 729*hjms**3/1281280 + 14499*hjms**2*hjph/116480 + 124659*hjms**2*hjps/640640 - 7695*hjms*hjph**2/23296 - 100521*hjms*hjph*hjps/160160 - 728271*hjms*hjps**2/1281280 + 5377*hjph**3/8960 + 490239*hjph**2*hjps/640640 + 19035*hjph*hjps**2/23296 + 235467*hjps**3/366080)   

    AMatE = uhe +  h3uxe
    
    return GMatE,AMatE


def FEM(hp,Gp,u0,u1,beta1,n,dx,nG):
    
    nb = 3*(n) + 1
    nbG = 3*nG + 1
    A = zeros((nb,nb))
    b = zeros(nb)
    
    # print(n,nbG,nb)

    
    for i in range(nG,n-nG):
    #    #elementwisematrices
        #j 
        hjmh = hp[4*i]
        hjms = hp[4*i+1]
        hjps = hp[4*i+2]
        hjph = hp[4*i+3]
        
        Gjmh = Gp[4*i]
        Gjms = Gp[4*i+1]
        Gjps = Gp[4*i+2]
        Gjph = Gp[4*i+3]   
        GMatE,AMatE = FEMElem(hjmh,hjms,hjps,hjph,Gjmh,Gjms,Gjps,Gjph,beta1,dx)
        
        # print(GMatE,AMatE)
        # print(nG, 3*i,3*i+4, [*range(3*i,3*i+4)])
        A[3*i  : 3*i+ 4,3*i : 3*i+ 4] = A[3*i  : 3*i+ 4,3*i : 3*i+ 4] + AMatE
        b[3*i  : 3*i+ 4]=  b[3*i : 3*i+ 4] + GMatE

 
    # # i = 0
    i=0
    A[nbG-1,nbG-1:nbG-1 + nbG] = zeros(nbG)
    A[nbG-1:nbG-1 + nbG,nbG-1] = zeros(nbG)
    A[0 : nbG,0 : nbG ] =  eye(nbG,nbG)
    b[0 : nbG ]=  u0*ones(nbG)
 
    # # i = 0
    i=nb-nbG
    A[nb-nbG,nb-nbG -nbG: nb-nbG] = zeros(nbG)
    A[nb-nbG -nbG: nb-nbG,nb-nbG] = zeros(nbG)
    A[nb-nbG : nb,nb-nbG : nb ] =  eye(nbG,nbG)
    b[nb-nbG: nb ]=  u0*ones(nbG)
    

    uN = solve(A,b)
    
    # print(A)
    # print(b)
    return A,b,uN 


def ContEdgesFromPW(q,xn,dx):
    n = len(xn)
    nq = zeros(3*n+1)
    nx = zeros(3*n+1)
    
    nq[0] = q[0]
    nx[0] = xn[0] - 0.5*dx
    for i in range(n):
        nq[3*i + 1] = q[4*i + 1]
        nq[3*i + 2] = q[4*i + 2]
        nq[3*i + 3] = q[4*i + 3]

        nx[3*i + 1] = xn[i] - dx/6.0
        nx[3*i + 2] = xn[i] + dx/6.0
        nx[3*i + 3] = xn[i] + dx/2.0  
    return nq,nx

def ContEdgesToPW(q,xn,dx):
    n = len(xn)
    nq = zeros(4*n)   
    for i in range(n):
        nq[4*i] = q[3*i]
        nq[4*i+1] = q[3*i+1]
        nq[4*i+2] = q[3*i+2]
        nq[4*i + 3] = q[3*i+3]
    return nq

def Derivatives(q,xn,dx):
    n = len(xn)
    dq = zeros(4*n)
    ddq = zeros(4*n)
    for i in range(n):
        qjmh = q[4*i]
        qjms = q[4*i+1]
        qjps = q[4*i+2]
        qjph = q[4*i+3]
        
        a3,a2,a1,a0 = PolyFromPoints(qjmh,qjms,qjps,qjph,dx)
        
        dqjmh = PointsFromPolyFromPoints(a1,2*a2,3*a3,0,-dx/2)
        dqjms = PointsFromPolyFromPoints(a1,2*a2,3*a3,0,-dx/6)
        dqjps = PointsFromPolyFromPoints(a1,2*a2,3*a3,0,dx/6)
        dqjph = PointsFromPolyFromPoints(a1,2*a2,3*a3,0,dx/2)
        
        dq[4*i] = dqjmh
        dq[4*i+1] = dqjms
        dq[4*i+2] = dqjps
        dq[4*i+3] = dqjph

        ddqjmh = PointsFromPolyFromPoints(2*a2,6*a3,0,0,-dx/2)
        ddqjms = PointsFromPolyFromPoints(2*a2,6*a3,0,0,-dx/6)
        ddqjps = PointsFromPolyFromPoints(2*a2,6*a3,0,0,dx/6)
        ddqjph = PointsFromPolyFromPoints(2*a2,6*a3,0,0,dx/2)

        ddq[4*i] = ddqjmh
        ddq[4*i+1] = ddqjms
        ddq[4*i+2] = ddqjps
        ddq[4*i+3] = ddqjph
    return dq,ddq
    


L2GNs = []
L2hNs = []
L2dhNs = []
L2ddhNs = []
L2uNs = []
L2uNNs = []
L2duNNs = []
L2dduNNs = []
L2b1Terms = []
L2b2Terms = []
dxs = []    
# for ij in range(5,6):
for ij in range(10):
    sx = -50
    ex = 50
    n= 10*2**ij
    # n=4
    nG = 6
    dx = (ex- sx)/n
    xn = arange(sx - nG*dx,ex+ (nG+ 1)*dx,dx)
    
    ga = 10.0
    a0 = 1.0
    a1 = 0.7
    t = 0
    # b1 = 0
    b1 = 2.0/3.0
    
    h,G,u,x,dh,ddh,du,ddu = InitialSolitonPointwise(xn,dx,t,ga,a0,a1,b1)
    
    hA = GetCellAverage(len(xn),h,dx)
    GA = GetCellAverage(len(xn),G,dx)
    
    hN = Recon(hA,xn,dx,nG)
    GN = Recon(GA,xn,dx,nG)
    
    hAN = GetCellAverage(len(xn),hN,dx)
    GAN = GetCellAverage(len(xn),GN,dx)
    
    hNN = Recon(hAN,xn,dx,nG)
    GNN = Recon(GAN,xn,dx,nG)    
    
    dhN,ddhN  = Derivatives(hN,xn,dx)
    
    uC,xC = ContEdgesFromPW(u,xn,dx)
    duC,xC = ContEdgesFromPW(du,xn,dx)
    
    
    u0 = 0
    u1 = 0
    # A,b = FEM(h,G,u0,u1,b1,len(xn),dx,nG )
    # A,b,uN = FEM(h,G,u0,u1,b1,len(xn),dx,nG )
    A,b,uNN = FEM(hN,GN,u0,u1,b1,len(xn),dx,nG )
    
    
    PWuNN = ContEdgesToPW(uNN,xn,dx)
    
    duNN,dduNN  = Derivatives(PWuNN,xn,dx)
    
    
    L2hN = norm(h - hN,ord=2)/norm(h,ord=2)
    L2dhN = norm(dhN - dh,ord=2)/norm(dh,ord=2)
    L2ddhN = norm(ddhN - ddh,ord=2)/norm(ddh,ord=2)
    
    L2GN = norm(G - GN,ord=2)/norm(G,ord=2)
    # L2uN = norm(uC - uN,ord=2)/norm(uC,ord=2)
    L2uNN = norm(u - PWuNN,ord=2)/norm(u,ord=2)
    L2duNN = norm(duNN - du,ord=2)/norm(du,ord=2)
    L2dduNN = norm(dduNN - ddu,ord=2)/norm(ddu,ord=2)
    
    L2b1Term = norm(b1*hN**3*duNN**2 - b1*h**3*du**2,ord=2 )/ norm(b1*h**3*du**2,ord=2 )
    # L2b1Term = norm(b1*hN**3*duNN**2 - b1*h**3*du**2,ord=2 )/ norm(b1*h**3*du**2,ord=2 )
    # L2b2Term = norm(b1/4*ga*hN**2*(2*hN*ddhN + dhN**2) - b1/4*ga*h**2*(2*h*ddh + dh**2),ord=2)  / norm(b1/4*ga*h**2*(2*h*ddh + dh**2),ord=2) 
    L2b2Term = norm(b1/4*ga*hN**2*(2*hN*ddhN**2 + dhN**2) - b1/4*ga*h**2*(2*h*ddh**2 + dh**2),ord=2)  / norm(b1/4*ga*h**2*(2*h*ddh**2 + dh**2),ord=2) 
    
    L2GNs.append(L2GN)
    L2hNs.append(L2hN)
    L2dhNs.append(L2dhN)
    L2ddhNs.append(L2ddhN)
    # L2uNs.append(L2uN)
    L2uNNs.append(L2uNN)
    L2duNNs.append(L2duNN)
    L2dduNNs.append(L2dduNN)
    dxs.append(dx)
    
    L2b1Terms.append(L2b1Term)
    L2b2Terms.append(L2b2Term)
    
# loglog(dxs,L2uNs,'.b')
loglog(dxs,L2uNNs,'.b')
loglog(dxs,L2duNNs,'.r')
loglog(dxs,L2dduNNs,'.g')
loglog(dxs,L2hNs,'+b')
loglog(dxs,L2dhNs,'+r')
loglog(dxs,L2ddhNs,'+g')
loglog(dxs,L2GNs,'*b')

loglog(dxs,L2b1Terms,'sb')
loglog(dxs,L2b2Terms,'sr')
loglog(dxs,(L2ddhNs[4]/dxs[4]**2) *array(dxs)**2)
loglog(dxs,(L2ddhNs[4]/dxs[4]**3) *array(dxs)**3)
loglog(dxs,(L2ddhNs[4]/dxs[4]**4) *array(dxs)**4)
# # # j= 10
# j = n//2
# be,Ae = FEMElem(Nh,NG,b1,dx,j)

# Nucell= Nu[4*j:4*(j+1)]
# NGcells= NG[4*(j-1):4*(j+2)]

# Gcell = dot(Ae,Nucell)

# Gs = dot(A,NuC)

# bs = lstsq(A,Gs,rcond=0.0001)


# Ga  = GetCellAverage(len(xn),G,dx)
# # un = FEM(b1,h,Ga,len(xn),x,dx)

# A,b = FEMforG(xn,h,Ga,u[0],u[-1],nG,dx)

# i = n//2
# ha0,ha1,ha2,ha3 = PolyFromPoints(h[4*i],h[4*i+1],h[4*i+2],h[4*i+3],dx)
# ua0,ua1,ua2,ua3 = PolyFromPoints(u[4*i],u[4*i+1],u[4*i+2],u[4*i+3],dx)
# Gjmh, Gjms, Gjps,  Gjph = CalculateG(b1,ha0,ha1,ha2,ha3,ua0,ua1,ua2,ua3 ,dx)
# print(Gjmh, Gjms, Gjps,  Gjph)
# print(G[4*i:4*(i+1)])
# G = FEMforG(xn,h,u,dx)