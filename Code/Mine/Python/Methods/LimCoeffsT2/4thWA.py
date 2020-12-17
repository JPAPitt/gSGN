
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from numpy.linalg import solve,norm
from matplotlib.pyplot import plot,loglog,title,xlabel,ylabel,legend,close


def minmodL(q):
    qmin = min(q)
    qmax = max(q)
    if qmin>0:
        return qmin
    if qmax<0 :
        return qmax
    else:
        return 0

def initialhDB(x,dx,h0,h1):
    n = len(x)
    hA = zeros(n)
    GA = zeros(n)

    for i in range(n):
        
        if x[i] < 0:
            hA[i] = h0
        else:
            hA[i] = h1

        
    return hA,GA


#Linear
def P1jtojp1(qaj,qajp1,dx):
    return (qajp1 - qaj)/ dx ,qaj

def P1jm1toj(qaj,qajm1,dx):
    return (qaj - qajm1)/ dx ,qaj

def P1jWeird(qajp1,qaj,qajm1,dx):
    return (qajp1 - qajm1)/ (2*dx) ,qaj

#Quadratic
def P2jtojp2(qaj,qajp1,qajp2,dx):
    a = (qajp2 - 2*qajp1 + qaj)/ (2*dx**2)
    b = (-qajp2 + 4*qajp1 - 3*qaj)/ (2*dx)
    c =  -qajp2/24 + qajp1/12 + 23*qaj/24
    return a,b,c

def P2jm1tojp1(qajm1,qaj,qajp1,dx):
    a = (qajp1 - 2*qaj + qajm1)/ (2*dx**2)
    b = (qajp1 - qajm1)/ (2*dx)
    c =  13*qaj/12 - qajp1/24 - qajm1/24
    return a,b,c

def P2jm2toj(qajm2,qajm1,qaj,dx):
    a = (qaj - 2*qajm1 + qajm2)/ (2*dx**2)
    b = (qajm2 - 4*qajm1 + 3*qaj)/ (2*dx)
    c =  -qajm2/24 + qajm1/12 + 23*qaj/24
    return a,b,c

#Cubic
def P3jtojp3(qaj,qajp1,qajp2,qajp3,dx):
    a=(3*qajp1 - 3*qajp2 + qajp3 - qaj)/(6*dx**3)
    b=(-5*qajp1 + 4*qajp2 - qajp3 + 2*qaj)/(2*dx**2)
    c=(69*qajp1 - 33*qajp2 + 7*qajp3 - 43*qaj)/(24*dx)
    d=5*qajp1/24 - qajp2/6 + qajp3/24 + 11*qaj/12
    return a,b,c,d


def P3jm1tojp2(qaj,qajm1,qajp1,qajp2,dx):
    a=(-3*qajp1 + qajp2 - qajm1 + 3*qaj)/(6*dx**3)
    b=(qajp1 + qajm1 - 2*qaj)/(2*dx**2)
    c=(27*qajp1 - 5*qajp2 - 7*qajm1 - 15*qaj)/(24*dx)
    d=-qajp1/24 - qajm1/24 + 13*qaj/12
    return a,b,c,d

def P3jm2tojp1(qaj,qajm1,qajm2,qajp1,dx):
    a=(qajp1 + 3*qajm1 - qajm2 - 3*qaj)/(6*dx**3)
    b=(qajp1 + qajm1 - 2*qaj)/(2*dx**2)
    c=(7*qajp1- 27*qajm1 + 5*qajm2 + 15*qaj)/(24*dx)
    d=-qajp1/24 - qajm1/24 + 13*qaj/12
    return a,b,c,d


def P3jm3toj(qaj,qajm1,qajm2,qajm3,dx):
    a=(-3*qajm1 + 3*qajm2 - qajm3 + qaj)/(6*dx**3)
    b=(-5*qajm1 + 4*qajm2 - qajm3 + 2*qaj)/(2*dx**2)
    c=(-69*qajm1 + 33*qajm2 - 7*qajm3 + 43*qaj)/(24*dx)
    d=5*qajm1/24 - qajm2/6 + qajm3/24 + 11*qaj/12
    return a,b,c,d




def Reconjmh(qA,j,theta,dx):
    
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
    na30 = qA[j]  - na32/3*(dx/2)**2 
    
    qjmh =  na33*(-dx/2)**3 +  na32*(-dx/2)**2 +  na31*(-dx/2) +  na30
    
    return qjmh
    
def Reconjph(qA,j,theta,dx):
    
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
    P1a11w11 = minmodL([P1r0a11,P1r1a11])
    
    #w2 = minmodL([P2r0a22,P2r1a2,P2r2a22])  / P2r1a22
    #P2a22 = P2r1a22
    P2a22w22 = minmodL([P2r0a22,P2r1a22,P2r2a22]) 
    P2a21w21 = minmodL([P2r0a21,P2r1a21,P2r2a21]) 
    
    #w3 = minmodL([P3r0a33,P3r1a33,P3r2a33,P3r3a33])  / P3r2a33
    #P3a33 = P3r2a33
    P3a33w33 = minmodL([P3r0a33,P3r1a33,P3r2a33,P3r3a33])
    P3a32w32 = minmodL([P3r0a32,P3r1a32,P3r2a32,P3r3a32])
    P3a31w31 = minmodL([P3r0a31,P3r1a31,P3r2a31,P3r3a31])
    
    #Weight on second coefficient
    if ( abs(P3r2a32) > 10.0**(-12.0)):
        w32 = P3a32w32  / P3r2a32
    else:
        w32 = 0
    #P3a32w3 = w3*P3r2a32
    P3a32w32Corr = P3a32w32 + (1 - w32)*P2a22w22
    
    
    #Weight on first coefficient, want it to reduce to Reconed P2a22
    if (abs(P2r1a21)  > 10.0**(-12.0)):
        w21 = P2a21w21  / P2r1a21
    else:
        w21 = 0

    if (abs(P3r2a31)  > 10.0**(-12.0)):
        w31 = P3a31w31  / P3r2a31
    else:
        w31 = 0
        
    P2a21w2Corr =  P2a21w21  + (1 - w21)*P1a11w11
    
    P3a31w31Corr = P3a31w31 + (1 - w31)*( P2a21w2Corr)
    
    na33 = P3a33w33
    na32 = P3a32w32Corr
    # na32 = P2r1a22
    na31 = P3a31w31Corr
    na30 = qA[j] - na32/3*(dx/2)**2 
    
    qjph =  na33*(dx/2)**3 +  na32*(dx/2)**2 +  na31*(dx/2) +  na30
    
    return qjph

def FVMSWWE(ha,Ga,nGcells,g,dx,dt):
    
    theta = 1
    n = len(ha)
    hap = zeros(n)
    Gap = zeros(n)
    
    j = nGcells-1
    
    # hir = Reconjph(ha[j-2],ha[j-1],ha[j],ha[j+1],ha[j+2], theta)
    # hip1l = Reconjmh(ha[j-1],ha[j],ha[j+1],ha[j+2],ha[j+3], theta)
    hir = Reconjph(ha,j,theta,dx)
    hip1l = Reconjmh(ha,j+1,theta,dx)
    
    Gir = Reconjph(Ga,j,theta,dx)
    Gip1l = Reconjmh(Ga,j+1,theta,dx)

    
    # Gir = Reconjph(Ga[j-2],Ga[j-1],Ga[j],Ga[j+1],Ga[j+2], theta)
    # Gip1l = Reconjmh(Ga[j-1],Ga[j],Ga[j+1],Ga[j+2],Ga[j+3], theta)
    
    
    uir = Gir/ hir
    uip1l = Gip1l/ hip1l
    
    sl = min(0,uir - sqrt(g*hir), uip1l  - sqrt(g*hip1l))
    sr = max(0,uir + sqrt(g*hir), uip1l + sqrt(g*hip1l))
    
    
    felh = Gir
    ferh = Gip1l
    
    felG = uir*Gir + g*hir*hir/2.0
    ferG = uip1l*Gip1l + g*hip1l*hip1l/2.0
    
    if sr == sl:
        isrmsl = 0
    else:
        isrmsl = 1.0/(sr - sl)
    
    foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir))
    foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir))
    

    fih = foh
    fiG = foG
    for j in range(nGcells,n-  nGcells):
        
        # hir = Reconjph(ha[j-2],ha[j-1],ha[j],ha[j+1],ha[j+2], theta)
        # hip1l = Reconjmh(ha[j-1],ha[j],ha[j+1],ha[j+2],ha[j+3], theta)
        
        # Gir = Reconjph(Ga[j-2],Ga[j-1],Ga[j],Ga[j+1],Ga[j+2], theta)
        # Gip1l = Reconjmh(Ga[j-1],Ga[j],Ga[j+1],Ga[j+2],Ga[j+3], theta)

        hir = Reconjph(ha,j,theta,dx)
        hip1l = Reconjmh(ha,j+1,theta,dx)
        
        Gir = Reconjph(Ga,j,theta,dx)
        Gip1l = Reconjmh(Ga,j+1,theta,dx)
        
        uir = Gir/ hir
        uip1l = Gip1l/ hip1l
        
        sl = min(0,uir - sqrt(g*hir), uip1l  - sqrt(g*hip1l))
        sr = max(0,uir + sqrt(g*hir), uip1l + sqrt(g*hip1l))
        
        
        felh = Gir
        ferh = Gip1l
        
        felG = uir*Gir + g*hir*hir/2.0
        ferG = uip1l*Gip1l + g*hip1l*hip1l/2.0
        
        if sr == sl:
            isrmsl = 0
        else:
            isrmsl = 1.0/(sr - sl)
        
        foh = isrmsl*(sr*felh - sl*ferh + sl*sr*(hip1l - hir))
        foG = isrmsl*(sr*felG - sl*ferG + sl*sr*(Gip1l - Gir))


        hap[j] =ha[j] - dt*(foh - fih)/dx
        Gap[j] =Ga[j] - dt*(foG - fiG)/dx
        
        # print(j + 1,hap[j] , ha[j])
        
        fih = foh
        fiG = foG
    
    for j in range(nGcells):
        hap[j] = ha[j]
        hap[n-1 - j]  = ha[n-1 - j]
        Gap[j] = Ga[j]
        Gap[n-1 - j]  = Ga[n-1 - j]
        
    return hap,Gap

def EulerStep(hA,GA,g,dx,nGcells,dt):
    
    #hcubN = CellAverageToCubic(hA,nGcells,dx)
    #GcubN = CellAverageToCubic(GA,nGcells,dx)
    hap,Gap = FVMSWWE(hA,GA,nGcells,g,dx,dt)
    return hap,Gap

def RKStep(hA,GA,g,dx,nGcells,dt):
    
    hAp,GAp = EulerStep(hA,GA,g,dx,nGcells,dt)
    hApp,GApp = EulerStep(hAp,GAp,g,dx,nGcells,dt)
    
    # hA = hA/2 + hApp/2
    # GA = GA/2 + GApp/2

    hAw = 3*hA/4 + hApp/4
    GAw = 3*GA/4 + GApp/4

    hAppp,GAppp = EulerStep(hAw,GAw,g,dx,nGcells,dt)
    
    hA = hA/3 + 2*hAppp/3
    GA = GA/3 + 2*GAppp/3
    return hA ,GA
    

def Solver(hA,GA,st,et,g,dx,dt,nGcells):
    
    ct = st
    while ct < st + et:
        
        hA,GA = RKStep(hA,GA,g,dx,nGcells,dt)
        # print(hA[len(hA) //2 -1],hA[len(hA) //2 ], hA[len(hA) //2 +1])
        ct = ct + dt
        print(ct)

    return hA,GA

nGcells = 4

expn = 10
lown = 10

g = 9.81

dxs = []
L2s = []

expi = 5
ncurr = lown*(2**expi)

sx = -50.0
ex = 50.0
dx = ((ex - sx))/(ncurr-1)
#x = arange(sx - nGcells*dx,ex+(nGcells + 0.1)*dx,dx)
x4 = arange(sx ,ex+(0.1)*dx,dx)
nx = len(x4)

hA,GA = initialhDB(x4,dx,2,1)


dt = 0.5*dx/sqrt(2*g)

hAN4,GAN4 = Solver(hA,GA,0,4,g,dx,dt,nGcells)

