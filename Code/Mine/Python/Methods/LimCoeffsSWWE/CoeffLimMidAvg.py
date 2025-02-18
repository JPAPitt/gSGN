
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from numpy.linalg import solve,norm
from matplotlib.pyplot import plot,loglog

def minmod(a,b,c):
    if a>0 and b>0 and c>0:
        return min(a,b,c)
    if a<0 and b<0 and c<0:
        return max(a,b,c)
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


def Reconjmh(qajm2,qajm1,qaj,qajp1,qajp2,theta):
    
    a0,b0,c0 =  P2jtojp2(qaj,qajp1,qajp2,dx)
    a1,b1,c1 =  P2jm1tojp1(qajm1,qaj,qajp1,dx)
    a2,b2,c2 =  P2jm2toj(qajm2,qajm1,qaj,dx)
    
    # nc = minmod(theta*(c0- qaj),c1- qaj,theta*(c2- qaj)) + qaj
    # nb = minmod(theta*b0,b1,theta*b2)
    na = minmod(theta*a0,a1,theta*a2)
            
    mfd = (qajp1 - qaj)/ dx
    mcd =(qajp1 - qajm1)/ 2*dx
    mbd =(qaj - qajm1)/ dx
    
    mmd = minmod(mfd,mcd,mbd)
    
    mcdd = (qajp1 - 2*qaj + qajm1) / (dx**2)
    mfdd = (-qajp2 + 4*qajp1 - 3*qaj) / (2*dx**2)
    mbdd = (qajm2 - 4*qajm1 + 3*qaj) / (2*dx**2)
    
    mmdd = minmod(mfdd, mcdd,mbdd)
    
    nc = minmod((c0-qaj),c1- qaj,(c2- qaj)) + qaj
    nb = minmod(b0-mmd,b1-mmd,b2-mmd) + mmd
    # na = minmod(a0-mmdd,a1-mmdd,a2-mmdd) + mmdd
    
    qjmh = na*(-dx/2)**2 + nb*(-dx/2) + nc
    
    return qjmh
    
 


def Reconjph(qajm2,qajm1,qaj,qajp1,qajp2,theta):
    
    a0,b0,c0 =  P2jtojp2(qaj,qajp1,qajp2,dx)
    a1,b1,c1 =  P2jm1tojp1(qajm1,qaj,qajp1,dx)
    a2,b2,c2 =  P2jm2toj(qajm2,qajm1,qaj,dx)
    
    nc = minmod(theta*(c0- qaj),c1- qaj,theta*(c2- qaj)) + qaj
    nb = minmod(theta*b0,b1,theta*b2)
    na = minmod(theta*a0,a1,theta*a2)

    mfd = (qajp1 - qaj)/ dx
    mcd =(qajp1 - qajm1)/ 2*dx
    mbd =(qaj - qajm1)/ dx
    
    mmd = minmod(mfd,mcd,mbd)
    
    mcdd = (qajp1 - 2*qaj + qajm1) / (dx**2)
    mfdd = (-qajp2 + 4*qajp1 - 3*qaj) / (2*dx**2)
    mbdd = (qajm2 - 4*qajm1 + 3*qaj) / (2*dx**2)
    
    mmdd = minmod(mfdd, mcdd,mbdd)
    
    nc = minmod((c0-qaj),c1- qaj,(c2- qaj)) + qaj
    nb = minmod(b0-mmd,b1-mmd,b2-mmd) + mmd
    # na = minmod(a0-mmdd,a1-mmdd,a2-mmdd) + mmdd
    
    qjph = na*(dx/2)**2 + nb*(dx/2) + nc
    
    
    return qjph

def FVMSWWE(ha,Ga,nGcells,g,dx,dt):
    
    theta = 1
    n = len(ha)
    hap = zeros(n)
    Gap = zeros(n)
    
    j = nGcells-1
    
    hir = Reconjph(ha[j-2],ha[j-1],ha[j],ha[j+1],ha[j+2], theta)
    hip1l = Reconjmh(ha[j-1],ha[j],ha[j+1],ha[j+2],ha[j+3], theta)
    
    Gir = Reconjph(Ga[j-2],Ga[j-1],Ga[j],Ga[j+1],Ga[j+2], theta)
    Gip1l = Reconjmh(Ga[j-1],Ga[j],Ga[j+1],Ga[j+2],Ga[j+3], theta)
    
    
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
        
        hir = Reconjph(ha[j-2],ha[j-1],ha[j],ha[j+1],ha[j+2], theta)
        hip1l = Reconjmh(ha[j-1],ha[j],ha[j+1],ha[j+2],ha[j+3], theta)
        
        Gir = Reconjph(Ga[j-2],Ga[j-1],Ga[j],Ga[j+1],Ga[j+2], theta)
        Gip1l = Reconjmh(Ga[j-1],Ga[j],Ga[j+1],Ga[j+2],Ga[j+3], theta)
        
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
    
    hA = 0.5*(hApp + hA)
    GA = 0.5*(GApp + GA)
    return hA ,GA
    

def Solver(hA,GA,st,et,g,dx,dt,nGcells):
    
    ct = st
    while ct < st + et:
        
        hA,GA = RKStep(hA,GA,g,dx,nGcells,dt)
        # print(hA[len(hA) //2 -1],hA[len(hA) //2 ], hA[len(hA) //2 +1])
        ct = ct + dt
        print(ct)

    return hA,GA

nGcells = 3

expn = 10
lown = 10

g = 9.81

dxs = []
L2s = []

expi = 7
ncurr = lown*(2**expi)

sx = -50.0
ex = 50.0
dx = ((ex - sx))/(ncurr-1)
#x = arange(sx - nGcells*dx,ex+(nGcells + 0.1)*dx,dx)
x = arange(sx ,ex+(0.1)*dx,dx)
nx = len(x)

hA,GA = initialhDB(x,dx,2,1)


dt = 0.5*dx/sqrt(2*g)

hAN,GAN = Solver(hA,GA,0,4,g,dx,dt,nGcells)

