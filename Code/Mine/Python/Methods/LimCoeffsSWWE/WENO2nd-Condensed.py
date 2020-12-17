
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from numpy.linalg import solve,norm
from matplotlib.pyplot import plot,loglog



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





def Reconjmh(qajm2,qajm1,qaj,qajp1,qajp2,eps):
    
    pjm2toja  =(qaj - 2*qajm1 + qajm2)/(2*dx**2)
    pjm2tojb  =(3*qaj - 4*qajm1 + qajm2)/(2*dx)
    pjm2tojc  =23*qaj/24 + qajm1/12 - qajm2/24
    
    pjm1tojp1a  = (-2*qaj + qajm1 + qajp1)/(2*dx**2)
    pjm1tojp1b  = (-qajm1 + qajp1)/(2*dx)
    pjm1tojp1c  = 13*qaj/12 - qajm1/24 - qajp1/24

    
    pjtojp2a  =(qaj - 2*qajp1 + qajp2)/(2*dx**2)
    pjtojp2b  =(-3*qaj + 4*qajp1 - qajp2)/(2*dx)
    pjtojp2c  =23*qaj/24 + qajp1/12 - qajp2/24
    
    Bjm2toj = 10*qaj**2/3 - 31*qaj*qajm1/3 + 11*qaj*qajm2/3 + 25*qajm1**2/3 - 19*qajm1*qajm2/3 + 4*qajm2**2/3
    Bjm1tojp1 =qaj**2/3 - qaj*qajm1/3 - qaj*qajp1/3 + qajm1**2/3 - qajm1*qajp1/3 + qajp1**2/3 + (-2*qaj + qajm1 + qajp1)**2
    Bjtojp2 = 10*qaj**2/3 - 31*qaj*qajp1/3 + 11*qaj*qajp2/3 + 25*qajp1**2/3 - 19*qajp1*qajp2/3 + 4*qajp2**2/3

   
    iw1 = ((3.0/10.0) / (eps +Bjm2toj )**2)
    iw2 = ((3.0/5.0) / (eps +Bjm1tojp1 )**2)
    iw3 = ((1.0/10.0)  / (eps +Bjtojp2 )**2)
    
    w1 = iw1 / (iw1 + iw2 + iw3)
    w2 = iw2 / (iw1 + iw2 + iw3)
    w3 = iw3 / (iw1 + iw2 + iw3)

    qa = w1*pjm2toja + w2*pjm1tojp1a + w3*pjtojp2a
    qb = w1*pjm2tojb + w2*pjm1tojp1b + w3*pjtojp2b
    qc = w1*pjm2tojc + w2*pjm1tojp1c + w3*pjtojp2c

    qjmh = qa*(-dx/2)**2 + qb*(-dx/2) + qc
    
    return qjmh


def Reconjph(qajm2,qajm1,qaj,qajp1,qajp2,eps):
    
    pjm2toja  =(qaj - 2*qajm1 + qajm2)/(2*dx**2)
    pjm2tojb  =(3*qaj - 4*qajm1 + qajm2)/(2*dx)
    pjm2tojc  =23*qaj/24 + qajm1/12 - qajm2/24
    
    pjm1tojp1a  = (-2*qaj + qajm1 + qajp1)/(2*dx**2)
    pjm1tojp1b  = (-qajm1 + qajp1)/(2*dx)
    pjm1tojp1c  = 13*qaj/12 - qajm1/24 - qajp1/24

    pjtojp2a  =(qaj - 2*qajp1 + qajp2)/(2*dx**2)
    pjtojp2b  =(-3*qaj + 4*qajp1 - qajp2)/(2*dx)
    pjtojp2c  =23*qaj/24 + qajp1/12 - qajp2/24
    
    Bjm2toj = 10*qaj**2/3 - 31*qaj*qajm1/3 + 11*qaj*qajm2/3 + 25*qajm1**2/3 - 19*qajm1*qajm2/3 + 4*qajm2**2/3
    Bjm1tojp1 =qaj**2/3 - qaj*qajm1/3 - qaj*qajp1/3 + qajm1**2/3 - qajm1*qajp1/3 + qajp1**2/3 + (-2*qaj + qajm1 + qajp1)**2
    Bjtojp2 = 10*qaj**2/3 - 31*qaj*qajp1/3 + 11*qaj*qajp2/3 + 25*qajp1**2/3 - 19*qajp1*qajp2/3 + 4*qajp2**2/3

   
    iw1 = ((1.0/10.0) / (eps +Bjm2toj )**2)
    iw2 = ((3.0/5.0) / (eps +Bjm1tojp1 )**2)
    iw3 = ((3.0/10.0)  / (eps +Bjtojp2 )**2)
    
    w1 = iw1 / (iw1 + iw2 + iw3)
    w2 = iw2 / (iw1 + iw2 + iw3)
    w3 = iw3 / (iw1 + iw2 + iw3)

    qa = w1*pjm2toja + w2*pjm1tojp1a + w3*pjtojp2a
    qb = w1*pjm2tojb + w2*pjm1tojp1b + w3*pjtojp2b
    qc = w1*pjm2tojc + w2*pjm1tojp1c + w3*pjtojp2c

    qjph = qa*(dx/2)**2 + qb*(dx/2) + qc
    
    return qjph

def FVMSWWE(ha,Ga,nGcells,g,dx,dt):
    
    eps = 10**(-10) 
    n = len(ha)
    hap = zeros(n)
    Gap = zeros(n)
    
    j = nGcells-1
    
    hir = Reconjph(ha[j-2],ha[j-1],ha[j],ha[j+1],ha[j+2],eps)
    hip1l = Reconjmh(ha[j-1],ha[j],ha[j+1],ha[j+2],ha[j+3],eps)
    
    Gir = Reconjph(Ga[j-2],Ga[j-1],Ga[j],Ga[j+1],Ga[j+2],eps)
    Gip1l = Reconjmh(Ga[j-1],Ga[j],Ga[j+1],Ga[j+2],Ga[j+3],eps)
    
    
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
        
        hir = Reconjph(ha[j-2],ha[j-1],ha[j],ha[j+1],ha[j+2],eps)
        hip1l = Reconjmh(ha[j-1],ha[j],ha[j+1],ha[j+2],ha[j+3],eps)
        
        Gir = Reconjph(Ga[j-2],Ga[j-1],Ga[j],Ga[j+1],Ga[j+2],eps)
        Gip1l = Reconjmh(Ga[j-1],Ga[j],Ga[j+1],Ga[j+2],Ga[j+3],eps)
        
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

