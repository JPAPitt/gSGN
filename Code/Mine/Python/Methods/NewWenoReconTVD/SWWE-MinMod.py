
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





def Reconmh(qa,j,eps):
    
    pjm2toja  =(qa[j] - 2*qa[j-1] + qa[j-2])/(2*dx**2)
    pjm2tojb  =(3*qa[j] - 4*qa[j-1] + qa[j-2])/(2*dx)
    pjm2tojc  =23*qa[j]/24 + qa[j-1]/12 - qa[j-2]/24
    
    pjm1tojp1a  = (-2*qa[j] + qa[j-1] + qa[j+1])/(2*dx**2)
    pjm1tojp1b  = (-qa[j-1] + qa[j+1])/(2*dx)
    pjm1tojp1c  = 13*qa[j]/12 - qa[j-1]/24 - qa[j+1]/24

    pjtojp2a  =(qa[j] - 2*qa[j+1] + qa[j+2])/(2*dx**2)
    pjtojp2b  =(-3*qa[j] + 4*qa[j+1] - qa[j+2])/(2*dx)
    pjtojp2c  =23*qa[j]/24 + qa[j+1]/12 - qa[j+2]/24
    
    Bjm2toj = 10*qa[j]**2/3 - 31*qa[j]*qa[j-1]/3 + 11*qa[j]*qa[j-2]/3 + 25*qa[j-1]**2/3 - 19*qa[j-1]*qa[j-2]/3 + 4*qa[j-2]**2/3
    Bjm1tojp1 =qa[j]**2/3 - qa[j]*qa[j-1]/3 - qa[j]*qa[j+1]/3 + qa[j-1]**2/3 - qa[j-1]*qa[j+1]/3 + qa[j+1]**2/3 + (-2*qa[j] + qa[j-1] + qa[j+1])**2
    Bjtojp2 = 10*qa[j]**2/3 - 31*qa[j]*qa[j+1]/3 + 11*qa[j]*qa[j+2]/3 + 25*qa[j+1]**2/3 - 19*qa[j+1]*qa[j+2]/3 + 4*qa[j+2]**2/3

    iw1 = ((3.0/10.0) / (eps +Bjm2toj )**2)
    iw2 = ((3.0/5.0) / (eps +Bjm1tojp1 )**2)
    iw3 = ((1.0/10.0)  / (eps +Bjtojp2 )**2)
    
    w1 = iw1 / (iw1 + iw2 + iw3)
    w2 = iw2 / (iw1 + iw2 + iw3)
    w3 = iw3 / (iw1 + iw2 + iw3)

    p2a = w1*pjm2toja + w2*pjm1tojp1a + w3*pjtojp2a
    p2b = w1*pjm2tojb + w2*pjm1tojp1b + w3*pjtojp2b
    p2c = w1*pjm2tojc + w2*pjm1tojp1c + w3*pjtojp2c   
    
    p2jmh = p2a*(-dx/2.0)**2 + p2b*(-dx/2.0) + p2c 
    
    TVDIndLeft = ((qa[j-1] >= qa[j]) and ( qa[j-1] >= p2jmh and p2jmh >= qa[j] )) \
        or ((qa[j-1] <= qa[j]) and ( qa[j-1] <= p2jmh and p2jmh <= qa[j] ))
        
    if (TVDIndLeft):
        qmjmh= p2jmh 
    else:
        pjm1toja  =(qa[j] - qa[j-1])/(dx)
        pjm1tojb  =qa[j]
        
        pjtojp1a  =(qa[j+1] - qa[j])/(dx)
        pjtojp1b  =qa[j]
        
        
        Bjm1toj = (qa[j] - qa[j-1])**2
        Bjtojp1 = (qa[j+1] - qa[j])**2
        
        iw1 = ((2.0/3.0) / (eps + Bjm1toj )**2)
        iw2 = ((1.0/3.0) / (eps + Bjtojp1 )**2)
        
        w1 = iw1 / (iw1 + iw2)
        w2 = iw2 / (iw1 + iw2)
    
        p1a = 0
        p1b = w1*pjm1toja + w2*pjtojp1a
        p1c = w1*pjm1tojb + w2*pjtojp1b
        
        p1jmh = p1b*(-dx/2.0) + p1c 
        p1jph = p1b*(dx/2.0) + p1c 
        
        TVDIndLeft = ((qa[j-1] >= qa[j]) and ( qa[j-1] >= p1jmh and p1jmh >= qa[j] )) \
        or ((qa[j-1] <= qa[j]) and ( qa[j-1] <= p1jmh and p1jmh <= qa[j] )) 

        if (TVDIndLeft):
            qmjmh= p1jmh            
        else:
    
            p0a = 0
            p0b = 0
            p0c = qa[j]
            
            p0jmh = p0c 
            p0jph = p0c 
            
            qmjmh= p0jmh 

    
    return qmjmh

def minmod(a,b,c):
    if (a > 0) and (b > 0) and (c>0):
        return min(a,b,c)
    elif (a < 0) and (b < 0) and (c<0):
        return max(a,b,c)
    else:
        return 0
    

def Recon(qa,j,eps):
    theta = 1.2
    bd = qa[j] - qa[j-1]
    md = (qa[j+1] - qa[j-1])/2.0
    fd = qa[j+1] - qa[j]
    
    s = minmod(theta*bd,md, theta*fd)
    
    qmh = qa[j] - s/2
    qph = qa[j] + s/2
    
    return qmh,qph



def FVMSWWE(ha,Ga,nGcells,g,dx,dt):
    
    eps = 10**(-10) 
    n = len(ha)
    hap = zeros(n)
    Gap = zeros(n)
    
    j = nGcells-1
    
    hil,hir = Recon(ha,j,eps)
    hip1l,hip1r = Recon(ha,j+1,eps)
    
    Gil,Gir = Recon(Ga,j,eps)
    Gip1l,Gip1r = Recon(Ga,j+1,eps)
    
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
        
        hil,hir = Recon(ha,j,eps)
        hip1l,hip1r = Recon(ha,j+1,eps)
        
        Gil,Gir = Recon(Ga,j,eps)
        Gip1l,Gip1r = Recon(Ga,j+1,eps)
        
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

expi = 6
ncurr = lown*(2**expi)

sx = -50.0
ex = 50.0
dx = ((ex - sx))/(ncurr-1)
#x = arange(sx - nGcells*dx,ex+(nGcells + 0.1)*dx,dx)
x = arange(sx ,ex+(0.1)*dx,dx)
nx = len(x)

hA,GA = initialhDB(x,dx,2,1)


dt = 0.5*dx/sqrt(2*g)

hAN,GAN = Solver(hA,GA,0,10,g,dx,dt,nGcells)

