
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from numpy.linalg import solve,norm
from matplotlib.pyplot import plot,loglog



def initialhDB(x,dx,h0,h1):
    n = len(x)
    hp = []
    Gp = []
    xp = []

    for i in range(n):
        xjmh = x[i] - dx/2
        xj = x[i]
        xjph = x[i] + dx/2
        
        if x[i] < 0:
            hjmh = h0
            hj = h0
            hjph = h0
        else:
            hjmh = h1
            hj = h1
            hjph = h1

        hcell = [hjmh,hj,hjph]
        xhcell = [xjmh,xj,xjph]
        
        hp = hp + hcell
        xp = xp + xhcell
        Gp = Gp + [0,0,0,0]
        
    return hp,xp,Gp





def GetCellAverage(n,q,dx):
    
    qn = zeros(n)
    #  quadratic passes through qjmh,qj,qjph
    for i in range(n):
        qjmh = q[3*i]
        qj = q[3*i + 1]
        qjph = q[3*i + 2]
        pa = 2*(qjph - 2*qj + qjmh)/ (dx**2)
        pb = (qjph - qjmh) / dx
        pc = qj

        qn[i] = 1.0/dx*(2*pa/3*(dx/2.0)**3 + 2*pc*(dx/2.0) )
    return qn

def CellAverageToCubic(qa,nGcells,dx):
    n = len(qa)
    qcubic = zeros(3*n)
    
    eps = 10.0**(-10)
    # j= 0
    j = 0
    qcubic[3*j: 3*(j+nGcells)] = qa[j]*ones(3*(nGcells))
    
    for j in range(nGcells,n-nGcells):
        
        pjm2toja  = (-2*qa[j-1] + qa[j-2] + qa[j])/(2*dx**2)
        pjm2tojb  =(-4*qa[j-1] + qa[j-2] + 3*qa[j])/(2*dx)
        pjm2tojc  = qa[j-1]/12 - qa[j-2]/24 + 23*qa[j]/24
        
        pjm1tojp1a  = (qa[j+1] + qa[j-1] - 2*qa[j])/(2*dx**2)
        pjm1tojp1b  =(qa[j+1] - qa[j-1])/(2*dx)
        pjm1tojp1c  =-qa[j+1]/24 - qa[j-1]/24 + 13*qa[j]/12
    
        pjtojp2a  = (-2*qa[j+1] + qa[j+2] + qa[j])/(2*dx**2)
        pjtojp2b  = (4*qa[j+1] - qa[j+2] - 3*qa[j])/(2*dx)
        pjtojp2c  = qa[j+1]/12 - qa[j+2]/24 + 23*qa[j]/24
        
        Bjm2toj = 25*qa[j-1]**2/3 - 19*qa[j-1]*qa[j-2]/3 - 31*qa[j-1]*qa[j]/3 + 4*qa[j-2]**2/3 + 11*qa[j-2]*qa[j]/3 + 10*qa[j]**2/3
        Bjm1tojp1 =qa[j+1]**2/3 - qa[j+1]*qa[j-1]/3 - qa[j+1]*qa[j]/3 + qa[j-1]**2/3 - qa[j-1]*qa[j]/3 + qa[j]**2/3 + (qa[j+1] + qa[j-1] - 2*qa[j])**2
        Bjtojp2 = 25*qa[j+1]**2/3 - 19*qa[j+1]*qa[j+2]/3 - 31*qa[j+1]*qa[j]/3 + 4*qa[j+2]**2/3 + 11*qa[j+2]*qa[j]/3 + 10*qa[j]**2/3

       
        iw1 = ((1.0/10.0) / (eps +Bjm2toj )**2)
        iw2 = ((3.0/5.0) / (eps +Bjm1tojp1 )**2)
        iw3 = ((3.0/10.0)  / (eps +Bjtojp2 )**2)

        iw1m = ((3.0/10.0) / (eps +Bjm2toj )**2)
        iw2m = ((3.0/5.0) / (eps +Bjm1tojp1 )**2)
        iw3m = ((1.0/10.0)  / (eps +Bjtojp2 )**2)
        
        w1 = iw1 / (iw1 + iw2 + iw3)
        w2 = iw2 / (iw1 + iw2 + iw3)
        w3 = iw3 / (iw1 + iw2 + iw3)
        
        w1m = iw1m / (iw1m + iw2m + iw3m)
        w2m = iw2m / (iw1m + iw2m + iw3m)
        w3m = iw3m / (iw1m + iw2m + iw3m)
            
        #print (w1,w2,w3,w4)
        
        qac = w1m*pjm2toja + w2m*pjm1tojp1a + w3m*pjtojp2a
        qbc = w1m*pjm2tojb + w2m*pjm1tojp1b + w3m*pjtojp2b
        qcc = w1m*pjm2tojc + w2m*pjm1tojp1c + w3m*pjtojp2c
        
        qcubic[3*j] = qac*(-dx/2)**2 + qbc*(-dx/2) + qcc
        
        
        qac = w1*pjm2toja + w2*pjm1tojp1a + w3*pjtojp2a
        qbc = w1*pjm2tojb + w2*pjm1tojp1b + w3*pjtojp2b
        qcc = w1*pjm2tojc + w2*pjm1tojp1c + w3*pjtojp2c
        
        qcubic[3*j + 1] = qcc
        qcubic[3*j + 2] = qac*(dx/2)**2 + qbc*(dx/2) + qcc


    # j= n-3
    j = n-nGcells
    qcubic[3*j: 3*(j+nGcells)] = qa[j]*ones(3*(nGcells)) 
    return qcubic

def FVMSWWE(ha,Ga,gcub,hcub,nGcells,g,dx,dt):
    
    n = len(ha)
    hap = zeros(n)
    Gap = zeros(n)
    
    j = nGcells-1
    
    hir = hcub[3*j + 2]
    hip1l = hcub[3*j + 3]
    
    Gir = gcub[3*j + 2]
    Gip1l = gcub[3*j + 3]
    
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
        
        hir = hcub[3*j + 2]
        hip1l = hcub[3*j + 3]
        
        Gir = gcub[3*j + 2]
        Gip1l = gcub[3*j + 3]
        
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
    
    hcubN = CellAverageToCubic(hA,nGcells,dx)
    GcubN = CellAverageToCubic(GA,nGcells,dx)
    print(max(hA),max( hcubN))
    hap,Gap = FVMSWWE(hA,GA,GcubN,hcubN,nGcells,g,dx,dt)
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

sx = -100.0
ex = 100.0
dx = ((ex - sx))/(ncurr-1)
#x = arange(sx - nGcells*dx,ex+(nGcells + 0.1)*dx,dx)
x = arange(sx ,ex+(0.1)*dx,dx)
nx = len(x)

hcub,xcub,Gcub = initialhDB(x,dx,2,1)

hA = GetCellAverage(nx,hcub,dx)
GA = GetCellAverage(nx,Gcub,dx)

dt = 0.5*dx/sqrt(2*g)

hAN,GAN = Solver(hA,GA,0,10,g,dx,dt,nGcells)

