
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





def Recon(qA,j,eps):
    
    eps = 10.0**(-12)
    theta = 0.5
    
    ErrorMeasure0 = ( abs(qA[j] -qA[j-1] ) + abs(qA[j] -qA[j+1] ) + abs(qA[j] -qA[j+2] ))/3.0
    P0jmh = qA[j] 
    P0jms = qA[j] 
    P0jps = qA[j] 
    P0jph = qA[j]  
    
    #second poly - (qA[j+1] - qA[j])/dx*(x - xj) +    qA[j]                    
    P1jmh = -(qA[j+1] - qA[j-1])/4.0 + qA[j]  
    P1jms = -(qA[j+1] - qA[j-1])/12.0 + qA[j] 
    P1jps = (qA[j+1] - qA[j-1])/12.0 + qA[j] 
    P1jph = (qA[j+1] - qA[j-1])/4.0 + qA[j]  
    
    Errjp2 = abs((qA[j] + qA[j+1] - qA[j-1]) - qA[j+2] )
    ErrorMeasure1 = (Errjp2)
    
    MML1 = (P1jmh - qA[j]) / (qA[j-1] - qA[j] + eps)
    MMR1 = (P1jph - qA[j]) / (qA[j+1] - qA[j] + eps )
    
    TVDIndL1 = (MML1 > 0 ) and (MML1 < theta )
    TVDIndR1 = (MMR1 > 0 ) and (MMR1 < theta )
    
    #third poly -   
    p2a  = (-2*qA[j] + qA[j-1] + qA[j+1])/(2*dx**2)
    p2b  = (-qA[j-1] + qA[j+1])/(2*dx)
    p2c  = 13*qA[j]/12 - qA[j-1]/24 - qA[j+1]/24
    
    P2jmh = p2a*(-dx/2)**2 + p2b*(-dx/2) + p2c
    P2jms = p2a*(-dx/6)**2 + p2b*(-dx/6) + p2c
    P2jps = p2a*(dx/6)**2 + p2b*(dx/6) + p2c
    P2jph = p2a*(dx/2)**2 + p2b*(dx/2) + p2c
    
    MML2 = (P2jmh - qA[j]) / (qA[j-1] - qA[j] + eps)
    MMR2 = (P2jph - qA[j]) / (qA[j+1] - qA[j] + eps )
    
    TVDIndL2 = (MML2 > 0 ) and (MML2 < theta )
    TVDIndR2 = (MMR2 > 0 ) and (MMR2 < theta )
        
    Errjp2 = abs((-3*qA[j] + 3*qA[j+1] + qA[j-1]) - qA[j+2] )
    ErrorMeasure2 =Errjp2
    
    #fourth poly  
    p3a = (-3*qA[j+1] + qA[j+2] - qA[j-1] + 3*qA[j])/(6*dx**3)
    p3b = (qA[j+1] + qA[j-1] - 2*qA[j])/(2*dx**2)
    p3c = (27*qA[j+1] - 5*qA[j+2] - 7*qA[j-1] - 15*qA[j])/(24*dx)
    p3d  = -qA[j+1]/24 - qA[j-1]/24 + 13*qA[j]/12
    
    P3jmh= p3a*(-dx/2)**3 + p3b*(-dx/2)**2 + p3c*(-dx/2) + p3d
    P3jms = p3a*(-dx/6)**3 + p3b*(-dx/6)**2 + p3c*(-dx/6) + p3d
    P3jps = p3a*(dx/6)**3 + p3b*(dx/6)**2 + p3c*(dx/6) + p3d
    P3jph = p3a*(dx/2)**3 + p3b*(dx/2)**2 + p3c*(dx/2) + p3d
            
    MML3 = (P3jmh - qA[j]) / (qA[j-1] - qA[j] + eps)
    MMR3 = (P3jph - qA[j]) / (qA[j+1] - qA[j] + eps )
    
    TVDIndL3 = (MML3 > 0 ) and (MML3 < theta )
    TVDIndR3 = (MMR3 > 0 ) and (MMR3 < theta )

    # print(j, ErrorMeasure0,ErrorMeasure1,ErrorMeasure2,TVDIndL1,TVDIndR1,TVDIndL2,TVDIndR2,TVDIndL3,TVDIndR3)
    if (ErrorMeasure1 > ErrorMeasure2) and TVDIndL3 and TVDIndR3 :
       qjmh = P3jmh
       qjph = P3jph
    elif (ErrorMeasure1 > ErrorMeasure2) and TVDIndL2 and TVDIndR2:
       qjmh = P2jmh
       qjph = P2jph
    elif (ErrorMeasure0 > ErrorMeasure1) and TVDIndL1 and TVDIndR1:
       qjmh = P1jmh
       qjph = P1jph
    else:
       qjmh = P0jmh
       qjph = P0jph
        
    return qjmh,qjph




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

hAN,GAN = Solver(hA,GA,0,4,g,dx,dt,nGcells)

