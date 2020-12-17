
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





def Recon(qa,j,eps):
    
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

    biw1 = (1.0 / (eps +Bjm2toj )**2)
    biw2 = (1.0 / (eps +Bjm1tojp1 )**2)
    biw3 = (1.0 / (eps +Bjtojp2  )**2)
        
    bw1 = biw1 / (biw1 + biw2 + biw3)
    bw2 = biw2 / (biw1 + biw2 + biw3)
    bw3 = biw3 / (biw1 + biw2 + biw3)
    
    opt1 = (1.0/10.0) 
    opt2 = (3.0/5.0)
    opt3 = (1.0/10.0) 
    

    p2a = opt1*bw1*pjm2toja + opt2*bw2*pjm1tojp1a + opt3*bw3*pjtojp2a
    p2b = opt1*bw1*pjm2tojb + opt2*bw2*pjm1tojp1b + opt3*bw3*pjtojp2b
    p2c = opt1*bw1*pjm2tojc + opt2*bw2*pjm1tojp1c + opt3*bw3*pjtojp2c
    
    p2am = opt3*bw1*pjm2toja + opt2*bw2*pjm1tojp1a + opt1*bw3*pjtojp2a
    p2bm = opt3*bw1*pjm2tojb + opt2*bw2*pjm1tojp1b + opt1*bw3*pjtojp2b
    p2cm = opt3*bw1*pjm2tojc + opt2*bw2*pjm1tojp1c + opt1*bw3*pjtojp2c  
    
    p2jmh = p2am*(-dx/2.0)**2 + p2bm*(-dx/2.0) + p2cm
    p2jph = p2a*(dx/2.0)**2 + p2b*(dx/2.0) + p2c
    
    TVDIndLeft = ((qa[j-1] >= qa[j]) and ( qa[j-1] >= p2jmh and p2jmh >= qa[j] )) \
        or ((qa[j-1] <= qa[j]) and ( qa[j-1] <= p2jmh and p2jmh <= qa[j] ))
    TVDIndRight = ((qa[j] >= qa[j+1]) and ( qa[j] >= p2jph and p2jph >= qa[j+1] )) \
        or ((qa[j] <= qa[j+1]) and ( qa[j] <= p2jph and p2jph <= qa[j+1] ))

        
    if (TVDIndLeft and TVDIndRight):
        qmjmh= p2jmh 
        qmjph = p2jph
    else:
        pjm1toja  =(qa[j] - qa[j-1])/(dx)
        pjm1tojb  =qa[j]
        
        pjtojp1a  =(qa[j+1] - qa[j])/(dx)
        pjtojp1b  =qa[j]
        
        
        Bjm1toj = (qa[j] - qa[j-1])**2
        Bjtojp1 = (qa[j+1] - qa[j])**2
        
        opt1 = 1.0/3.0
        op2 = 2.0/3.0
        
        biw1 = (1.0 / (eps + Bjm1toj )**2)
        biw2 = (1.0 / (eps + Bjtojp1 )**2)
        
        bw1 = biw1 / (biw1 + biw2)
        bw2 = biw2 / (biw1 + biw2)
    
        p1a = 0
        p1b = opt1*bw1*pjm1toja + opt2*bw2*pjtojp1a
        p1c = opt1*bw1*pjm1tojb + opt2*bw2*pjtojp1b

        p1am = 0
        p1bm = opt2*bw1*pjm1toja + opt1*bw2*pjtojp1a
        p1cm = opt2*bw1*pjm1tojb + opt1*bw2*pjtojp1b
        
        
        p1jmh = p1bm*(-dx/2.0) + p1cm 
        p1jph = p1b*(dx/2.0) + p1c 
        
        TVDIndLeft = ((qa[j-1] >= qa[j]) and ( qa[j-1] >= p1jmh and p1jmh >= qa[j] )) \
        or ((qa[j-1] <= qa[j]) and ( qa[j-1] <= p1jmh and p1jmh <= qa[j] )) 
        TVDIndRight = ((qa[j] >= qa[j+1]) and ( qa[j] >= p1jph and p1jph >= qa[j+1] )) \
        or ((qa[j] <= qa[j+1]) and ( qa[j] <= p1jph and p1jph <= qa[j+1] ))

        if (TVDIndLeft and TVDIndRight):
            qmjmh= p1jmh 
            qmjph = p1jph
            
        
        else:
    
            p0a = 0
            p0b = 0
            p0c = qa[j]
            
            p0jmh = p0c 
            p0jph = p0c 
            
            qmjmh= p0jmh 
            qmjph = p0jph

    
    return qmjmh,qmjph




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

