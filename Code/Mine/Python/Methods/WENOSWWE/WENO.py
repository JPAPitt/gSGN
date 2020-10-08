
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
        xjms = x[i] - dx/6
        xj = x[i]
        xjps = x[i] + dx/6
        xjph = x[i] + dx/2
        
        if x[i] < 0:
            hjmh = h0
            hjms = h0
            hjps = h0
            hjph = h0
        else:
            hjmh = h1
            hjms = h1
            hjps = h1
            hjph = h1

        hcell = [hjmh,hjms,hjps,hjph]
        xhcell = [xjmh,xjms,xjps,xjph]
        
        hp = hp + hcell
        xp = xp + xhcell
        Gp = Gp + [0,0,0,0]
        
    return hp,xp,Gp





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

def CellAverageToCubic(qA,nGcells,dx):
    n = len(qA)
    qcubic = zeros(4*n)
    
    eps = 10.0**(-10)
    # j= 0
    j = 0
    qcubic[4*j: 4*(j+nGcells)] = qA[j]*ones(4*(nGcells))
    
    for j in range(nGcells,n-nGcells):
        
        pjm3toja  = (-3*qA[j-1] + 3*qA[j-2] - qA[j-3] + qA[j])/(6*dx**3)
        pjm3tojb  = (-5*qA[j-1] + 4*qA[j-2] - qA[j-3] + 2*qA[j])/(2*dx**2)
        pjm3tojc  = (-69*qA[j-1] + 33*qA[j-2] - 7*qA[j-3] + 43*qA[j])/(24*dx)
        pjm3tojd  =  5*qA[j-1]/24 - qA[j-2]/6 + qA[j-3]/24 + 11*qA[j]/12
    
        pjm2tojp1a  = (qA[j+1] + 3*qA[j-1] - qA[j-2] - 3*qA[j])/(6*dx**3)
        pjm2tojp1b  = (qA[j+1] + qA[j-1] - 2*qA[j])/(2*dx**2)
        pjm2tojp1c  = (7*qA[j+1] - 27*qA[j-1] + 5*qA[j-2] + 15*qA[j])/(24*dx)
        pjm2tojp1d  = -qA[j+1]/24 - qA[j-1]/24 + 13*qA[j]/12
        
        pjm1tojp2a = (-3*qA[j+1] + qA[j+2] - qA[j-1] + 3*qA[j])/(6*dx**3)
        pjm1tojp2b = (qA[j+1] + qA[j-1] - 2*qA[j])/(2*dx**2)
        pjm1tojp2c = (27*qA[j+1] - 5*qA[j+2] - 7*qA[j-1] - 15*qA[j])/(24*dx)
        pjm1tojp2d  = -qA[j+1]/24 - qA[j-1]/24 + 13*qA[j]/12
 
        pjtojp3a = (3*qA[j+1] - 3*qA[j+2] + qA[j+3] - qA[j])/(6*dx**3)
        pjtojp3b = (-5*qA[j+1] + 4*qA[j+2] - qA[j+3] + 2*qA[j])/(2*dx**2)
        pjtojp3c = (69*qA[j+1] - 33*qA[j+2] + 7*qA[j+3] - 43*qA[j])/(24*dx)
        pjtojp3d  = 5*qA[j+1]/24 - qA[j+2]/6 + qA[j+3]/24 + 11*qA[j]/12 
                
        # Bjm3toj = 11003*qA[j-1]**2/240 - 8623*qA[j-1]*qA[j-2]/120 + 2321*qA[j-1]*qA[j-3]/120 - 1567*qA[j-1]*qA[j]/40 + 7043*qA[j-2]**2/240 - 647*qA[j-2]*qA[j-3]/40 + 3521*qA[j-2]*qA[j]/120 + 547*qA[j-3]**2/240 - 309*qA[j-3]*qA[j]/40 + 2107*qA[j]**2/240
        # Bjm2tojp1 = 547*qA[j+1]**2/240 + 961*qA[j+1]*qA[j-1]/120 - 247*qA[j+1]*qA[j-2]/120 - 1261*qA[j+1]*qA[j]/120 + 2843*qA[j-1]**2/240 - 821*qA[j-1]*qA[j-2]/120 - 2983*qA[j-1]*qA[j]/120 + 89*qA[j-2]**2/80 + 267*qA[j-2]*qA[j]/40 + 3443*qA[j]**2/240
        # Bjm1tojp2 = 2843*qA[j+1]**2/240 - 821*qA[j+1]*qA[j+2]/120 + 961*qA[j+1]*qA[j-1]/120 - 2983*qA[j+1]*qA[j]/120 + 89*qA[j+2]**2/80 - 247*qA[j+2]*qA[j-1]/120 + 267*qA[j+2]*qA[j]/40 + 547*qA[j-1]**2/240 - 1261*qA[j-1]*qA[j]/120 + 3443*qA[j]**2/240
        # Bjtojp3 = 11003*qA[j+1]**2/240 - 8623*qA[j+1]*qA[j+2]/120 + 2321*qA[j+1]*qA[j+3]/120 - 1567*qA[j+1]*qA[j]/40 + 7043*qA[j+2]**2/240 - 647*qA[j+2]*qA[j+3]/40 + 3521*qA[j+2]*qA[j]/120 + 547*qA[j+3]**2/240 - 309*qA[j+3]*qA[j]/40 + 2107*qA[j]**2/240        
 

        Bjm3toj = 8843*qA[j-1]**2/240 - 6463*qA[j-1]*qA[j-2]/120 + 1601*qA[j-1]*qA[j-3]/120 - 1327*qA[j-1]*qA[j]/40 + 4883*qA[j-2]**2/240 - 407*qA[j-2]*qA[j-3]/40 + 2801*qA[j-2]*qA[j]/120 + 307*qA[j-3]**2/240 - 229*qA[j-3]*qA[j]/40 + 1867*qA[j]**2/240
        Bjm2tojp1 = 307*qA[j+1]**2/240 + 241*qA[j+1]*qA[j-1]/120 - 7*qA[j+1]*qA[j-2]/120 - 541*qA[j+1]*qA[j]/120 + 683*qA[j-1]**2/240 - 101*qA[j-1]*qA[j-2]/120 - 823*qA[j-1]*qA[j]/120 + 9*qA[j-2]**2/80 + 27*qA[j-2]*qA[j]/40 + 1283*qA[j]**2/240
        Bjm1tojp2 = 683*qA[j+1]**2/240 - 101*qA[j+1]*qA[j+2]/120 + 241*qA[j+1]*qA[j-1]/120 - 823*qA[j+1]*qA[j]/120 + 9*qA[j+2]**2/80 - 7*qA[j+2]*qA[j-1]/120 + 27*qA[j+2]*qA[j]/40 + 307*qA[j-1]**2/240 - 541*qA[j-1]*qA[j]/120 + 1283*qA[j]**2/240
        Bjtojp3 =8843*qA[j+1]**2/240 - 6463*qA[j+1]*qA[j+2]/120 + 1601*qA[j+1]*qA[j+3]/120 - 1327*qA[j+1]*qA[j]/40 + 4883*qA[j+2]**2/240 - 407*qA[j+2]*qA[j+3]/40 + 2801*qA[j+2]*qA[j]/120 + 307*qA[j+3]**2/240 - 229*qA[j+3]*qA[j]/40 + 1867*qA[j]**2/240       

       
        iw1 = ((1.0/35.0) / (eps +Bjm3toj )**2)
        iw2 = ((12.0/35.0) / (eps +Bjm2tojp1 )**2)
        iw3 = ((18.0/35.0)  / (eps +Bjm1tojp2 )**2)
        iw4 = ((4.0/35.0)  / (eps +Bjtojp3 )**2)

        
        w1 = iw1 / (iw1 + iw2 + iw3 + iw4 )
        w2 = iw2 / (iw1 + iw2 + iw3 + iw4 )
        w3 = iw3 / (iw1 + iw2 + iw3 + iw4 )
        w4 = iw4 / (iw1 + iw2 + iw3 + iw4 )
        
        if w1 < eps:
            w1 = 0
        if w2 < eps:
            w2 = 0       
        if w3 < eps:
            w3 = 0
        if w4 < eps:
            w4 = 0
            
        #print (w1,w2,w3,w4)
        
        qa = w1*pjm3toja + w2*pjm2tojp1a + w3*pjm1tojp2a + w4*pjtojp3a
        qb = w1*pjm3tojb + w2*pjm2tojp1b + w3*pjm1tojp2b + w4*pjtojp3b
        qc = w1*pjm3tojc + w2*pjm2tojp1c + w3*pjm1tojp2c + w4*pjtojp3c
        qd = w1*pjm3tojd + w2*pjm2tojp1d + w3*pjm1tojp2d + w4*pjtojp3d
        
        qcubic[4*j] = qa*(-dx/2)**3 + qb*(-dx/2)**2 + qc*(-dx/2) + qd
        qcubic[4*j + 1] = qa*(-dx/6)**3 + qb*(-dx/6)**2 + qc*(-dx/6) + qd
        qcubic[4*j + 2] = qa*(dx/6)**3 + qb*(dx/6)**2 + qc*(dx/6) + qd
        qcubic[4*j + 3] = qa*(dx/2)**3 + qb*(dx/2)**2 + qc*(dx/2) + qd


    # j= n-3
    j = n-nGcells
    qcubic[4*j: 4*(j+nGcells)] = qA[j]*ones(4*(nGcells)) 
    return qcubic

def FVMSWWE(ha,Ga,gcub,hcub,nGcells,g,dx,dt):
    
    n = len(ha)
    hap = zeros(n)
    Gap = zeros(n)
    
    j = nGcells-1
    
    hir = hcub[4*j + 3]
    hip1l = hcub[4*j + 4]
    
    Gir = gcub[4*j + 3]
    Gip1l = gcub[4*j + 4]
    
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
        
        hir = hcub[4*j + 3]
        hip1l = hcub[4*j + 4]
        
        Gir = gcub[4*j + 3]
        Gip1l = gcub[4*j + 4]
        
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

expi = 5
ncurr = lown*(2**expi)

sx = -50.0
ex = 50.0
dx = ((ex - sx))/(ncurr-1)
#x = arange(sx - nGcells*dx,ex+(nGcells + 0.1)*dx,dx)
x = arange(sx ,ex+(0.1)*dx,dx)
nx = len(x)

hcub,xcub,Gcub = initialhDB(x,dx,2,1)

hA = GetCellAverage(nx,hcub,dx)
GA = GetCellAverage(nx,Gcub,dx)

dt = 0.5*dx/sqrt(2*g)

hAN,GAN = Solver(hA,GA,0,10,g,dx,dt,nGcells)

