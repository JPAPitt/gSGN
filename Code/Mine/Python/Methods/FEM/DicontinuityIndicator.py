
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from numpy.linalg import solve,norm
from matplotlib.pyplot import plot,loglog


def initialhu(x,dx,a1):
    n = len(x)
    hp = []
    up = []
    dup = []
    xp = []

    c = sqrt(1 + a1)
        
    for i in range(n):
        xjmh = x[i] - dx/2
        xjms = x[i] - dx/6
        xj = x[i]
        xjps = x[i] + dx/6
        xjph = x[i] + dx/2
        
        #h
        hjmh = 1 + a1*exp(-sqrt(3)*abs(xjmh))
        hjms = 1 + a1*exp(-sqrt(3)*abs(xjms))
        hjps = 1 + a1*exp(-sqrt(3)*abs(xjps))
        hjph = 1 + a1*exp(-sqrt(3)*abs(xjph))
        
        nxmh = -sqrt(3)*a1*sign(xjmh)*exp(-sqrt(3)*abs(xjmh))
        nxms = -sqrt(3)*a1*sign(xjms)*exp(-sqrt(3)*abs(xjms))
        nxps = -sqrt(3)*a1*sign(xjps)*exp(-sqrt(3)*abs(xjps))
        nxph = -sqrt(3)*a1*sign(xjph)*exp(-sqrt(3)*abs(xjph))
        

        hcell = [hjmh,hjms,hjps,hjph]
        xhcell = [xjmh,xjms,xjps,xjph]
        hp = hp + hcell
        xp = xp + xhcell
        

        ujmh = c*(hjmh - 1)/hjmh 
        ujms = c*(hjms - 1)/hjms 
        ujps = c*(hjps - 1)/hjps 
        ujph = c*(hjph - 1)/hjph 
        
        uxjmh = c*(nxmh)/hjmh**2
        uxjms = c*(nxms)/hjms**2
        uxjps = c*(nxps)/hjps**2
        uxjph = c*(nxph)/hjph**2
        
        ucell = [ujmh,ujms,ujps,ujph]
        uxcell = [uxjmh,uxjms,uxjps,uxjph]
        up = up + ucell
        dup = dup + uxcell
        
    
        
    return hp,up,dup,xp

def GetCellMidpoint(n,q,dx):
    
    qn = zeros(n)
    # cubic passes through y1, y2, y3, y4
    # pa = -9 Gjmh + 27 Gms - 27 Gps + 9 Gph / 2*dx3
    # pb = 9 Gjmh - 9 Gms - 9 Gps + 9 Gph / 4*dx2
    # pc = Gjmh - 27 Gms + 27 Gps - Gph / 8*dx
    # pd = -Gjmh + 9 Gms + 9 Gps - Gph / 16
    for i in range(n):
        qmh = q[4*i]
        qms = q[4*i + 1]
        qps = q[4*i + 2]
        qph = q[4*i + 3]
        qn[i] = (-qmh + 9*qms + 9*qps - qph) / 16.0
    
    return qn

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

def CellAverageToCubic(qA,x,dx):
    n = len(qA)
    qcubic = zeros(4*n)
    dqcubic = zeros(4*n)
    
    eps = 10.0**(-10)
    # j= 0
    j = 0
    qcubic[4*j: 4*(j+4)] = qA[j]*ones(16)
    dqcubic[4*j: 4*(j+4)] = zeros(16)
    for j in range(4,n-4):
        
        pjm3toja  = (qA[j-1] - 3*qA[j-2] + 3*qA[j-3] - qA[j-4])/(6*dx**3)
        pjm3tojb  = (3*qA[j-1] - 8*qA[j-2] + 7*qA[j-3] - 2*qA[j-4])/(2*dx**2)
        pjm3tojc  = (103*qA[j-1] - 225*qA[j-2] + 165*qA[j-3] - 43*qA[j-4])/(24*dx)
        pjm3tojd  =  31*qA[j-1]/8 - 17*qA[j-2]/3 + 89*qA[j-3]/24 - 11*qA[j-4]/12

 
        pjtojp3a = (-qA[j+1] + 3*qA[j+2] - 3*qA[j+3] + qA[j+4])/(6*dx**3)
        pjtojp3b = (3*qA[j+1] - 8*qA[j+2] + 7*qA[j+3] - 2*qA[j+4])/(2*dx**2)
        pjtojp3c =  (-103*qA[j+1] + 225*qA[j+2] - 165*qA[j+3] + 43*qA[j+4])/(24*dx)
        pjtojp3d  = 31*qA[j+1]/8 - 17*qA[j+2]/3 + 89*qA[j+3]/24 - 11*qA[j+4]/12

                
        # Bjm3toj = 11003*qA[j-1]**2/240 - 8623*qA[j-1]*qA[j-2]/120 + 2321*qA[j-1]*qA[j-3]/120 - 1567*qA[j-1]*qA[j]/40 + 7043*qA[j-2]**2/240 - 647*qA[j-2]*qA[j-3]/40 + 3521*qA[j-2]*qA[j]/120 + 547*qA[j-3]**2/240 - 309*qA[j-3]*qA[j]/40 + 2107*qA[j]**2/240       
        # Bjm2tojp1 = 547*qA[j+1]**2/240 + 961*qA[j+1]*qA[j-1]/120 - 247*qA[j+1]*qA[j-2]/120 - 1261*qA[j+1]*qA[j]/120 + 2843*qA[j-1]**2/240 - 821*qA[j-1]*qA[j-2]/120 - 2983*qA[j-1]*qA[j]/120 + 89*qA[j-2]**2/80 + 267*qA[j-2]*qA[j]/40 + 3443*qA[j]**2/240      
        # Bjm1tojp2 = 2843*qA[j+1]**2/240 - 821*qA[j+1]*qA[j+2]/120 + 961*qA[j+1]*qA[j-1]/120 - 2983*qA[j+1]*qA[j]/120 + 89*qA[j+2]**2/80 - 247*qA[j+2]*qA[j-1]/120 + 267*qA[j+2]*qA[j]/40 + 547*qA[j-1]**2/240 - 1261*qA[j-1]*qA[j]/120 + 3443*qA[j]**2/240         
        # Bjtojp3 = 11003*qA[j+1]**2/240 - 8623*qA[j+1]*qA[j+2]/120 + 2321*qA[j+1]*qA[j+3]/120 - 1567*qA[j+1]*qA[j]/40 + 7043*qA[j+2]**2/240 - 647*qA[j+2]*qA[j+3]/40 + 3521*qA[j+2]*qA[j]/120 + 547*qA[j+3]**2/240 - 309*qA[j+3]*qA[j]/40 + 2107*qA[j]**2/240
            
        # iw1 = (0.15 / (eps +Bjm3toj )**2)**(1)
        # iw2 = (0.35 / (eps +Bjm2tojp1 )**2)**(1)
        # iw3 = (0.35 / (eps +Bjm1tojp2 )**2)**(1)
        # iw4 = (0.15 / (eps +Bjtojp3 )**2)**(1)
        
        if x[4*j] < 0:
            qa = pjm3toja
            qb = pjm3tojb
            qc = pjm3tojc
            qd = pjm3tojd
            
        else:
            qa = pjtojp3a
            qb = pjtojp3b
            qc = pjtojp3c
            qd = pjtojp3d
        
        qcubic[4*j] = qa*(-dx/2)**3 + qb*(-dx/2)**2 + qc*(-dx/2) + qd
        dqcubic[4*j] = 3*qa*(-dx/2)**2 + 2*qb*(-dx/2) + qc
        
        if x[4*j+1] < 0:
            qa = pjm3toja
            qb = pjm3tojb
            qc = pjm3tojc
            qd = pjm3tojd
            
        else:
            qa = pjtojp3a
            qb = pjtojp3b
            qc = pjtojp3c
            qd = pjtojp3d
        
        qcubic[4*j + 1] = qa*(-dx/6)**3 + qb*(-dx/6)**2 + qc*(-dx/6) + qd
        dqcubic[4*j + 1]  = 3*qa*(-dx/6)**2 + 2*qb*(-dx/6) + qc
        if x[4*j+2] < 0:
            qa = pjm3toja
            qb = pjm3tojb
            qc = pjm3tojc
            qd = pjm3tojd
            
        else:
            qa = pjtojp3a
            qb = pjtojp3b
            qc = pjtojp3c
            qd = pjtojp3d
        
        qcubic[4*j + 2] = qa*(dx/6)**3 + qb*(dx/6)**2 + qc*(dx/6) + qd
        dqcubic[4*j + 2]  = 3*qa*(dx/6)**2 + 2*qb*(dx/6) + qc
        if x[4*j+3] < 0:
            qa = pjm3toja
            qb = pjm3tojb
            qc = pjm3tojc
            qd = pjm3tojd
            
        else:
            qa = pjtojp3a
            qb = pjtojp3b
            qc = pjtojp3c
            qd = pjtojp3d
        
        qcubic[4*j + 3] = qa*(dx/2)**3 + qb*(dx/2)**2 + qc*(dx/2) + qd
        dqcubic[4*j + 3]  = 3*qa*(dx/2)**2 + 2*qb*(dx/2) + qc
        
        # w1 = iw1 / (iw1 + iw2 + iw3 + iw4 )
        # w2 = iw2 / (iw1 + iw2 + iw3 + iw4 )
        # w3 = iw3 / (iw1 + iw2 + iw3 + iw4 )
        # w4 = iw4 / (iw1 + iw2 + iw3 + iw4 )
        # # if j > n/2 - 10 and j < n/2 + 10:
        # #    print(w1,w2,w3,w4)
        # qa = w1*pjm3toja + w2*pjm2tojp1a + w3*pjm1tojp2a + w4*pjtojp3a
        # qb = w1*pjm3tojb + w2*pjm2tojp1b + w3*pjm1tojp2b + w4*pjtojp3b
        # qc = w1*pjm3tojc + w2*pjm2tojp1c + w3*pjm1tojp2c + w4*pjtojp3c
        # qd = w1*pjm3tojd + w2*pjm2tojp1d + w3*pjm1tojp2d + w4*pjtojp3d
        
        # qcubic[4*j] = qa*(-dx/2)**3 + qb*(-dx/2)**2 + qc*(-dx/2) + qd
        # qcubic[4*j + 1] = qa*(-dx/6)**3 + qb*(-dx/6)**2 + qc*(-dx/6) + qd
        # qcubic[4*j + 2] = qa*(dx/6)**3 + qb*(dx/6)**2 + qc*(dx/6) + qd
        # qcubic[4*j + 3] = qa*(dx/2)**3 + qb*(dx/2)**2 + qc*(dx/2) + qd


    # j= n-3
    j = n-4
    qcubic[4*j: 4*(j+4)] = qA[j]*ones(16)
    dqcubic[4*j: 4*(j+4)] = zeros(16)  
    return qcubic,dqcubic


expn = 10
lown = 100

dxs = []
L2s = []
L2sdu = []

for expi in range(expn):
    ncurr = lown*(1.5**expi)
    
    sx = -10.0
    ex = 10.0
    dx = ((ex - sx))/ncurr
    xn = arange(sx ,ex+ dx,dx)
    nxn = len(xn)
    a1 = 0.5
    
    h,u,du,x = initialhu(xn,dx,a1)
    
    # AG,bG = FEMforG(xn,h,u,dx)

    # Gc = solve(AG, bG.dot(u))
    
    # Au,bu = FEMforu(xn, h,Gc,dx)

    # uc = solve(Au, bu.dot(Gc))
    
    uA = GetCellAverage(nxn,u,dx)
    
    ucN,ducN = CellAverageToCubic(uA,x,dx)
    
    L2 = norm(u - ucN ,ord=2)/norm(u ,ord=2)
    L2du = norm(du - ducN ,ord=2)/norm(du ,ord=2)
    
    dxs.append(dx)
    L2s.append(L2)
    L2sdu.append(L2du)
    
    
loglog(dxs,L2s,'.')
loglog(dxs, L2sdu,'.')
loglog(dxs,array(dxs)**2,'-')
loglog(dxs,array(dxs)**3,'-')
loglog(dxs,array(dxs)**4,'-')


# n= 1000
# sx = -10.0
# ex = 10.0
# dx = ((ex - sx))/n
# xn = arange(sx ,ex+ dx,dx)
# nxn = len(xn)
# a1 = 0.5


# h,u,x = initialhu(xn,dx,a1)

# Guh = array(h)*array(u)

# AG,bG = FEMforG(xn,h,u,dx)

# Gc = solve(AG, bG.dot(u))

# # Gm = GetCellMidpoint(nxn,Gc,dx)
# # GA = GetCellAverage(nxn,Gc,dx)


# Au,bu = FEMforu(xn, h,Gc,dx)

# uc = solve(Au, bu.dot(Gc))

# uA = GetCellAverage(nxn,uc,dx)

# ucN = CellAverageToCubic(uA,dx)

