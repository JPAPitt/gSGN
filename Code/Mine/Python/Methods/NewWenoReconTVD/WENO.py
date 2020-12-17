
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from numpy.linalg import solve,norm
from matplotlib.pyplot import plot,loglog




def initialhuPeak(x,dx,a1):
    n = len(x)
    hp = []
    up = []
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
        

        hcell = [hjmh,hjms,hjps,hjph]
        xhcell = [xjmh,xjms,xjps,xjph]
        hp = hp + hcell
        xp = xp + xhcell
        

        ujmh = c*(hjmh - 1)/hjmh 
        ujms = c*(hjms - 1)/hjms 
        ujps = c*(hjps - 1)/hjps 
        ujph = c*(hjph - 1)/hjph 
        
        ucell = [ujmh,ujms,ujps,ujph]
        up = up + ucell
        
        
    return hp,up,xp


def initialhDB(x,dx,h0,h1):
    n = len(x)
    hp = []
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
        
        
    return hp,xp


def initialhSmooth(x,dx,h0,h1):
    n = len(x)
    hp = []
    xp = []

    for i in range(n):
        xjmh = x[i] - dx/2
        xjms = x[i] - dx/6
        xj = x[i]
        xjps = x[i] + dx/6
        xjph = x[i] + dx/2
        
        hjmh = h0 + h1*exp( - (xjmh)**2)
        hjms = h0 + h1*exp( - (xjms)**2)
        hjps = h0 + h1*exp( - (xjps)**2)
        hjph = h0 + h1*exp( - (xjph)**2)


        hcell = [hjmh,hjms,hjps,hjph]
        xhcell = [xjmh,xjms,xjps,xjph]
        
        hp = hp + hcell
        xp = xp + xhcell
        
        
    return hp,xp


def initialhDBS(x,dx,h0):
    n = len(x)
    hp = []
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
            hjmh = h0 - xjmh
            hjms = h0 - xjms
            hjps = h0 - xjps
            hjph = h0 - xjph
        

        hcell = [hjmh,hjms,hjps,hjph]
        xhcell = [xjmh,xjms,xjps,xjph]
        
        hp = hp + hcell
        xp = xp + xhcell
        
        
    return hp,xp


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
    
    w1s = zeros(n)
    w2s = zeros(n)
    w3s = zeros(n)
    w4s = zeros(n)
    
    eps = 10.0**(-12)
    # j= 0
    j = 0
    qcubic[4*j: 4*(j+3)] = qA[j]*ones(12)
    
    for j in range(3,n-3):
        
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
                
        Bjm3toj = 11003*qA[j-1]**2/240 - 8623*qA[j-1]*qA[j-2]/120 + 2321*qA[j-1]*qA[j-3]/120 - 1567*qA[j-1]*qA[j]/40 + 7043*qA[j-2]**2/240 - 647*qA[j-2]*qA[j-3]/40 + 3521*qA[j-2]*qA[j]/120 + 547*qA[j-3]**2/240 - 309*qA[j-3]*qA[j]/40 + 2107*qA[j]**2/240
        Bjm2tojp1 = 547*qA[j+1]**2/240 + 961*qA[j+1]*qA[j-1]/120 - 247*qA[j+1]*qA[j-2]/120 - 1261*qA[j+1]*qA[j]/120 + 2843*qA[j-1]**2/240 - 821*qA[j-1]*qA[j-2]/120 - 2983*qA[j-1]*qA[j]/120 + 89*qA[j-2]**2/80 + 267*qA[j-2]*qA[j]/40 + 3443*qA[j]**2/240
        Bjm1tojp2 = 2843*qA[j+1]**2/240 - 821*qA[j+1]*qA[j+2]/120 + 961*qA[j+1]*qA[j-1]/120 - 2983*qA[j+1]*qA[j]/120 + 89*qA[j+2]**2/80 - 247*qA[j+2]*qA[j-1]/120 + 267*qA[j+2]*qA[j]/40 + 547*qA[j-1]**2/240 - 1261*qA[j-1]*qA[j]/120 + 3443*qA[j]**2/240
        Bjtojp3 = 11003*qA[j+1]**2/240 - 8623*qA[j+1]*qA[j+2]/120 + 2321*qA[j+1]*qA[j+3]/120 - 1567*qA[j+1]*qA[j]/40 + 7043*qA[j+2]**2/240 - 647*qA[j+2]*qA[j+3]/40 + 3521*qA[j+2]*qA[j]/120 + 547*qA[j+3]**2/240 - 309*qA[j+3]*qA[j]/40 + 2107*qA[j]**2/240        
 

        iw1 = ((1.0/35.0) / (eps +Bjm3toj )**2)
        iw2 = ((12.0/35.0) / (eps +Bjm2tojp1 )**2)
        iw3 = ((18.0/35.0)  / (eps +Bjm1tojp2 )**2)
        iw4 = ((4.0/35.0)  / (eps +Bjtojp3 )**2)

        
        w1 = iw1 / (iw1 + iw2 + iw3 + iw4 )
        w2 = iw2 / (iw1 + iw2 + iw3 + iw4 )
        w3 = iw3 / (iw1 + iw2 + iw3 + iw4 )
        w4 = iw4 / (iw1 + iw2 + iw3 + iw4 )

        w1s[j] = w1
        w2s[j] = w2
        w3s[j] = w3
        w4s[j] = w4
        
        qa = w1*pjm3toja + w2*pjm2tojp1a + w3*pjm1tojp2a + w4*pjtojp3a
        qb = w1*pjm3tojb + w2*pjm2tojp1b + w3*pjm1tojp2b + w4*pjtojp3b
        qc = w1*pjm3tojc + w2*pjm2tojp1c + w3*pjm1tojp2c + w4*pjtojp3c
        qd = w1*pjm3tojd + w2*pjm2tojp1d + w3*pjm1tojp2d + w4*pjtojp3d
        
        qcubic[4*j] = qa*(-dx/2)**3 + qb*(-dx/2)**2 + qc*(-dx/2) + qd
        qcubic[4*j + 1] = qa*(-dx/6)**3 + qb*(-dx/6)**2 + qc*(-dx/6) + qd
        qcubic[4*j + 2] = qa*(dx/6)**3 + qb*(dx/6)**2 + qc*(dx/6) + qd
        qcubic[4*j + 3] = qa*(dx/2)**3 + qb*(dx/2)**2 + qc*(dx/2) + qd


    # j= n-3
    j = n-3
    qcubic[4*j: 4*(j+3)] = qA[j]*ones(12)
        
    return qcubic,w1s,w2s,w3s,w4s


expn = 10
lown = 10

dxs = []
L2s = []

for expi in range(expn):
    ncurr = lown*(2**expi)
    
    sx = -10.0
    ex = 10.0
    dx = ((ex - sx))/ncurr
    xn = arange(sx,ex+ dx,dx)
    nxn = len(xn)
    a1 = 0.5
    
    hcub,xcub = initialhDBS(xn,dx,2)
    #hcub,xcub = initialhDB(xn,dx,2,1)
    #hcub,xcub = initialhSmooth(xn,dx,1,0.1)
    
    hA = GetCellAverage(nxn,hcub,dx)
    
    hcubN,w1s,w2s,w3s,w4s = CellAverageToCubic(hA,xn,dx)

    L2 = norm(hcub - hcubN ,ord=2)/norm(hcub ,ord=2)

    
    dxs.append(dx)
    L2s.append(L2)
    
# loglog(dxs,L2s,'.')
# loglog(dxs,0.1*array(dxs)**2,'-')
# loglog(dxs,0.1*array(dxs)**4,'-')

plot(xcub,hcub,'-')
plot(xcub,hcubN,'.')
