
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




def initialhSmooth(x,dx,h0,h1):
    n = len(x)
    hp = []
    xp = []
    hA = zeros(n)

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
        
        hAjmh = h0*(xjmh) + 1/2*sqrt(pi)*h1*math.erf(xjmh)
        hAjph = h0*(xjph) + 1/2*sqrt(pi)*h1*math.erf(xjph)
        
        hA[i] = (hAjph - hAjmh)/dx 
        
        
    return hp,xp,hA


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
        
        if xjmh <= 0:
            hjmh = h0
        else:
            hjmh = h0 - xjmh
            
        if xjms <= 0:
            hjms = h0
        else:
            hjms = h0 - xjms
 
        if xjps <= 0:
            hjps = h0
        else:
            hjps = h0 - xjps

        if xjph <= 0:
            hjph = h0
        else:
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
    qsix = zeros(4*n)
    
    w1s = zeros(n)
    w2s = zeros(n)
    w3s = zeros(n)
    w4s = zeros(n)
    
    eps = 10.0**(-12)
    # j= 0
    j = 0
    qcubic[4*j: 4*(j+3)] = qA[j]*ones(12)
    
    qsix[4*j: 4*(j+3)] = qA[j]*ones(12)
    for j in range(3,n-3):

        
        qas =(15*qA[j+1] - 6*qA[j+2] + qA[j+3] + 15*qA[j-1] - 6*qA[j-2] + qA[j-3] - 20*qA[j])/(720*dx**6)
        qbs =(5*qA[j+1] - 4*qA[j+2] + qA[j+3] - 5*qA[j-1] + 4*qA[j-2] - qA[j-3])/(240*dx**5)
        qcs =(-171*qA[j+1] + 54*qA[j+2] - 5*qA[j+3] - 171*qA[j-1] + 54*qA[j-2] - 5*qA[j-3] + 244*qA[j])/(576*dx**4)
        qds =(-83*qA[j+1] + 52*qA[j+2] - 7*qA[j+3] + 83*qA[j-1] - 52*qA[j-2] + 7*qA[j-3])/(288*dx**3)
        qes =(3435*qA[j+1] - 462*qA[j+2] + 37*qA[j+3] + 3435*qA[j-1] - 462*qA[j-2] + 37*qA[j-3] - 6020*qA[j])/(3840*dx**2)
        qfs =(9455*qA[j+1] - 2236*qA[j+2] + 259*qA[j+3] - 9455*qA[j-1] + 2236*qA[j-2] - 259*qA[j-3])/(11520*dx)
        qgs =-7621*qA[j+1]/107520 + 159*qA[j+2]/17920 - 5*qA[j+3]/7168 - 7621*qA[j-1]/107520 + 159*qA[j-2]/17920 - 5*qA[j-3]/7168 + 30251*qA[j]/26880
                
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

        
        # print(j,Bjm3toj,Bjm2tojp1,Bjm1tojp2,Bjtojp3)
        # print(j, iw1,iw2,iw3,iw4)
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
        
        #qsix[4*j] = qas*(-dx/2)**6 + qbs*(-dx/2)**5 + qcs*(-dx/2)**4 + qds*(-dx/2)**3 + qes*(-dx/2)**2 + qfs*(-dx/2) + qgs
        #qsix[4*j + 1] = qas*(-dx/6)**6 + qbs*(-dx/6)**5 + qcs*(-dx/6)**4 + qds*(-dx/6)**3 + qes*(-dx/6)**2 + qfs*(-dx/6) + qgs
        #qsix[4*j + 2] = qas*(dx/6)**6 + qbs*(dx/6)**5 + qcs*(dx/6)**4 + qds*(dx/6)**3 + qes*(dx/6)**2 + qfs*(dx/6) + qgs
        qsix[4*j + 3] = qas*(dx/2)**6 + qbs*(dx/2)**5 + qcs*(dx/2)**4 + qds*(dx/2)**3 + qes*(dx/2)**2 + qfs*(dx/2) + qgs
        # qcubic[4*j + 1] = qa*(-dx/6)**3 + qb*(-dx/6)**2 + qc*(-dx/6) + qd
        # qcubic[4*j + 2] = qa*(dx/6)**3 + qb*(dx/6)**2 + qc*(dx/6) + qd
        # qcubic[4*j + 3] = qa*(dx/2)**3 + qb*(dx/2)**2 + qc*(dx/2) + qd


    # j= n-3
    j = n-3
    qcubic[4*j: 4*(j+3)] = qA[j]*ones(12)
    qsix[4*j: 4*(j+3)] = qA[j]*ones(12)
        
    return qcubic,qsix,w1s,w2s,w3s,w4s


expn = 10
lown = 2

dxs = []
L2s = []
L2phs = []
L2phsixs = []

for expi in range(expn):
    ncurr = lown*(2**expi)
    
    sx = -10.0
    ex = 10.0
    dx = ((ex - sx))/ncurr
    xn = arange(sx,ex+ dx,dx)
    nxn = len(xn)
    a1 = 0.5
    
    #hcub,xcub = initialhDB(xn,dx,2,1)
    #hcub,xcub = initialhDBS(xn,dx,2)
    hcub,xcub,hA = initialhSmooth(xn,dx,1,0.1)
    
    #hA = GetCellAverage(nxn,hcub,dx)
    
    hcubN,hsix,w1s,w2s,w3s,w4s = CellAverageToCubic(hA,xn,dx)

    L2 = norm(hcub- hcubN ,ord=2)/norm(hcub,ord=2)
    L2ph = norm(hcub[3::4]- hcubN[3::4] ,ord=2)/norm(hcub[3::4] ,ord=2)
    L2phsix = norm(hcub[3::4] - hsix[3::4]  ,ord=2)/norm(hcub[3::4]  ,ord=2)

    
    dxs.append(dx)
    L2s.append(L2)
    L2phs.append(L2ph)
    L2phsixs.append(L2phsix)
    
loglog(dxs,L2s,'.')
loglog(dxs,L2phs,'o')
loglog(dxs,L2phsixs,'+')
loglog(dxs,0.8*array(dxs)**2,'-')
loglog(dxs,0.8*array(dxs)**4,'-')
loglog(dxs,0.8*array(dxs)**6,'-')
loglog(dxs,0.8*array(dxs)**7,'-')

