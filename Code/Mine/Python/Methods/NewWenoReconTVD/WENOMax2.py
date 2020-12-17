
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from numpy.linalg import solve,norm
from matplotlib.pyplot import plot,loglog



def initialhDB(x,dx,h0,h1):
    n = len(x)
    hp = []
    xp = []
    hA = zeros(n)

    for i in range(n):
        
        if x[i] < 0:
            hA[i] = h0
        else:
            hA[i] = h1
            
        xjmh = x[i] - dx/2
        xj = x[i]
        xjph = x[i] + dx/2
        
        if xjmh < 0:
            hjmh = h0
        else:
            hjmh = h1

        if xj < 0:
            hj = h0
        else:
            hj = h1

        if xjph < 0:
            hjph = h0
        else:
            hjph = h1
            
        hcell = [hjmh,hj,hjph]
        xhcell = [xjmh,xj,xjph]
        hp = hp + hcell
        xp = xp + xhcell
            
    return hp,xp,hA


def initialPeak(x,dx,a1):
    n = len(x)
    hp = []
    xp = []
    
    hA = zeros(n)
    c = sqrt(1 + a1)
        
    for i in range(n):
        xjmh = x[i] - dx/2
        xj = x[i]
        xjph = x[i] + dx/2
        
        #h
        hjmh = 1 + a1*exp(-sqrt(3)*abs(xjmh))
        hj = 1 + a1*exp(-sqrt(3)*abs(xj))
        hjph = 1 + a1*exp(-sqrt(3)*abs(xjph))
                      
        
        hcell = [hjmh,hj,hjph]
        xhcell = [xjmh,xj,xjph]
        hp = hp + hcell
        xp = xp + xhcell
        
        
        if (0 < x[i] - 0.5*dx):
            hA[i] = 1 + a1 /sqrt(3) *(-exp(-sqrt(3)*(x[i] + 0.5*dx)) + exp(-sqrt(3)*(x[i] - 0.5*dx)))/dx
        if(0 > x[i] - 0.5*dx ) and (0 < x[i] + 0.5*dx):
            hA[i] = 1 + a1/sqrt(3)*(exp(0) - exp(sqrt(3)*(x[i] - 0.5*dx)) - exp(-sqrt(3)*(x[i] + 0.5*dx)) + exp(0) )/dx
        if (0 > x[i] + 0.5*dx):
            hA[i] = 1 + a1/sqrt(3)*(exp(sqrt(3)*(x[i] + 0.5*dx)) -exp(sqrt(3)*(x[i] - 0.5*dx)))/dx
        
        
        
    return hp,xp,hA

def initialhSmooth(x,dx,h0,h1):
    n = len(x)
    hp = []
    xp = []
    hA = zeros(n)

    for i in range(n):
        xjmh = x[i] - dx/2
        xj = x[i]
        xjph = x[i] + dx/2
        
        hjmh = h0 + h1*exp( - (xjmh)**2)
        hj = h0 + h1*exp( - (xj)**2)
        hjph = h0 + h1*exp( - (xjph)**2)
        

        hcell = [hjmh,hj,hjph]
        xhcell = [xjmh,xj,xjph]
        
        hp = hp + hcell
        xp = xp + xhcell
        
        hAjmh = h0*(xjmh) + 1/2*sqrt(pi)*h1*math.erf(xjmh)
        hAjph = h0*(xjph) + 1/2*sqrt(pi)*h1*math.erf(xjph)
        
        hA[i] = (hAjph - hAjmh)/dx 
        
        
    return hp,xp,hA


def ReconPoly(qa,eps,nG,npc,dx):
    n = len(qa)
    qm = zeros(npc*n)
    
    #Boundaries
    j = 0
    qm[npc*j: npc*(j +nG)] = qa[j]*ones(npc*nG)
    j = n-nG
    qm[npc*j: npc*(j+nG)] = qa[j]*ones(npc*nG)
    
    for j in range(1,n-nG):
        
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

        iw1 = ((1.0/10.0) / (eps +Bjm2toj )**2)
        iw2 = ((3.0/5.0) / (eps +Bjm1tojp1 )**2)
        iw3 = ((3.0/10.0)  / (eps +Bjtojp2 )**2)
        
        w1 = iw1 / (iw1 + iw2 + iw3)
        w2 = iw2 / (iw1 + iw2 + iw3)
        w3 = iw3 / (iw1 + iw2 + iw3)
    
        p2a = w1*pjm2toja + w2*pjm1tojp1a + w3*pjtojp2a
        p2b = w1*pjm2tojb + w2*pjm1tojp1b + w3*pjtojp2b
        p2c = w1*pjm2tojc + w2*pjm1tojp1c + w3*pjtojp2c   
        
        p2jmh = p2a*(-dx/2.0)**2 + p2b*(-dx/2.0) + p2c 
        p2jph = p2a*(dx/2.0)**2 + p2b*(dx/2.0) + p2c
        
        TVDIndLeft = ((qa[j-1] >= qa[j]) and ( qa[j-1] >= p2jmh and p2jmh >= qa[j] )) \
            or ((qa[j-1] <= qa[j]) and ( qa[j-1] <= p2jmh and p2jmh <= qa[j] ))
        TVDIndRight = ((qa[j] >= qa[j+1]) and ( qa[j] >= p2jph and p2jph >= qa[j+1] )) \
            or ((qa[j] <= qa[j+1]) and ( qa[j] <= p2jph and p2jph <= qa[j+1] ))

            
        if (TVDIndLeft and TVDIndRight):
            qm[npc*j] = p2jmh 
            qm[npc*j+1] = p2c
            qm[npc*j+2] = p2jph
        else:
            pjm1toja  =(qa[j] - qa[j-1])/(dx)
            pjm1tojb  =qa[j]
            
            pjtojp1a  =(qa[j+1] - qa[j])/(dx)
            pjtojp1b  =qa[j]
            
            
            Bjm1toj = (qa[j] - qa[j-1])**2
            Bjtojp1 = (qa[j+1] - qa[j])**2
            
            iw1 = ((1.0/3.0) / (eps + Bjm1toj )**2)
            iw2 = ((2.0/3.0) / (eps + Bjtojp1 )**2)
            
            w1 = iw1 / (iw1 + iw2)
            w2 = iw2 / (iw1 + iw2)
        
            p1a = 0
            p1b = w1*pjm1toja + w2*pjtojp1a
            p1c = w1*pjm1tojb + w2*pjtojp1b
            
            p1jmh = p1b*(-dx/2.0) + p1c 
            p1jph = p1b*(dx/2.0) + p1c 
            
            TVDIndLeft = ((qa[j-1] >= qa[j]) and ( qa[j-1] >= p1jmh and p1jmh >= qa[j] )) \
            or ((qa[j-1] <= qa[j]) and ( qa[j-1] <= p1jmh and p1jmh <= qa[j] )) 
            TVDIndRight = ((qa[j] >= qa[j+1]) and ( qa[j] >= p1jph and p1jph >= qa[j+1] )) \
            or ((qa[j] <= qa[j+1]) and ( qa[j] <= p1jph and p1jph <= qa[j+1] ))
    
            if (TVDIndLeft and TVDIndRight):
                qm[npc*j] = p1jmh 
                qm[npc*j+1] = p1c
                qm[npc*j+2] = p1jph
            
            else:
        
                p0a = 0
                p0b = 0
                p0c = qa[j]
                
                p0jmh = p0c 
                p0jph = p0c 
                
                qm[npc*j] = p0c
                qm[npc*j+1] = p0c
                qm[npc*j+2] = p0c


    
    return qm



nG = 3
npc = 3
eps = 10.0**(-10)

expn = 10
lown = 10

g = 9.81

dxs = []
L2s = []
L2phs = []


for expi in range(expn):
    ncurr = lown*(2**expi)
    
    sx = -10.0
    ex = 10.0
    dx = ((ex - sx))/ncurr
    x = arange(sx,ex+ dx,dx)


    hMSA,xS,hAS = initialhDB(x,dx,2,1)
    hMS = ReconPoly(hAS,eps,nG,npc,dx)
    
    # hMSA,xS,hAS= initialPeak(x,dx,0.5)
    # hMS = ReconPoly(hAS,eps,nG,npc,dx)
    
    
    # hMSA,xS,hAS=initialhSmooth(x,dx,2,1)
    # hMS = ReconPoly(hAS,eps,nG,npc,dx)
    
    L2 = norm(hMSA - hMS ,ord=2)/norm(hMS,ord=2)
    L2ph = norm(hMSA[2::3]- hMS[2::3],ord=2)/norm(hMSA[2::3],ord=2)
    
    dxs.append(dx)
    L2s.append(L2)
    L2phs.append(L2ph)

loglog(dxs,L2s,'.b')
loglog(dxs,L2phs,'or')
loglog(dxs,0.1*array(dxs)**1,'-')
loglog(dxs,0.1*array(dxs)**2,'-')
