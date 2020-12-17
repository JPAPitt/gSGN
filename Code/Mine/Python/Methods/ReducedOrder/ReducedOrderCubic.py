
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
        
    eps = 10.0**(-10)
    # j= 0
    j = 0
    qcubic[4*j: 4*(j+3)] = qA[j]*ones(12)
    
    for j in range(3,n-3):
        
        #lets pick a stencil, and try and determine best order
        #should also be a monotonic reconstruction
        
        # stencil is qA[j-1], qA[j], qA[j+1], qA[j+2]    
        #first poly -  qA[j]
        
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
        
        Errj = abs((qA[j] + qA[j+1] - qA[j-1]) - qA[j+2] )
        ErrorMeasure1 = (Errj)
        
        TVDIndL1 = ((qA[j-1] >= qA[j] ) and (qA[j-1] >= P1jmh and P1jmh >= qA[j]  )) \
        or ((qA[j-1] <= qA[j] ) and (qA[j-1] <= P1jmh and P1jmh <= qA[j]  )) 
        
        TVDIndR1 = ((qA[j] >= qA[j+1] ) and (qA[j] >= P1jph and P1jph >= qA[j+1]  )) \
        or ((qA[j] <= qA[j+1] ) and (qA[j] <= P1jph and P1jph <= qA[j+1]  )) 
        
        #third poly -   
        p2a  = (-2*qA[j] + qA[j-1] + qA[j+1])/(2*dx**2)
        p2b  = (-qA[j-1] + qA[j+1])/(2*dx)
        p2c  = 13*qA[j]/12 - qA[j-1]/24 - qA[j+1]/24
        
        P2jmh = p2a*(-dx/2)**2 + p2b*(-dx/2) + p2c
        P2jms = p2a*(-dx/6)**2 + p2b*(-dx/6) + p2c
        P2jps = p2a*(dx/6)**2 + p2b*(dx/6) + p2c
        P2jph = p2a*(dx/2)**2 + p2b*(dx/2) + p2c
        
        TVDIndL2 = ((qA[j-1] >= qA[j] ) and (qA[j-1] >= P2jmh and P2jmh >= qA[j]  )) \
            or ((qA[j-1] <= qA[j] ) and (qA[j-1] <= P2jmh and P2jmh <= qA[j]  )) 

        TVDIndR2 = ((qA[j] >= qA[j+1] ) and (qA[j] >= P2jph  and P2jph  >= qA[j+1]  )) \
            or ((qA[j] <= qA[j+1] ) and (qA[j] <= P2jph  and P2jph  <= qA[j+1]  )) 
            
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
                
        TVDIndL3 = ((qA[j-1] >= qA[j] ) and (qA[j-1] >= P3jmh and P3jmh>= qA[j]  )) \
        or ((qA[j-1] <= qA[j] ) and (qA[j-1] <= P3jmh and P3jmh <= qA[j]  )) 

        TVDIndR3 = ((qA[j] >= qA[j+1] ) and (qA[j] >= P3jph and P3jph >= qA[j+1]  )) \
        or ((qA[j] <= qA[j+1] ) and (qA[j] <= P3jph and P3jph <= qA[j+1]  )) 
     
        qcubic[4*j] = P3jmh
        qcubic[4*j+1] = P3jms
        qcubic[4*j+2] = P3jps
        qcubic[4*j+3] = P3jph
        
        # print(j, ErrorMeasure0,ErrorMeasure1,ErrorMeasure2,TVDIndL1,TVDIndR1,TVDIndL2,TVDIndR2,TVDIndL3,TVDIndR3)
        # if (ErrorMeasure1 > ErrorMeasure2) and TVDIndL3 and TVDIndR3 :
        #    qcubic[4*j] = P3jmh
        #    qcubic[4*j+1] = P3jms
        #    qcubic[4*j+2] = P3jps
        #    qcubic[4*j+3] = P3jph
        # elif (ErrorMeasure1 > ErrorMeasure2) and TVDIndL2 and TVDIndR2:
        #    qcubic[4*j] = P2jmh
        #    qcubic[4*j+1] = P2jms
        #    qcubic[4*j+2] = P2jps
        #    qcubic[4*j+3] = P2jph
        # elif (ErrorMeasure0 > ErrorMeasure1) and TVDIndL1 and TVDIndR1:
        #    qcubic[4*j] = P1jmh
        #    qcubic[4*j+1] = P1jms
        #    qcubic[4*j+2] = P1jps
        #    qcubic[4*j+3] = P1jph
        # else:
        #    qcubic[4*j] = P0jmh
        #    qcubic[4*j+1] = P0jms
        #    qcubic[4*j+2] = P0jps
        #    qcubic[4*j+3] = P0jph
        
            
            



    # j= n-3
    j = n-3
    qcubic[4*j: 4*(j+3)] = qA[j]*ones(12)
        
    return qcubic


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
    
    # hcub,xcub = initialhDBS(xn,dx,2)
    # hcub,xcub = initialhDB(xn,dx,2,1)
    hcub,xcub = initialhSmooth(xn,dx,1,0.1)
    
    hA = GetCellAverage(nxn,hcub,dx)
    
    hcubN = CellAverageToCubic(hA,xn,dx)

    L2 = norm(hcub - hcubN ,ord=2)/norm(hcub ,ord=2)

    
    dxs.append(dx)
    L2s.append(L2)
    
loglog(dxs,L2s,'.')
loglog(dxs,0.1*array(dxs)**2,'-')
loglog(dxs,0.1*array(dxs)**4,'-')

# plot(xcub,hcub,'-')
# plot(xcub,hcubN,'.')
