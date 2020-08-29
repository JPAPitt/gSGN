
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from numpy.linalg import solve,norm
from scipy.optimize import bisect 
from matplotlib.pyplot import plot,loglog,show,legend,title,figure


def InitialCellAvgFlatLinear(x,x0,x1,a0,a1,dx):
    
    n = len(x)
    qA = zeros(n)
    for i in range(n):
        
        if x[i] < x0:
            qA[i] = a0 
            
        elif x[i] > x1:
            qA[i] = a1
        else:
            qA[i] = ((a1 - a0) / (x1 - x0)) * (x[i] - x0) + a0


        
    return qA
    
    
    
def InitialCellAvgSmoothPeak(x,a0,a1,dx):
    
    n = len(x)
    qA = zeros(n)
    for i in range(n):

        xph = x[i] + 0.5*dx
        xmh = x[i] - 0.5*dx
        qAph = (a0*xph + a1/3*xph**3)
        qAmh = (a0*xmh + a1/3*xmh**3)
        qA[i] = (qAph - qAmh) /dx

        
    return qA

def InitialCellAvgSpike(x,a0,a1,dx):
    
    n = len(x)
    qA = zeros(n)
    for i in range(n):
        
        
        if (x[i]> 0):
            qA[i] = a0 
        elif (x[i] < 0 and x[i] > -dx):
            qA[i] = a0 + a1
        else:
            qA[i] = a0 
        
    return qA

def InitialCellAvgQNonFlatJump(x,a0,a1,dx):
    
    n = len(x)
    qA = zeros(n)
    for i in range(n):
        
        
        if (x[i]> 0):
            qA[i] = a0 + a1*x[i]
        else:
            qA[i] = 2*a0 + a1*x[i]
        
    return qA

def InitialCellAvgQCusp(x,a0,a1,dx):
    
    n = len(x)
    qA = zeros(n)
    for i in range(n):
        
        
        
        if (0 < x[i] - 0.5*dx):
            qA[i] = a0 + a1*(-exp(-(x[i] + 0.5*dx)) + exp(-(x[i] - 0.5*dx)))/dx
        if(0 > x[i] - 0.5*dx ) and (0 < x[i] + 0.5*dx):
            qA[i] = a0 + a1*(exp(0) - exp(x[i] - 0.5*dx) - exp(-(x[i] + 0.5*dx)) + exp(0) )/dx
        if (x0 > x[i] + 0.5*dx):
            qA[i] = a0 + a1*(exp(x[i] + 0.5*dx) -exp(x[i] - 0.5*dx))/dx
        
    return qA


def InitialCellAvgQStep(x,x0,a0,a1,dx):
    
    n = len(x)
    qA = zeros(n)
    for i in range(n):
        
        if (x0 < x[i] - 0.5*dx):
            qA[i] = a1
        if(x0 > x[i] - 0.5*dx ) and (x0 < x[i] + 0.5*dx):
            qA[i] = ( a0*(x0 - (x[i] - 0.5*dx) )  + a1* ((x[i] + 0.5*dx)  - x0) )/dx
            
        if (x0 > x[i] + 0.5*dx):
            qA[i] = a0
        
    return qA

def ReconPerfectKnowledge(qA,x,x0,nBC,dx):
    n = len(qA)
    qe = zeros(n-2*nBC)
    xe = zeros(n-2*nBC)
    
    for j in range(nBC,n-nBC):
        
        i = j - nBC
        
        xe[i] = x[j] + 0.5*dx

        
        pjm3toja  = (qA[j-1] - 3*qA[j-2] + 3*qA[j-3] - qA[j-4])/(6*dx**3)
        pjm3tojb  = (3*qA[j-1] - 8*qA[j-2] + 7*qA[j-3] - 2*qA[j-4])/(2*dx**2)
        pjm3tojc  = (103*qA[j-1] - 225*qA[j-2] + 165*qA[j-3] - 43*qA[j-4])/(24*dx)
        pjm3tojd  =  31*qA[j-1]/8 - 17*qA[j-2]/3 + 89*qA[j-3]/24 - 11*qA[j-4]/12

 
        pjtojp3a = (-qA[j+1] + 3*qA[j+2] - 3*qA[j+3] + qA[j+4])/(6*dx**3)
        pjtojp3b = (3*qA[j+1] - 8*qA[j+2] + 7*qA[j+3] - 2*qA[j+4])/(2*dx**2)
        pjtojp3c =  (-103*qA[j+1] + 225*qA[j+2] - 165*qA[j+3] + 43*qA[j+4])/(24*dx)
        pjtojp3d  = 31*qA[j+1]/8 - 17*qA[j+2]/3 + 89*qA[j+3]/24 - 11*qA[j+4]/12

        #Perfect knowledge of discontinuity locations
        #Just need to know if we are left or right of it
        if xe[i] < x0:
            qa = pjm3toja
            qb = pjm3tojb
            qc = pjm3tojc
            qd = pjm3tojd
            
        else:
            qa = pjtojp3a
            qb = pjtojp3b
            qc = pjtojp3c
            qd = pjtojp3d
        
        qe[i] = qa*(dx/2)**3 + qb*(dx/2)**2 + qc*(dx/2) + qd

    return qe,xe


def ReconComp(qA,x,x0,nBC,dx):
    n = len(qA)
    qem = zeros(n-2*nBC)
    qec = zeros(n-2*nBC)
    qep = zeros(n-2*nBC)    
    for j in range(nBC,n-nBC):
        
        i = j - nBC
                
        pjm3toja  = (qA[j-1] - 3*qA[j-2] + 3*qA[j-3] - qA[j-4])/(6*dx**3)
        pjm3tojb  = (3*qA[j-1] - 8*qA[j-2] + 7*qA[j-3] - 2*qA[j-4])/(2*dx**2)
        pjm3tojc  = (103*qA[j-1] - 225*qA[j-2] + 165*qA[j-3] - 43*qA[j-4])/(24*dx)
        pjm3tojd  =  31*qA[j-1]/8 - 17*qA[j-2]/3 + 89*qA[j-3]/24 - 11*qA[j-4]/12
        
        pjm1tojp2a = (-3*qA[j+1] + qA[j+2] - qA[j-1] + 3*qA[j])/(6*dx**3)
        pjm1tojp2b = (qA[j+1] + qA[j-1] - 2*qA[j])/(2*dx**2)
        pjm1tojp2c = (27*qA[j+1] - 5*qA[j+2] - 7*qA[j-1] - 15*qA[j])/(24*dx)
        pjm1tojp2d =-qA[j+1]/24 - qA[j-1]/24 + 13*qA[j]/12
        
        pjtojp3a = (-qA[j+1] + 3*qA[j+2] - 3*qA[j+3] + qA[j+4])/(6*dx**3)
        pjtojp3b = (3*qA[j+1] - 8*qA[j+2] + 7*qA[j+3] - 2*qA[j+4])/(2*dx**2)
        pjtojp3c =  (-103*qA[j+1] + 225*qA[j+2] - 165*qA[j+3] + 43*qA[j+4])/(24*dx)
        pjtojp3d  = 31*qA[j+1]/8 - 17*qA[j+2]/3 + 89*qA[j+3]/24 - 11*qA[j+4]/12
        
        qa = pjm3toja
        qb = pjm3tojb
        qc = pjm3tojc
        qd = pjm3tojd
        
        qem[i] = qa*(dx/2)**3 + qb*(dx/2)**2 + qc*(dx/2) + qd
        
        qa = pjm1tojp2a
        qb = pjm1tojp2b
        qc = pjm1tojp2c
        qd = pjm1tojp2d      
        qec[i] = qa*(dx/2)**3 + qb*(dx/2)**2 + qc*(dx/2) + qd
 
        qa = pjtojp3a
        qb = pjtojp3b
        qc = pjtojp3c
        qd = pjtojp3d       
        qep[i] = qa*(dx/2)**3 + qb*(dx/2)**2 + qc*(dx/2) + qd
            
    return qem,qec,qep

def DiscontinuityWeighting(qA,x,x0,nBC,dx):
    n = len(qA)
    qfd = zeros(n-2*nBC)
    qbd = zeros(n-2*nBC)
    qcd = zeros(n-2*nBC)
    
    dqem = zeros(n-2*nBC)
    dqec = zeros(n-2*nBC)
    dqep = zeros(n-2*nBC)
    
    ddqem = zeros(n-2*nBC)
    ddqec = zeros(n-2*nBC)
    ddqep = zeros(n-2*nBC)
    
    
    for j in range(nBC,n-nBC):
        
        i = j - nBC
        

        qfd[i] = (qA[j+1] - qA[j]) / dx
        qcd[i] = (qA[j+1] - qA[j-1]) / (2*dx)
        qbd[i] = (qA[j] - qA[j-1]) / dx
        
        pjm3toja  = (qA[j-1] - 3*qA[j-2] + 3*qA[j-3] - qA[j-4])/(6*dx**3)
        pjm3tojb  = (3*qA[j-1] - 8*qA[j-2] + 7*qA[j-3] - 2*qA[j-4])/(2*dx**2)
        pjm3tojc  = (103*qA[j-1] - 225*qA[j-2] + 165*qA[j-3] - 43*qA[j-4])/(24*dx)
        pjm3tojd  =  31*qA[j-1]/8 - 17*qA[j-2]/3 + 89*qA[j-3]/24 - 11*qA[j-4]/12
        
        pjm1tojp2a = (-3*qA[j+1] + qA[j+2] - qA[j-1] + 3*qA[j])/(6*dx**3)
        pjm1tojp2b = (qA[j+1] + qA[j-1] - 2*qA[j])/(2*dx**2)
        pjm1tojp2c = (27*qA[j+1] - 5*qA[j+2] - 7*qA[j-1] - 15*qA[j])/(24*dx)
        pjm1tojp2d =-qA[j+1]/24 - qA[j-1]/24 + 13*qA[j]/12
        
        pjtojp3a = (-qA[j+1] + 3*qA[j+2] - 3*qA[j+3] + qA[j+4])/(6*dx**3)
        pjtojp3b = (3*qA[j+1] - 8*qA[j+2] + 7*qA[j+3] - 2*qA[j+4])/(2*dx**2)
        pjtojp3c =  (-103*qA[j+1] + 225*qA[j+2] - 165*qA[j+3] + 43*qA[j+4])/(24*dx)
        pjtojp3d  = 31*qA[j+1]/8 - 17*qA[j+2]/3 + 89*qA[j+3]/24 - 11*qA[j+4]/12
        
        qa = pjm3toja
        qb = pjm3tojb
        qc = pjm3tojc
        qd = pjm3tojd
        
        dqem[i] = 3*qa*(dx/2)**2 + 2*qb*(dx/2) + qc
        ddqem[i] = 6*qa*(dx/2) + 2*qb
        
        qa = pjm1tojp2a
        qb = pjm1tojp2b
        qc = pjm1tojp2c
        qd = pjm1tojp2d      
        dqec[i] =3*qa*(dx/2)**2 + 2*qb*(dx/2) + qc
        ddqec[i] = 6*qa*(dx/2) + 2*qb
        
        qa = pjtojp3a
        qb = pjtojp3b
        qc = pjtojp3c
        qd = pjtojp3d       
        dqep[i] = 3*qa*(dx/2)**2 + 2*qb*(dx/2) + qc
        ddqep[i] = 6*qa*(dx/2) + 2*qb
            
    return dqem,dqec,dqep,ddqem,ddqec,ddqep,qfd,qcd,qbd

def Poly(x,pa,pb,pc,pd):
    
    return pa*x**3 + pb*x**2 + pc*x + pd

def IntersectionPoly(qA,nBC,dx):
    n = len(qA)
    xlocs = []
    xlinlocs = []
    
    indval_cubic = zeros(n-2*nBC)
    indval_lin  = zeros(n-2*nBC)
    
    for j in range(nBC,n-nBC):
        
        i = j - nBC
        
        linpap = (qA[j+2] - qA[j+1]) / dx
        linpbp = qA[j+1]
        linpam = (qA[j-1] - qA[j-2]) / dx
        linpbm = qA[j-1]
        
        
        if (linpap - linpam != 0):
            xintloc = x[j] + (((linpbp - dx*linpap) - (linpbm + dx*linpam))/ (linpap - linpam) )
            if(xintloc > x[j] + 0.5*dx):
                indval_lin[i] = 1.0
            else:
                indval_lin[i] = -1.0
        
        pjm3toja  = (qA[j-1] - 3*qA[j-2] + 3*qA[j-3] - qA[j-4])/(6*dx**3)
        pjm3tojb  = (3*qA[j-1] - 8*qA[j-2] + 7*qA[j-3] - 2*qA[j-4])/(2*dx**2)
        pjm3tojc  = (103*qA[j-1] - 225*qA[j-2] + 165*qA[j-3] - 43*qA[j-4])/(24*dx)
        pjm3tojd  =  31*qA[j-1]/8 - 17*qA[j-2]/3 + 89*qA[j-3]/24 - 11*qA[j-4]/12
        
        
        pjtojp3a = (-qA[j+1] + 3*qA[j+2] - 3*qA[j+3] + qA[j+4])/(6*dx**3)
        pjtojp3b = (3*qA[j+1] - 8*qA[j+2] + 7*qA[j+3] - 2*qA[j+4])/(2*dx**2)
        pjtojp3c =  (-103*qA[j+1] + 225*qA[j+2] - 165*qA[j+3] + 43*qA[j+4])/(24*dx)
        pjtojp3d  = 31*qA[j+1]/8 - 17*qA[j+2]/3 + 89*qA[j+3]/24 - 11*qA[j+4]/12
        
        loc1 = 2.5*dx
        loc2 = -2.5*dx
        diffqhp = (pjm3toja - pjtojp3a)*(loc1)**3 + (pjm3tojb - pjtojp3b)*(loc1)**2 + (pjm3tojc - pjtojp3c)*(loc1 ) + (pjm3tojd - pjtojp3d )
        diffqmp = (pjm3toja - pjtojp3a)*(loc2)**3 + (pjm3tojb - pjtojp3b)*(loc2)**2 + (pjm3tojc - pjtojp3c)*(loc2) + (pjm3tojd - pjtojp3d )
        if (sign(diffqhp) != sign(diffqmp) ):
            currxloc = bisect(Poly,loc2,loc1,args=(pjm3toja - pjtojp3a ,pjm3tojb - pjtojp3b, pjm3tojc - pjtojp3c, pjm3tojd - pjtojp3d     ))
        
            xlocs.append(x[j]+ currxloc)
            
            
    return xlocs,xlinlocs,indval_lin


def SlopeLogic(signslopejm1b ,signslopejf,signslopejb,signslopejp1f ):
    
    if( (signslopejm1b ==signslopejf) and (signslopejf ==  signslopejb) and (signslopejb == signslopejp1f)):
        return 0
    else:
        sign01 = (signslopejm1b == signslopejb )
        sign23 = (signslopejf ==  signslopejp1f)
        sign12 = (signslopejb == signslopejf )
        print(sign01,sign23,sign12)
        if (sign01 and sign12 and not sign23):
            return -1
        elif (sign01 and not sign12 and sign23):
            return 1
        elif (not sign01 and  sign12 and sign23):
            return 1
        elif (not sign01 and not sign12 and sign23):
            return 1
        elif (not sign01 and sign12 and not sign23):
            return 1        
        elif (sign01 and not sign12 and not sign23):
            return -1                 
        elif (not sign01 and not sign12 and not sign23):
            return 1  
    

def OppositeSlope(qA,nBC,dx):
    n = len(qA)
    xlocs = []
    xlinlocs = []

    indval_lin  = zeros(n-2*nBC)
    
    for j in range(nBC,n-nBC):
        
        i = j - nBC

        signslopejm1b = sign((qA[j-1] - qA[j-2]))
        signslopejf = sign((qA[j+1] - qA[j]))
        signslopejc = sign((qA[j+1] - qA[j-1]))
        signslopejb = sign((qA[j] - qA[j-1]))
        signslopejp1f = sign((qA[j+2] - qA[j+1]))
        
        
        signcube1  = sign(qA[j-1] - 3*qA[j-2] + 3*qA[j-3] - qA[j-4])   
        signcube2 = sign(-3*qA[j+1] + qA[j+2] - qA[j-1] + 3*qA[j])
        signcube3 = sign(-qA[j+1] + 3*qA[j+2] - 3*qA[j+3] + qA[j+4])
        

        indval_lin[i] = SlopeLogic(signslopejm1b ,signslopejf,signslopejb,signslopejp1f )
          
    return indval_lin



def BadRegion(qA,nBC,dx):
    n = len(qA)
    indval_lin  = zeros(n-2*nBC)
    secondindval = zeros(n-2*nBC)
    
    for j in range(nBC,n-nBC):
        
        i = j - nBC

        signslopejm1b = sign((qA[j-1] - qA[j-2]))
        signslopejb = sign((qA[j] - qA[j-1]))
        signslopejf = sign((qA[j+1] - qA[j]))
        signslopejp1f = sign((qA[j+2] - qA[j+1]))
        
        if( (signslopejm1b ==signslopejf) and (signslopejf ==  signslopejb) and (signslopejb == signslopejp1f)):
            indval_lin[i] = 0
            secondindval[i] = 0
            #pick central one
        else:
            #bad region
            indval_lin[i] = 1  
            
            signslopejm2b = sign((qA[j-2] - qA[j-3]))
            signslopejp2f = sign((qA[j+3] - qA[j+2]))
            if( ((signslopejm1b == signslopejm2b) and (signslopejm1b == signslopejb)) and not ((signslopejp1f == signslopejf) and (signslopejp2f == signslopejp1f) ) ):
                #left good/ right bad
                secondindval[i] = -1
            elif( not ((signslopejm1b == signslopejm2b) and (signslopejm1b == signslopejb)) and ((signslopejp1f == signslopejf) and (signslopejp2f == signslopejp1f) ) ):

                 #right good/ left bad
                secondindval[i] = 1
                
            elif( ((signslopejm1b == signslopejm2b) and (signslopejm1b == signslopejb)) and ((signslopejp1f == signslopejf) and (signslopejp2f == signslopejp1f) ) ):
                #both left and right good
                
                linpap = (qA[j+2] - qA[j+1]) / dx
                linpbp = qA[j+1]
                linpam = (qA[j-1] - qA[j-2]) / dx
                linpbm = qA[j-1]
            
                xintloc = x[j] + (((linpbp - dx*linpap) - (linpbm + dx*linpam))/ (linpap - linpam) )
                if (xintloc  < x[j] + 0.5*dx):
                    secondindval[i] = 1
                else:
                    secondindval[i] = -1
            else:
                #both left and right bad
                linpap = (qA[j+2] - qA[j+1]) / dx
                linpbp = qA[j+1]
                linpam = (qA[j-1] - qA[j-2]) / dx
                linpbm = qA[j-1]

                xintloc = x[j] + (((linpbp - dx*linpap) - (linpbm + dx*linpam))/ (linpap - linpam) )
                if (xintloc  < x[j] + 0.5*dx):
                    secondindval[i] = 1
                else:
                    secondindval[i] = -1
              
            print(secondindval[i],signslopejm2b,signslopejm1b,signslopejb,signslopejf,signslopejp1f)
          
    return indval_lin,secondindval



def NewtonDivDiff(u,j,k,dx):
    
    #
    if (k ==0):
        return u[j]
    else:
        return  (NewtonDivDiff(u,j+1,k-1,dx) - NewtonDivDiff(u,j,k-1,dx) )
        
def UndividedDiff(uA,i,j):
    
    if (j==1):
        # print(i-1/2, i+1/2, uA[i])
        return uA[i]
    else:
        P1 = UndividedDiff(uA,i+1,j-1)
        P2 = UndividedDiff(uA,i,j-1)

        # print(i,j,i + 1 - 1/2, i +1 + (j-1) -1/2 , i  - 1/2, i + (j-1) -1/2, P1, P2)
        return P1 - P2
   


def ENOStencil(qA,nBC,dx):
    n = len(qA)
    LeftVal  = zeros(n-2*nBC)
    qe  = zeros(n-2*nBC)
    k =4
    
    for j in range(nBC,n-nBC):
        
        i = j - nBC
        
        ####
        # Choose Stencil
        ####
        
        #Initial Stencil
        LeftLoc = 0
        #Number of points to add - 2
        
        # print( ' ')
        # print( ' ')
        for l in range(2,k+1):
            # print(i,l)
            # print('Current')
            # print(j + LeftLoc+1 , j + LeftLoc+1 + l)
            DivDiffCurrent = UndividedDiff(qA, j + LeftLoc+1,l)
            # print('')
            # print('Left')
            # print(j + LeftLoc, j + LeftLoc + l)
            DivDiffLeft = UndividedDiff(qA,j + LeftLoc,l)
            
            if LeftLoc <= -1:
            
                if (1*abs(DivDiffLeft) < abs((DivDiffCurrent))):
                    LeftLoc = LeftLoc -1
                    
            else:
                
                if (abs(DivDiffLeft) < abs((DivDiffCurrent))):
                    LeftLoc = LeftLoc -1
        
        LeftVal[i] =  LeftLoc 
        
        # print( ' ')
        # print( ' ')
        ####
        # Reconstruct Over Stencil
        ####
        
        if LeftLoc == -3:            
            qa =(-3*qA[j-1] + 3*qA[j-2] - qA[j-3] + qA[j])/(6*dx**3)
            qb =(-5*qA[j-1] + 4*qA[j-2] - qA[j-3] + 2*qA[j])/(2*dx**2)
            qc =(-69*qA[j-1] + 33*qA[j-2] - 7*qA[j-3] + 43*qA[j])/(24*dx)
            qd = 5*qA[j-1]/24 - qA[j-2]/6 + qA[j-3]/24 + 11*qA[j]/12
            
        elif LeftLoc == -2:
            
            qa = (qA[j+1] + 3*qA[j-1] - qA[j-2] - 3*qA[j])/(6*dx**3)
            qb = (qA[j+1] + qA[j-1] - 2*qA[j])/(2*dx**2)
            qc = (7*qA[j+1] - 27*qA[j-1] + 5*qA[j-2] + 15*qA[j])/(24*dx)
            qd = -qA[j+1]/24 - qA[j-1]/24 + 13*qA[j]/12   
        
        elif LeftLoc == -1:
            
            qa = (-3*qA[j+1] + qA[j+2] - qA[j-1] + 3*qA[j])/(6*dx**3)
            qb = (qA[j+1] + qA[j-1] - 2*qA[j])/(2*dx**2)
            qc = (27*qA[j+1] - 5*qA[j+2] - 7*qA[j-1] - 15*qA[j])/(24*dx)
            qd = -qA[j+1]/24 - qA[j-1]/24 + 13*qA[j]/12
                        
        else:
            qa = (3*qA[j+1] - 3*qA[j+2] + qA[j+3] - qA[j])/(6*dx**3)
            qb = (-5*qA[j+1] + 4*qA[j+2] - qA[j+3] + 2*qA[j])/(2*dx**2)
            qc = (69*qA[j+1] - 33*qA[j+2] + 7*qA[j+3] - 43*qA[j])/(24*dx)
            qd = 5*qA[j+1]/24 - qA[j+2]/6 + qA[j+3]/24 + 11*qA[j]/12
            
        qe[i] = qa*(dx/2)**3 + qb*(dx/2)**2 + qc*(dx/2) + qd
        
    return LeftVal,qe
        

    # return indval_lin,secondindval


nBC = 5
n = 25
sx = -1.0
ex = 1.0
dx = ((ex - sx))/n
x = arange(sx,ex+ dx,dx)
# x = arange(sx,ex+ dx,dx)

a0 = 2.1
a1 = 1
x0 = 0

# qA = InitialCellAvgFlatLinear(x,-0.3,0.3,a0,a1,dx)
# qA =InitialCellAvgSmoothPeak(x,a0,a1,dx)


# qA =  InitialCellAvgSpike(x,a0,a1,dx)
# qA = InitialCellAvgQStep(x,x0,a0,a1,dx)

# qA = InitialCellAvgQNonFlatJump(x,a0,a1,dx)

qA =InitialCellAvgQCusp(x,1,1,dx)


LeftS,qe = ENOStencil(qA,nBC,dx)
qe_pk,xe_pk = ReconPerfectKnowledge(qA,x,x0,nBC,dx)
# qem,qec, qep = ReconComp(qA,x,x0,nBC,dx)

# dqem,dqec,dqep,ddqem,ddqec,ddqep,qfd,qcd,qbd = DiscontinuityWeighting(qA,x,x0,nBC,dx)

# xind ,x2ind= BadRegion(qA,nBC,dx)

figure()
plot(x,qA,'ok', label = 'averages' )
plot(xe_pk,qe,'*', label = 'Recon' )
plot(xe_pk,LeftS,'.', label = 'Left' )
# plot(xe_pk,RightS,'.', label = 'Right' )
legend()
show()
