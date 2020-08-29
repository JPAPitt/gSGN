
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from numpy.linalg import solve,norm
from scipy.optimize import bisect 
from matplotlib.pyplot import plot,loglog,show,legend,title,figure


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


def ReconSlopeLogic(qA,x,x0,nBC,dx):
    n = len(qA)
    qe = zeros(n-2*nBC)
    xe = zeros(n-2*nBC)
    
    for j in range(nBC,n-nBC):
        
        i = j - nBC
        
        xe[i] = x[j] + 0.5*dx
        
        signslopejm1b = sign((qA[j-1] - qA[j-2]))
        signslopejf = sign((qA[j+1] - qA[j]))
        signslopejb = sign((qA[j] - qA[j-1]))
        signslopejp1f = sign((qA[j+2] - qA[j+1]))

        Bias = SlopeLogic(signslopejm1b ,signslopejf,signslopejb,signslopejp1f )
        
        #left
        if Bias == -1:
            qa  = (qA[j-1] - 3*qA[j-2] + 3*qA[j-3] - qA[j-4])/(6*dx**3)
            qb  = (3*qA[j-1] - 8*qA[j-2] + 7*qA[j-3] - 2*qA[j-4])/(2*dx**2)
            qc  = (103*qA[j-1] - 225*qA[j-2] + 165*qA[j-3] - 43*qA[j-4])/(24*dx)
            qd  =  31*qA[j-1]/8 - 17*qA[j-2]/3 + 89*qA[j-3]/24 - 11*qA[j-4]/12
        
        #right
        elif Bias == 1:
            qa = (-qA[j+1] + 3*qA[j+2] - 3*qA[j+3] + qA[j+4])/(6*dx**3)
            qb = (3*qA[j+1] - 8*qA[j+2] + 7*qA[j+3] - 2*qA[j+4])/(2*dx**2)
            qc =  (-103*qA[j+1] + 225*qA[j+2] - 165*qA[j+3] + 43*qA[j+4])/(24*dx)
            qd  = 31*qA[j+1]/8 - 17*qA[j+2]/3 + 89*qA[j+3]/24 - 11*qA[j+4]/12
        
        #none
        else:
            qa = (-3*qA[j+1] + qA[j+2] - qA[j-1] + 3*qA[j])/(6*dx**3)
            qb = (qA[j+1] + qA[j-1] - 2*qA[j])/(2*dx**2)
            qc = (27*qA[j+1] - 5*qA[j+2] - 7*qA[j-1] - 15*qA[j])/(24*dx)
            qd  = -qA[j+1]/24 - qA[j-1]/24 + 13*qA[j]/12
        
     
        qe[i] = qa*(dx/2)**3 + qb*(dx/2)**2 + qc*(dx/2) + qd

    return qe,xe




def SlopeLogic(signslopejm1b ,signslopejf,signslopejb,signslopejp1f ):
    
    
    if( (signslopejm1b ==signslopejf) and (signslopejf ==  signslopejb) and (signslopejb == signslopejp1f)):
        return 0
    else:
        sign01 = (signslopejm1b == signslopejb )
        sign23 = (signslopejf ==  signslopejp1f)
        sign12 = (signslopejb == signslopejf )
        
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
        # signslopejc = sign((qA[j+1] - qA[j-1]))
        signslopejb = sign((qA[j] - qA[j-1]))
        signslopejp1f = sign((qA[j+2] - qA[j+1]))
        
        signcube1  = sign(qA[j-1] - 3*qA[j-2] + 3*qA[j-3] - qA[j-4])   
        signcube2 = sign(-3*qA[j+1] + qA[j+2] - qA[j-1] + 3*qA[j])
        signcube3 = sign(-qA[j+1] + 3*qA[j+2] - 3*qA[j+3] + qA[j+4])
        

        indval_lin[i] = SlopeLogic(signslopejm1b ,signslopejf,signslopejb,signslopejp1f )
          
    return indval_lin


nBC = 4
n = 25
sx = -1.0
ex = 1.0
dx = ((ex - sx))/n
x = arange(sx,ex+ dx,dx)

a0 = 2.1
a1 = 1
x0 = 0


# qA =  InitialCellAvgSpike(x,a0,a1,dx)
# qA = InitialCellAvgQStep(x,x0,a0,a1,dx)

# qA = InitialCellAvgQNonFlatJump(x,a0,a1,dx)

qA =InitialCellAvgQCusp(x,1,1,dx)

qe_pk,xe_pk = ReconPerfectKnowledge(qA,x,x0,nBC,dx)

qe_sl,xe_sl = ReconSlopeLogic(qA,x,x0,nBC,dx)

plot(x,qA,'.',label='Input')
plot(xe_pk, qe_pk,'s',label='Perfect Knowledge')
plot(xe_pk, qe_sl,'*',label='Slope Logic')
legend()
show()



