
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from numpy.linalg import solve,norm
from scipy.optimize import bisect 
from matplotlib.pyplot import plot,loglog,show,legend,title,figure

tol = 10.0**12
    
def InitialCellAvgSmoothPeak(x,a0,a1,dx):
    
    n = len(x)
    qA = zeros(n)
    qE = zeros(n)
    for i in range(n):
        xph = x[i] + 0.5*dx
        xmh = x[i] - 0.5*dx
        qE[i] = a0 + a1*(xph)**2
        qAph = a0*xph  + a1/3*xph**3
        qAmh = a0*xmh  + a1/3*xmh**3
        qA[i] = (qAph - qAmh)/dx


        
    return qA,qE


def InitialCellAvgQCusp(x,a0,a1,dx):
    
    n = len(x)
    qA = zeros(n)
    qE = zeros(n)
    for i in range(n):
        
        qE[i] = a0 + a1*exp(-abs(x[i] + 0.5*dx))
        if (0 < x[i] - 0.5*dx):
            qA[i] = a0 + a1*(-exp(-(x[i] + 0.5*dx)) + exp(-(x[i] - 0.5*dx)))/dx
        if(0 > x[i] - 0.5*dx ) and (0 < x[i] + 0.5*dx):
            qA[i] = a0 + a1*(exp(0) - exp(x[i] - 0.5*dx) - exp(-(x[i] + 0.5*dx)) + exp(0) )/dx
        if (x0 > x[i] + 0.5*dx):
            qA[i] = a0 + a1*(exp(x[i] + 0.5*dx) -exp(x[i] - 0.5*dx))/dx
        
    return qA,qE


def InitialCellAvgQStep(x,x0,a0,a1,dx):
    
    n = len(x)
    qA = zeros(n)
    qE = zeros(n)
    for i in range(n):
        
        qE[i] = a0*(x[i] + 0.5*dx <= x0) + a1*(x[i]+ 0.5*dx > x0)
        
        if (x0 < x[i] - 0.5*dx):
            qA[i] = a1
        if(x0 > x[i] - 0.5*dx ) and (x0 < x[i] + 0.5*dx):
            qA[i] = ( a0*(x0 - (x[i] - 0.5*dx) )  + a1* ((x[i] + 0.5*dx)  - x0) )/dx
            
        if (x0 > x[i] + 0.5*dx):
            qA[i] = a0
        
    return qA,qE

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


# def checkEqual(lst):
#    return lst[1:] == lst[:-1]

# def BadRegionLogic(qA,nBC,dx,j):
    
#     signslopejm2b = sign(qA[j-2] - qA[j-3])
#     signslopejm1b = sign((qA[j-1] - qA[j-2]))
#     signslopejb = sign((qA[j] - qA[j-1]))
#     signslopejf = sign((qA[j+1] - qA[j]))
#     signslopejp1f = sign((qA[j+2] - qA[j+1]))
#     signslopejp2f = sign(qA[j+3] - qA[j+2])
    
#     AllSame = checkEqual([signslopejm2b,signslopejm1b,signslopejb,signslopejf,signslopejp1f,signslopejp2f ])
    
#     if( AllSame):
#         return 0
#     else:        
#         signslopejm2b = sign(qA[j-2] - qA[j-3])
#         signslopejp2f = sign(qA[j+3] - qA[j+2])
#         if( ((signslopejm1b == signslopejm2b) and (signslopejm1b == signslopejb)) and not ((signslopejp1f == signslopejf) and (signslopejp2f == signslopejp1f) ) ):
#             #left good/ right bad
#             return -1
#         elif( not ((signslopejm1b == signslopejm2b) and (signslopejm1b == signslopejb)) and ((signslopejp1f == signslopejf) and (signslopejp2f == signslopejp1f) ) ):

#              #right good/ left bad
#             return 1
            
#         else:
#             #both left and right bad or both good
#             linpap = (qA[j+2] - qA[j+1]) / dx
#             linpbp = qA[j+1]
#             linpam = (qA[j-1] - qA[j-2]) / dx
#             linpbm = qA[j-1]
            

#             xintloc = (((linpbp - dx*linpap) - (linpbm + dx*linpam))/ (linpap - linpam) )


#             if (xintloc  < 0.5*dx - tol):
#                 return 1
#             else:
#                 return -1


def checkEqual(lst):
   return lst[1:] == lst[:-1]

def BadRegionLogic(qA,nBC,dx,j):
    
    signslopejm2b = sign(qA[j-2] - qA[j-3])
    signslopejm1b = sign((qA[j-1] - qA[j-2]))
    signslopejb = sign((qA[j] - qA[j-1]))
    signslopejf = sign((qA[j+1] - qA[j]))
    signslopejp1f = sign((qA[j+2] - qA[j+1]))
    signslopejp2f = sign(qA[j+3] - qA[j+2])
    
    AllSame = checkEqual([signslopejm2b,signslopejm1b,signslopejb,signslopejf,signslopejp1f,signslopejp2f ])
    
    if( AllSame):
        return 0
    else:        
        signslopejm2b = sign(qA[j-2] - qA[j-3])
        signslopejp2f = sign(qA[j+3] - qA[j+2])
        if( ((signslopejm1b == signslopejm2b) and (signslopejm1b == signslopejb)) and not ((signslopejp1f == signslopejf) and (signslopejp2f == signslopejp1f) ) ):
            #left good/ right bad
            return -1
        elif( not ((signslopejm1b == signslopejm2b) and (signslopejm1b == signslopejb)) and ((signslopejp1f == signslopejf) and (signslopejp2f == signslopejp1f) ) ):

             #right good/ left bad
            return 1
            
        else:
            #both left and right bad or both good
            linpap = (qA[j+2] - qA[j+1]) / dx
            linpbp = qA[j+1]
            linpam = (qA[j-1] - qA[j-2]) / dx
            linpbm = qA[j-1]
        
            if (abs(linpap - linpam) > tol):
                xintloc = (((linpbp - dx*linpap) - (linpbm + dx*linpam))/ (linpap - linpam) )
                if (xintloc  <  0.5*dx):
                    return 1
                else:
                    return -1
            else:
                return 1

def BadRegionRecon(qA,x,x0,nBC,dx):
    n = len(qA)
    qe = zeros(n-2*nBC)
    xe = zeros(n-2*nBC)
    
    for j in range(nBC,n-nBC):
        
        i = j - nBC
        
        xe[i] = x[j] + 0.5*dx
        
        Bias = BadRegionLogic(qA,nBC,dx,j)
        
        # if Bias == 2:
        #     qam  = (qA[j-1] - 3*qA[j-2] + 3*qA[j-3] - qA[j-4])/(6*dx**3)
        #     qbm  = (3*qA[j-1] - 8*qA[j-2] + 7*qA[j-3] - 2*qA[j-4])/(2*dx**2)
        #     qcm  = (103*qA[j-1] - 225*qA[j-2] + 165*qA[j-3] - 43*qA[j-4])/(24*dx)
        #     qdm  =  31*qA[j-1]/8 - 17*qA[j-2]/3 + 89*qA[j-3]/24 - 11*qA[j-4]/12
            
        #     qap = (-qA[j+1] + 3*qA[j+2] - 3*qA[j+3] + qA[j+4])/(6*dx**3)
        #     qbp = (3*qA[j+1] - 8*qA[j+2] + 7*qA[j+3] - 2*qA[j+4])/(2*dx**2)
        #     qcp =  (-103*qA[j+1] + 225*qA[j+2] - 165*qA[j+3] + 43*qA[j+4])/(24*dx)
        #     qdp  = 31*qA[j+1]/8 - 17*qA[j+2]/3 + 89*qA[j+3]/24 - 11*qA[j+4]/12
            
        #     qem = qam*(dx/2)**3 + qbm*(dx/2)**2 + qcm*(dx/2) + qdm
        #     qep = qap*(dx/2)**3 + qbp*(dx/2)**2 + qcp*(dx/2) + qdp
            
        #     qe[i] = 0.5*(qem + qep)
        # else:
            
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

expn = 11
lown = 10

dxs = []
L2err_pks = []
L2err_sls = []
L2err_slpks = []
for i in range(expn):
    nBC = 4
    n = lown*2**i
    sx = -1.0
    ex = 1.0
    dx = ((ex - sx))/n
    # x = arange(sx,ex+ dx,dx)
    x = arange(sx+0.5*dx,ex+ dx,dx)
    
    a0 = 2.1
    a1 = 1
    x0 = 0
    
    # qA,qE = InitialCellAvgQStep(x,x0,a0,a1,dx)
    # qA,qE = InitialCellAvgQCusp(x,1,1,dx)
    qA,qE = InitialCellAvgSmoothPeak(x,a0,a1,dx)
    
    qE_nodes = qE[nBC:n-nBC + 1]
    qe_pk,xe_pk = ReconPerfectKnowledge(qA,x,x0,nBC,dx)
    qe_sl,xe_sl = BadRegionRecon(qA,x,x0,nBC,dx)
    
    L2err_pk = norm(qe_pk - qE_nodes , ord =2 )/ norm(qE_nodes , ord =2 )
    L2err_sl = norm(qe_sl - qE_nodes , ord =2 )/ norm(qE_nodes , ord =2 )
    L2err_pk_sl = norm(qe_sl - qe_pk , ord =2 )/ norm(qe_pk , ord =2 )
    
    dxs.append(dx)
    L2err_pks.append(L2err_pk)
    L2err_sls.append(L2err_sl)
    L2err_slpks.append(L2err_pk_sl)

dxs = array(dxs)
L2err_pks = array(L2err_pks)
L2err_sls = array(L2err_sls)
L2err_slpks = array(L2err_slpks)

figure()
title('Convergence')
loglog(dxs,L2err_pks, '*', label='Perfect knowledge')    
loglog(dxs,L2err_sls , '.', label='Slope Logic')    
loglog(dxs,dxs**2 , '-', label='2nd order')    
loglog(dxs,dxs**4 , '-', label='4th order')    
legend()
show()

# figure()
# title('Plot Comparison')
# plot(x,qA,'.',label='Input')
# plot(x + 0.5*dx, qE,'s',label='Function')
# plot(xe_pk, qe_pk,'s',label='Perfect Knowledge')
# plot(xe_pk, qe_sl,'*',label='Slope Logic')
# legend()
# show()
    


# figure()
# title('Plot Comparison')
# plot(xe_pk, qe_pk -qe_sl ,'s',label='Perfect Knowledge')
# legend()
# show()