
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
    dqE = zeros(n)
    for i in range(n):
        xph = x[i] + 0.5*dx
        xmh = x[i] - 0.5*dx
        qE[i] = a0 + a1*(xph)**2
        dqE[i] = 2*a1*(xph)
        qAph = a0*xph  + a1/3*xph**3
        qAmh = a0*xmh  + a1/3*xmh**3
        qA[i] = (qAph - qAmh)/dx


        
    return qA,qE,dqE


def InitialCellAvgQCusp(x,a0,a1,dx):
    
    n = len(x)
    qA = zeros(n)
    qE = zeros(n)
    dqE = zeros(n)
    for i in range(n):
        
        qE[i] = a0 + a1*exp(-abs(x[i] + 0.5*dx))
        
        if (x[i] + 0.5*dx > 0):
             dqE[i] = -a1*exp(-(x[i] + 0.5*dx))
        else:
            dqE[i] = a1*exp((x[i] + 0.5*dx))
            
        
        if (0 < x[i] - 0.5*dx):
            qA[i] = a0 + a1*(-exp(-(x[i] + 0.5*dx)) + exp(-(x[i] - 0.5*dx)))/dx
        if(0 > x[i] - 0.5*dx ) and (0 < x[i] + 0.5*dx):
            qA[i] = a0 + a1*(exp(0) - exp(x[i] - 0.5*dx) - exp(-(x[i] + 0.5*dx)) + exp(0) )/dx
        if (x0 > x[i] + 0.5*dx):
            qA[i] = a0 + a1*(exp(x[i] + 0.5*dx) -exp(x[i] - 0.5*dx))/dx
        
    return qA,qE,dqE


def InitialCellAvgQStep(x,x0,a0,a1,dx):
    
    n = len(x)
    qA = zeros(n)
    qE = zeros(n)
    dqE= zeros(n)
    for i in range(n):
        
        qE[i] = a0*(x[i] + 0.5*dx <= x0) + a1*(x[i]+ 0.5*dx > x0)
        if (x0 < x[i] - 0.5*dx):
            qA[i] = a1
        if(x0 > x[i] - 0.5*dx ) and (x0 < x[i] + 0.5*dx):
            qA[i] = ( a0*(x0 - (x[i] - 0.5*dx) )  + a1* ((x[i] + 0.5*dx)  - x0) )/dx
            
        if (x0 > x[i] + 0.5*dx):
            qA[i] = a0
        
    return qA,qE,dqE

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
    dqe  = zeros(n-2*nBC)
    k =4
    eps = 10.0**-12
    for j in range(nBC,n-nBC):
        
        i = j - nBC
        
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
        # if j > n/2 - 10 and j < n/2 + 10:
        #    print(w1,w2,w3,w4)
        qa = w1*pjm3toja + w2*pjm2tojp1a + w3*pjm1tojp2a + w4*pjtojp3a
        qb = w1*pjm3tojb + w2*pjm2tojp1b + w3*pjm1tojp2b + w4*pjtojp3b
        qc = w1*pjm3tojc + w2*pjm2tojp1c + w3*pjm1tojp2c + w4*pjtojp3c
        qd = w1*pjm3tojd + w2*pjm2tojp1d + w3*pjm1tojp2d + w4*pjtojp3d
        
        if (j > n/2 - 10) and (j <n/2 + 10):
            print(i,w1,w2,w3,w4)
            print(i,qa,qb,qc,qd)
            
        qe[i] = qa*(dx/2)**3 + qb*(dx/2)**2 + qc*(dx/2) + qd
        dqe[i] = 3*qa*(dx/2)**2 + 2*qb*(dx/2) + qc
        
    return qe,dqe 

expn = 10
lown = 10

dxs = []
L2err_pks = []
L2err_sls = []
L2err_dsls = []
L2err_slpks = []
for i in range(expn):
    nBC = 5
    n = lown*2**i
    sx = -1.0
    ex = 1.0
    dx = ((ex - sx))/n
    # x = arange(sx,ex+ dx,dx)
    x = arange(sx+0.1*dx,ex+ dx,dx)
    
    a0 = 2.1
    a1 = 1
    x0 = 0
    
    qA,qE,dqE = InitialCellAvgQStep(x,x0,a0,a1,dx)
    # qA,qE,dqE = InitialCellAvgQCusp(x,1,1,dx)
    # qA,qE,dqE = InitialCellAvgSmoothPeak(x,a0,a1,dx)
    
    qE_nodes = qE[nBC:n-nBC + 1]
    dqE_nodes = dqE[nBC:n-nBC + 1]
    
    qe_pk,xe_pk = ReconPerfectKnowledge(qA,x,x0,nBC,dx)
    qe_sl,dqe_sl = ENOStencil(qA,nBC,dx)
    
    L2err_dsl = norm(dqe_sl - dqE_nodes , ord =2 )/ norm(dqE_nodes , ord =2 )
    L2err_pk = norm(qe_pk - qE_nodes , ord =2 )/ norm(qE_nodes , ord =2 )
    L2err_sl = norm(qe_sl - qE_nodes , ord =2 )/ norm(qE_nodes , ord =2 )
    L2err_pk_sl = norm(qe_sl - qe_pk , ord =2 )/ norm(qe_pk , ord =2 )
    
    dxs.append(dx)
    L2err_pks.append(L2err_pk)
    L2err_sls.append(L2err_sl)
    L2err_dsls.append(L2err_dsl)
    L2err_slpks.append(L2err_pk_sl)

dxs = array(dxs)
L2err_pks = array(L2err_pks)
L2err_sls = array(L2err_sls)
L2err_dsls = array(L2err_dsls)
L2err_slpks = array(L2err_slpks)

figure()
title('Convergence')
loglog(dxs,L2err_pks, '*', label='Perfect knowledge')    
loglog(dxs,L2err_sls , '.', label='WENO')   
loglog(dxs,L2err_dsls , '.', label='derivative of WENO')   
loglog(dxs,dxs , '-', label='1st order')    
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