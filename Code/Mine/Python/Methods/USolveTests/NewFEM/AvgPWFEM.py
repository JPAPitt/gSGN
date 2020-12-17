
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from matplotlib.pyplot import plot, loglog
from numpy.linalg import solve,lstsq,cond,norm

def Recon(qA,x,dx,nG):
    n = len(qA)
    qcubic = zeros(4*n)
    
    # j= 0
    j = 0
    qcubic[4*j: 4*(j+nG)] = qA[j]*ones(4*nG)
    
    for j in range(nG,n-nG):

        P3r0a33,P3r0a32,P3r0a31,P3r0a30 = P3jtojp3(qA[j],qA[j+1],qA[j+2],qA[j+3],dx)
        P3r1a33,P3r1a32,P3r1a31,P3r1a30 = P3jm1tojp2(qA[j],qA[j-1],qA[j+1],qA[j+2],dx)
        P3r2a33,P3r2a32,P3r2a31,P3r2a30 = P3jm2tojp1(qA[j],qA[j-1],qA[j-2],qA[j+1],dx)
        P3r3a33,P3r3a32,P3r3a31,P3r3a30 =  P3jm3toj(qA[j],qA[j-1],qA[j-2],qA[j-3],dx)
        
        P2r0a22,P2r0a21,P2r0a20=   P2jtojp2(qA[j],qA[j+1],qA[j+2],dx)
        P2r1a22,P2r1a21,P2r1a20  =  P2jm1tojp1(qA[j-1],qA[j],qA[j+1],dx)
        P2r2a22,P2r2a21,P2r2a20  =  P2jm2toj(qA[j-2],qA[j-1],qA[j],dx)
        
        P1r0a11,P1r0a10  =  P1jtojp1(qA[j],qA[j+1],dx)
        P1r1a11,P1r1a10  =  P1jm1toj(qA[j],qA[j-1],dx)
        
        P1Wa11,P1Wa10  = P1jWeird(qA[j+1],qA[j],qA[j-1],dx)

        # w1 = minmodL([P1r0a11,P1r1a11,P1Wa11])  / P1Wa11
        # P1a11 = P1Wa11
        P1a11w1 = minmodL([P1r0a11,P1r1a11,P1Wa11])
        
        #w2 = minmodL([P2r0a22,P2r1a2,P2r2a22])  / P2r1a22
        #P2a22 = P2r1a22
        P2a22w2 = minmodL([P2r0a22,P2r1a22,P2r2a22]) 
        
        #w3 = minmodL([P3r0a33,P3r1a33,P3r2a33,P3r3a33])  / P3r2a33
        #P3a33 = P3r2a33
        P3a33w3 = minmodL([P3r0a33,P3r1a33,P3r2a33,P3r3a33])
        
        #Weight on second coefficient
        if ( abs(P3r2a33) > 10.0**(-12.0)):
            w3 = minmodL([P3r0a33,P3r1a33,P3r2a33,P3r3a33])  / P3r2a33
        else:
            w3 = 0
        #P3a32w3 = w3*P3r2a32
        P3a32w3Corr = w3*(P3r2a32 - P2a22w2) + P2a22w2
        
        
        #Weight on first coefficient, want it to reduce to Reconed P2a22
        if (abs(P2r1a22)  > 10.0**(-12.0)):
            w2 = minmodL([P2r0a22,P2r1a22,P2r2a22])  / P2r1a22
        else:
            w2 = 0
        P2a21w2 = w2 *P2r1a21
        
        P2a21w2Corr = w2*(P2r1a21 - P1a11w1) + P1a11w1
        
        P3a31w3Corr = w3*(P3r2a31 - P2a21w2Corr) + P2a21w2Corr
        
        na33 = P3a33w3
        na32 = P3a32w3Corr
 
        na31 = P3a31w3Corr
        na30 = qA[j] - na32/3*(dx/2)**2 
        

        qcubic[4*j] = na33*(-dx/2)**3 + na32*(-dx/2)**2 + na31*(-dx/2) + na30
        qcubic[4*j + 1] = na33*(-dx/6)**3 + na32*(-dx/6)**2 + na31*(-dx/6) + na30
        qcubic[4*j + 2] = na33*(dx/6)**3 + na32*(dx/6)**2 + na31*(dx/6) + na30
        qcubic[4*j + 3] = na33*(dx/2)**3 + na32*(dx/2)**2 + na31*(dx/2) + na30


    # j= n-3
    j = n-nG
    qcubic[4*j: 4*(j+nG)] = qA[j]*ones(4*nG)
        
    return qcubic


def PointsFromPolyFromPoints(a0,a1,a2,a3,x):
    return a0 + a1*x + a2*x**2 + a3*x**3


def PolyFromPoints(yjmh,yjms,yjps,yjph,dx):
    a3 = (-9*yjmh + 27*yjms - 27*yjps + 9*yjph)/ (2*dx**3)
    a2 = (9*yjmh - 9*yjms - 9*yjps + 9*yjph  )/ (4*dx**2)
    a1 = (yjmh  - 27*yjms  + 27*yjps  - yjph  )/ (8*dx)
    a0 = (-yjmh+ 9*yjms + 9*yjps  - yjph)/ 16
    return a3,a2,a1,a0

def PolyFromRep(ujmh,ujph,dujmh,dujph,dx):
    a3 = (dx*dujph + dx*dujmh - 2*ujph + 2*ujmh)/dx**3
    a2 = (dujph- dujmh)/(2*dx)
    a1 = (-dx*dujph - dx*dujmh + 6*ujph - 6*ujmh)/(4*dx)
    a0 = -dx*dujph/8 + dx*dujmh/8 + ujph/2 + ujmh/2
    
    return a3,a2,a1,a0

def GetCellAverage(n,q,dx):
    
    qn = zeros(n)
    #  quadratic passes through qjmh,qj,qjph
    for i in range(n):
        qjmh = q[3*i]
        qj = q[3*i + 1]
        qjph = q[3*i + 2]
        pa = 2*(qjph - 2*qj + qjmh)/ (dx**2)
        pb = (qjph - qjmh) / dx
        pc = qj

        qn[i] = 1.0/dx*(2*pa/3*(dx/2.0)**3 + 2*pc*(dx/2.0) )
    return qn

def InitialSolitonPointwise(x,dx,t,ga,a0,a1,b1):
    n = len(x)
    hp = []
    Gp = []
    up = []
    xp = []
    
    k = sqrt(3*a1) / (2*a0*sqrt(a0 + a1))
    c = sqrt(ga*(a0 + a1))
    for i in range(n):
        
        
        xjmh = x[i] - dx/2
        xjms = x[i] - dx/6
        xjps = x[i] + dx/6
        xjph = x[i] + dx/2
        
        #xjmh
        phi = xjmh - c*t
        sechkphi = 1.0 / cosh(k*phi)
        hjmh = a0 + a1*sechkphi**2
        ujmh = c*(1 - a0 /  hjmh)
        Gjmh = ujmh*hjmh + b1*a0*a1*c*(k**2)*sechkphi**4*hjmh - 2*b1*a0*a1**2*c*k**2*sechkphi**4*tanh(k*phi)**2  - 2*b1*a0*a1*c*k**2*sechkphi**2*hjmh*tanh(k*phi)**2 

        phi = xjms - c*t
        sechkphi = 1.0 / cosh(k*phi)
        hjms = a0 + a1*sechkphi**2
        ujms = c*(1 - a0 /  hjms)
        Gjms = ujms*hjms + b1*a0*a1*c*(k**2)*sechkphi**4*hjms - 2*b1*a0*a1**2*c*k**2*sechkphi**4*tanh(k*phi)**2  - 2*b1*a0*a1*c*k**2*sechkphi**2*hjms*tanh(k*phi)**2       
        

        phi = xjps - c*t
        sechkphi = 1.0 / cosh(k*phi)
        hjps = a0 + a1*sechkphi**2
        ujps = c*(1 - a0 /  hjps)
        Gjps = ujps*hjps + b1*a0*a1*c*(k**2)*sechkphi**4*hjps - 2*b1*a0*a1**2*c*k**2*sechkphi**4*tanh(k*phi)**2  - 2*b1*a0*a1*c*k**2*sechkphi**2*hjps*tanh(k*phi)**2  
 
        phi = xjph - c*t
        sechkphi = 1.0 / cosh(k*phi)
        hjph = a0 + a1*sechkphi**2
        ujph = c*(1 - a0 /  hjph)
        Gjph = ujph*hjph + b1*a0*a1*c*(k**2)*sechkphi**4*hjph - 2*b1*a0*a1**2*c*k**2*sechkphi**4*tanh(k*phi)**2  - 2*b1*a0*a1*c*k**2*sechkphi**2*hjph*tanh(k*phi)**2  
       

        hcell = [hjmh,hjms,hjps,hjph]
        ucell = [ujmh,ujms,ujps,ujph]
        Gcell = [Gjmh,Gjms,Gjps,Gjph]
        xhcell = [xjmh,xjms,xjps,xjph]
        hp = hp + hcell
        xp = xp + xhcell
        up = up + ucell
        Gp = Gp + Gcell
        
        
    return array(hp),array(Gp),array(up),array(xp)



def FEMElem(hp,Gp,beta1,dx,j):
    
    GMatE = zeros((1,4))
    AMatE = zeros((4,4))
    
    #j 
    hjmh = hp[4*j]
    hjms = hp[4*j+1]
    hjps = hp[4*j+2]
    hjph = hp[4*j+3]
    
    Gjmh = Gp[4*j]
    Gjms = Gp[4*j+1]
    Gjps = Gp[4*j+2]
    Gjph = Gp[4*j+3]   


    #G matrix
    GMatE[0][0] = Gjmh
    GMatE[0][1] = Gjms
    GMatE[0][2] = Gjps
    GMatE[0][3] = Gjph
    
    # A matrix - row 0
    GjmhujmhC = -435*beta1*hjmh**3/(8*dx**2) + 297*beta1*hjmh**2*hjms/(4*dx**2) + 33*beta1*hjmh**2*hjph/(4*dx**2) - 297*beta1*hjmh**2*hjps/(8*dx**2) + hjmh
    GjmhujmsC = 387*beta1*hjmh**3/(4*dx**2) - 243*beta1*hjmh**2*hjms/(2*dx**2) - 27*beta1*hjmh**2*hjph/(2*dx**2) + 243*beta1*hjmh**2*hjps/(4*dx**2)
    GjmhujpsC = -441*beta1*hjmh**3/(8*dx**2) + 243*beta1*hjmh**2*hjms/(4*dx**2) + 27*beta1*hjmh**2*hjph/(4*dx**2) - 243*beta1*hjmh**2*hjps/(8*dx**2)
    GjmhujphC = 51*beta1*hjmh**3/(4*dx**2) - 27*beta1*hjmh**2*hjms/(2*dx**2) - 3*beta1*hjmh**2*hjph/(2*dx**2) + 27*beta1*hjmh**2*hjps/(4*dx**2)
    
    AMatE[0][0] = GjmhujmhC
    AMatE[0][1] = GjmhujmsC
    AMatE[0][2] = GjmhujpsC
    AMatE[0][3] = GjmhujphC  

    # A matrix - row 1
    GjmsujmhC = -3*beta1*hjmh*hjms**2/(2*dx**2) - 27*beta1*hjms**3/(4*dx**2) - 3*beta1*hjms**2*hjph/(4*dx**2) + 9*beta1*hjms**2*hjps/(2*dx**2)
    GjmsujmsC = -9*beta1*hjmh*hjms**2/(4*dx**2) + 45*beta1*hjms**3/(8*dx**2) - 9*beta1*hjms**2*hjph/(8*dx**2) + 27*beta1*hjms**2*hjps/(4*dx**2) + hjms
    GjmsujpsC = 9*beta1*hjmh*hjms**2/(2*dx**2) + 9*beta1*hjms**3/(4*dx**2) + 9*beta1*hjms**2*hjph/(4*dx**2) - 27*beta1*hjms**2*hjps/(2*dx**2)
    GjmsujphC = -3*beta1*hjmh*hjms**2/(4*dx**2) - 9*beta1*hjms**3/(8*dx**2) - 3*beta1*hjms**2*hjph/(8*dx**2) + 9*beta1*hjms**2*hjps/(4*dx**2)
    
    AMatE[1][0] = GjmsujmhC
    AMatE[1][1] = GjmsujmsC
    AMatE[1][2] = GjmsujpsC
    AMatE[1][3] = GjmsujphC  

    # A matrix - row 2
    GjpsujmhC = -3*beta1*hjmh*hjps**2/(8*dx**2) + 9*beta1*hjms*hjps**2/(4*dx**2) - 3*beta1*hjph*hjps**2/(4*dx**2) - 9*beta1*hjps**3/(8*dx**2)
    GjpsujmsC = 9*beta1*hjmh*hjps**2/(4*dx**2) - 27*beta1*hjms*hjps**2/(2*dx**2) + 9*beta1*hjph*hjps**2/(2*dx**2) + 9*beta1*hjps**3/(4*dx**2)
    GjpsujpsC = -9*beta1*hjmh*hjps**2/(8*dx**2) + 27*beta1*hjms*hjps**2/(4*dx**2) - 9*beta1*hjph*hjps**2/(4*dx**2) + 45*beta1*hjps**3/(8*dx**2) + hjps
    GjpsujphC = -3*beta1*hjmh*hjps**2/(4*dx**2) + 9*beta1*hjms*hjps**2/(2*dx**2) - 3*beta1*hjph*hjps**2/(2*dx**2) - 27*beta1*hjps**3/(4*dx**2)
    
    AMatE[2][0] = GjpsujmhC
    AMatE[2][1] = GjpsujmsC
    AMatE[2][2] = GjpsujpsC
    AMatE[2][3] = GjpsujphC  

    # A matrix - row 3
    GjphujmhC = -3*beta1*hjmh*hjph**2/(2*dx**2) + 27*beta1*hjms*hjph**2/(4*dx**2) + 51*beta1*hjph**3/(4*dx**2) - 27*beta1*hjph**2*hjps/(2*dx**2)
    GjphujmsC = 27*beta1*hjmh*hjph**2/(4*dx**2) - 243*beta1*hjms*hjph**2/(8*dx**2) - 441*beta1*hjph**3/(8*dx**2) + 243*beta1*hjph**2*hjps/(4*dx**2)
    GjphujpsC = -27*beta1*hjmh*hjph**2/(2*dx**2) + 243*beta1*hjms*hjph**2/(4*dx**2) + 387*beta1*hjph**3/(4*dx**2) - 243*beta1*hjph**2*hjps/(2*dx**2)
    GjphujphC = 33*beta1*hjmh*hjph**2/(4*dx**2) - 297*beta1*hjms*hjph**2/(8*dx**2) - 435*beta1*hjph**3/(8*dx**2) + 297*beta1*hjph**2*hjps/(4*dx**2) + hjph
    
    AMatE[3][0] = GjphujmhC
    AMatE[3][1] = GjphujmsC
    AMatE[3][2] = GjphujpsC
    AMatE[3][3] = GjphujphC  
    
    return GMatE,AMatE


def FEM(hp,Gp,u0,u1,beta1,n,dx,nG):
    
    nbG = 4*(nG )
    nb = 4*(n)
    A = eye(nb,nb)
    b = zeros(nb)
    
    # print(n,nbG,nb)

    
    for i in range(nG,n-nG):
    #    #elementwisematrices
        GMatE,AMatE = FEMElem(hp,Gp,beta1,dx,i)
        
        A[4*i : 4*(i+1), 4*i : 4*(i+1)] = AMatE
        b[4*i : 4*(i+1) ]= GMatE

 
    # # i = 0
    i=0
    # A[i:i+nbG, i:i+nbG] =  eye(nbG,nbG)
    b[i:i+nbG] =  u0*ones(nbG)
 
    # # i = 0
    i=nb-nbG
    # A[i:i+nbG, i:i+nbG] =  eye(nbG,nbG)
    b[i:i+nbG] =  u1*ones(nbG)
    # # print(A)
    # print(b)
    

    uN = solve(A,b)
    
    # print(A)
    # print(b)
    return A,b,uN 

def AverageEdges(q,xn,dx):
    n = len(xn)
    nq = zeros(2*(n+1))
    nx = zeros(n+1)
    
    nq[0] = q[0]
    nq[1] = q[1]
    
    nx[0] = xn[0] - 0.5*dx
    
    for i in range(n-1):
        nq[2*i + 2] = 0.5*(q[4*i + 2] + q[4*i + 4])
        nq[2*i + 3] = 0.5*(q[4*i + 3] + q[4*i +5])
        nx[i+1] = xn[i+1] - 0.5*dx

    #n
    i =n-1
    nq[2*i + 2] = 0.5*(q[4*i + 2])
    nq[2*i + 3] = 0.5*(q[4*i + 3])
    nx[i+1] = xn[i] + 0.5*dx
        
    return nq,nx

L2us = []
L2dus = []
dxs = []    
# for ij in range(10):
for ij in range(0,1):
    sx = -50
    ex = 50
    n= 10*2**ij
    # n=4
    nG = 3
    dx = (ex- sx)/n
    xn = arange(sx - nG*dx,ex+ (nG+ 1)*dx,dx)
    
    ga = 10.0
    a0 = 1.0
    a1 = 0.1
    t = 0
    # b1 = 0
    b1 = 2.0/3.0
    
    h,G,u,x = InitialSolitonPointwise(xn,dx,t,ga,a0,a1,b1)
    
    
    u0 = 0
    du0 = 0
    u1 = 0
    du1 = 0
    
    A,b,NuN = FEM(h,G,u0,u1,b1,len(xn),dx,nG )
    
    L2u = norm(u - NuN,ord=2)/norm(u,ord=2)
    
    L2us.append(L2u)
    dxs.append(dx)
    
loglog(dxs,L2us,'.b')
loglog(dxs,0.01*array(dxs)**2)
# # # j= 10
# j = n//2
# be,Ae = FEMElem(Nh,NG,b1,dx,j)

# Nucell= Nu[4*j:4*(j+1)]
# NGcells= NG[4*(j-1):4*(j+2)]

# Gcell = dot(Ae,Nucell)

# Gs = dot(A,NuC)

# bs = lstsq(A,Gs,rcond=0.0001)


# Ga  = GetCellAverage(len(xn),G,dx)
# # un = FEM(b1,h,Ga,len(xn),x,dx)

# A,b = FEMforG(xn,h,Ga,u[0],u[-1],nG,dx)

# i = n//2
# ha0,ha1,ha2,ha3 = PolyFromPoints(h[4*i],h[4*i+1],h[4*i+2],h[4*i+3],dx)
# ua0,ua1,ua2,ua3 = PolyFromPoints(u[4*i],u[4*i+1],u[4*i+2],u[4*i+3],dx)
# Gjmh, Gjms, Gjps,  Gjph = CalculateG(b1,ha0,ha1,ha2,ha3,ua0,ua1,ua2,ua3 ,dx)
# print(Gjmh, Gjms, Gjps,  Gjph)
# print(G[4*i:4*(i+1)])
# G = FEMforG(xn,h,u,dx)