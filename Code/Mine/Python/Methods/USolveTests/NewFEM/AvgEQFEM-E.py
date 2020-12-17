
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from matplotlib.pyplot import plot
from numpy.linalg import solve,lstsq,cond

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


def PointwiseToNew(q,xn,dx):
    
    n = len(xn)
    qn = len(q)
    qn = zeros(qn)
    
    for i in range(n):
        
       a3,a2,a1,a0 = PolyFromPoints(q[4*i],q[4*i + 1],q[4*i + 2],q[4*i + 3],dx)
       
       qjmh = q[4*i]
       dqjmh = 3*a3*(-dx/2)**2 + 2*a2*(-dx/2) + a1
       qjph = q[4*i+3]
       dqjph = 3*a3*(dx/2)**2 + 2*a2*(dx/2) + a1
       
       qn[4*i] =qjmh 
       qn[4*i+1] =dqjmh 
       qn[4*i+2] =qjph 
       qn[4*i+3] =dqjph 
       
    return qn

def FEMElem(hR,GR,beta1,dx,j):
    
    GMatE = zeros((1,2))
    AMatE = zeros((2,6))
    
    #j 
    hjmh = hR[4*j]
    dhjmh = hR[4*j + 1]
    hjph = hR[4*j + 2]
    dhjph = hR[4*j + 3]


    Gjmh = GR[4*j]
    dGjmh = GR[4*j + 1]
    Gjph = GR[4*j + 2]
    dGjph = GR[4*j + 3]
    
    aj3,aj2,aj1,aj0 = PolyFromRep(Gjmh,Gjph,dGjmh,dGjph,dx)

    #j+1
    hjp1mh = hR[4*(j+1)]
    dhjp1mh = hR[4*(j+1) + 1]
    hjp1ph = hR[4*(j+1) + 2]
    dhjp1ph = hR[4*(j+1) + 3]


    Gjp1mh = GR[4*(j+1)]
    dGjp1mh = GR[4*(j+1) + 1]
    Gjp1ph = GR[4*(j+1) + 2]
    dGjp1ph = GR[4*(j+1) + 3]
    
    ajp13,ajp12,ajp11,ajp10 = PolyFromRep(Gjp1mh,Gjp1ph,dGjp1mh,dGjp1ph,dx)
    
    Gaj = 2*aj2/3*(dx/2.0)**3 + 2*aj0*(dx/2.0) 
    Gajp1 = 2*ajp12/3*(dx/2.0)**3 + 2*ajp10*(dx/2.0) 
    
    dGaj = Gjph - Gjmh 
    dGajp1 = Gjp1ph - Gjp1mh 

    GMatE[0][0] = (Gaj + Gajp1)
    GMatE[0][1] =  dGaj +dGajp1
    
    # G average Equation
    GajequjmhC = dhjmh*dx**2/12 + dx*hjmh/2
    GajeqdujmhC = beta1*hjmh**3/2 + dx**2*hjmh/12
    GajequjphC = -dhjph*dx**2/12 + dx*hjph/2
    GajeqdujphC = -beta1*hjph**3/2 - dx**2*hjph/12

    Gajp1equjmhC = dhjp1mh*dx**2/12 + dx*hjp1mh/2
    Gajp1eqdujmhC = beta1*hjp1mh**3/2 + dx**2*hjp1mh/12
    Gajp1equjphC = -dhjp1ph*dx**2/12 + dx*hjp1ph/2
    Gajp1eqdujphC = -beta1*hjp1ph**3/2 - dx**2*hjp1ph/12


    print('jmh')
    print(GajequjmhC,GajeqdujmhC)  
    print('jph')
    print(GajequjphC,Gajp1equjmhC)
    print(GajeqdujphC,Gajp1eqdujmhC)
    print('jp1ph')
    print(Gajp1equjphC,Gajp1eqdujphC )
    
    AMatE[0][0] = GajequjmhC
    AMatE[0][1] = GajeqdujmhC
    AMatE[0][2] = GajequjphC +  Gajp1equjmhC
    AMatE[0][3] = GajeqdujphC +  Gajp1eqdujmhC
    AMatE[0][4] = Gajp1equjphC
    AMatE[0][5] = Gajp1eqdujphC 

    # dG/dx average Equation
    dGajequjmhC = -3*beta1*hjmh**3/dx**2 - 3*beta1*hjph**3/dx**2 - hjmh
    dGajeqdujmhC = 3*beta1*dhjmh*hjmh**2/2 - 2*beta1*hjmh**3/dx - beta1*hjph**3/dx
    dGajequjphC = 3*beta1*hjmh**3/dx**2 + 3*beta1*hjph**3/dx**2 + hjph
    dGajeqdujphC = -3*beta1*dhjph*hjph**2/2 - beta1*hjmh**3/dx - 2*beta1*hjph**3/dx

    dGajp1equjmhC = -3*beta1*hjp1mh**3/dx**2 - 3*beta1*hjp1ph**3/dx**2 - hjp1mh
    dGajp1eqdujmhC = 3*beta1*dhjp1mh*hjmh**2/2 - 2*beta1*hjp1mh**3/dx - beta1*hjp1ph**3/dx
    dGajp1equjphC = 3*beta1*hjp1mh**3/dx**2 + 3*beta1*hjp1ph**3/dx**2 + hjp1ph
    dGajp1eqdujphC = -3*beta1*dhjp1ph*hjp1ph**2/2 - beta1*hjp1mh**3/dx - 2*beta1*hjp1ph**3/dx
    
    
    AMatE[1][0] = dGajequjmhC
    AMatE[1][1] = dGajeqdujmhC
    AMatE[1][2] = dGajequjphC  + dGajp1equjmhC
    AMatE[1][3] = dGajeqdujphC  + dGajp1eqdujmhC
    AMatE[1][4] = dGajp1equjphC
    AMatE[1][5] = dGajp1eqdujphC     

    print('djmh')
    print(dGajequjmhC,dGajeqdujmhC)  
    print('djph')
    print(dGajequjphC,dGajp1equjmhC)
    print(dGajeqdujphC,dGajp1eqdujmhC)
    print('djp1ph')
    print(dGajp1equjphC,dGajp1eqdujphC )    
    
    
    # AMatE[0][0] = 0.0
    # AMatE[0][1] = 0.1
    # AMatE[0][2] = 0.2
    # AMatE[0][3] = 0.3
    # AMatE[0][4] = 0.4
    # AMatE[0][5] = 0.5

    
    # AMatE[1][0] = 1.0
    # AMatE[1][1] = 1.1
    # AMatE[1][2] = 1.2
    # AMatE[1][3] = 1.3
    # AMatE[1][4] = 1.4
    # AMatE[1][5] = 1.5    
    
    # GMatE[0][0] =1.0
    # GMatE[0][1] =1.1
    # print(GMatE)
    # print(AMatE)

    
    return GMatE,AMatE


def FEM(hR,GR,uR,u0,du0,u1,du1,beta1,n,dx,nG):
    
    n = len(xn)
    nbG = 2*(nG + 1)
    nb = 2*(n+1)
    A = eye(nb,nb)
    b = zeros(nb)
    
    # print(n,nbG,nb)

    
    for i in range(nG-1,n-nG-1):
    #    #elementwisematrices
        print(i,2*i)
        GMatE,AMatE = FEMElem(hR,GR,beta1,dx,i)
        
        A[2*i : 2*i + 2, 2*i-2 : 2*i + 4] = AMatE#A[2*i : 2*i + 2, 2*i-2 : 2*i + 4]  +   AMatE
        b[2*i: 2*(i+1)] = GMatE #b[2*i : 2*(i+1)]+   GMatE

 
    # # i = 0
    i=0
    # A[i:i+nbG, i:i+nbG] =  eye(nbG,nbG)
    b[i:i+nbG:2] =  u0*ones(nbG//2)
    b[i+1:i+nbG:2] =  du0*ones(nbG//2)
 
    # # i = 0
    i=nb-nbG
    # A[i:i+nbG, i:i+nbG] =  eye(nbG,nbG)
    b[i:i+nbG:2] =  u1*ones(nbG//2)
    b[i+1:i+nbG:2] =  du1*ones(nbG//2) 
    # # print(A)
    # print(b)
    

    uN = solve(A,b)
    
    # print(A)
    # print(b)
    return A,b,uN 

def AverageEdges(q,xn):
    n = len(xn)
    nq = zeros(2*(n+1))
    
    nq[2*0] = q[4*0]
    nq[2*0 + 1] = q[4*0 + 1]
    for i in range(1,n):
        nq[2*i] = 0.5*(q[4*i] + q[4*i - 2])
        nq[2*i+1] = 0.5*(q[4*i+1] + q[4*i - 1])

    
    nq[2*(n+1) -2] = q[4*n -2]
    nq[2*(n+1) -1]  = q[4*n-1]
        
    return nq
    

sx = -50
ex = 50
n= 1000
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
Nh = PointwiseToNew(h,xn,dx)
NG = PointwiseToNew(G,xn,dx)
Nu = PointwiseToNew(u,xn,dx)

NuC = AverageEdges(Nu,xn)
NGC = AverageEdges(NG,xn)

u0 = 0
du0 = 0
u1 = 0
du1 = 0
# A,b = FEM(Nh,NG,Nu,u0,du0,u1,du1,b1,len(xn),dx,nG )
A,b,NuN = FEM(Nh,NG,Nu,u0,du0,u1,du1,b1,len(xn),dx,nG )

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