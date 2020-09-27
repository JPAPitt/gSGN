
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from numpy.linalg import solve,norm
from matplotlib.pyplot import plot,loglog

def sech(x):
    return 2.0/ (exp(x) + exp(-x))

def tanh(x):
    return  (exp(x) - exp(-x))/ (exp(x) + exp(-x))

def Soliton(x,t,dx,g,a0,a1):
    n = len(x)
    h = zeros(n)
    G = zeros(n)
    u = zeros(n)

    c = sqrt(g*(a0 + a1))
    k = sqrt(3*a1) / (2*a0*sqrt(a0 + a1))    
    for i in range(n):

        phi = x[i] - c*t
        sechkphi = sech(k*phi)
        #h
        h[i] = a0 + a1*sechkphi**2
        u[i] = c*(1 - a0 / h[i])
        #u   
        G[i] = u[i]*h[i] + 2.0/3*a0*a1*c*(k**2)*sechkphi**4*h[i] - 4.0/3*a0*a1**2*c*k**2*sechkphi**4*tanh(k*phi)**2 - 4.0/3*a0*a1*c*k**2*sechkphi**2*h[i]*tanh(k*phi)**2 
        
    return h,u,G

def Peakon(x,dx,a0):
    n = len(x)
    h = zeros(n)
    G = zeros(n)
    u = zeros(n)
    
    c = sqrt(1 + a1)
    beta1 = 2.0/3.0
    for i in range(n):
        #h
        h[i] = 1 + a0*exp(-sqrt(3)*abs(x[i]))
        u[i] = c*(1 - a0 / h[i])
        #u   
        
        hjph = 1 + a0*exp(-sqrt(3)*abs(x[i] + 0.5*dx))
        hjmh = 1 + a0*exp(-sqrt(3)*abs(x[i] - 0.5*dx))
        dhjph = -a0*sqrt(3)*sign(x[i] + 0.5*dx)*exp(-sqrt(3)*abs(x[i] + 0.5*dx))
        dhjmh = -a0*sqrt(3)*sign(x[i] - 0.5*dx)*exp(-sqrt(3)*abs(x[i] - 0.5*dx))
        dujph = c*a0*dhjph/ hjph**2
        dujmh = c*a0*dhjmh/ hjmh**2
        G[i] = u[i]*h[i] - beta1/(2.0 *dx)*( (hjph**3)*dujph - (hjmh**3)*dujmh)
        
    return h,u,G

def FDHGforu(h,G,dx,beta1,u0,h0,G0,u1,h1,G1):
    n = len(h)
    A = zeros((n,n))
    RHS = zeros(n)
    
    i = 0
    dhc = (h[i+1] - h0) / (2*dx)
    dh1 = (beta1/2.0)* (h[i]**3)/(dx*dx)
    dh2 = 3.0/2.0*(beta1)* (h[i]**2)*dhc/(2*dx)
    
    A[i,i+1] = 0.0
    A[i,i] = 1.0
    
    RHS[i] = u0
    
    for i in range(1,n-1):
        dhc = (h[i+1] - h[i-1]) / (2*dx)
        dh1 = (beta1/2.0)* (h[i]**3)/(dx*dx)
        dh2 = 3.0/2.0*(beta1)* (h[i]**2)*dhc/(2*dx)
        
        A[i,i-1] = -dh1 + dh2
        A[i,i] = h[i] + 2*dh1
        A[i,i+1] = -dh1 - dh2
    
        RHS[i] = G[i]
    
    i = n-1;
    dhc = (h1 - h[i-1]) / (2*dx)
    dh1 = (beta1/2.0)* (h[i]**3)/(dx**2)
    dh2 = 3.0/2.0*(beta1)* (h[i]**2)*dhc/(2*dx)
    
    A[i,i-1] = 0.0
    A[i,i] = 1.0
    
    RHS[i] = u1
    
    u = solve(A,RHS)
    
    return u


def FDHGforuAvg(h,G,dx,beta1,u0,h0,G0,u1,h1,G1):
    n = len(h)
    A = zeros((n,n))
    RHS = zeros(n)
    
    #Using cell average G equation with
    # int(u,h) = int(u)*int(h)
    # and (du/dx)_{j+1/2} = (u[j+1] - u[j])/ dx
    
    i = 0
    
    A[i,i+1] = 0.0
    A[i,i] = 1.0
    
    RHS[i] = u0
    
    for i in range(1,n-1):
        hjph = 0.5*(h[i] + h[i+1])
        hjmh = 0.5*(h[i] + h[i-1])
        
        A[i,i-1] = - beta1/(2.0*dx)*(hjmh**3 / dx )
        A[i,i] = h[i]  - beta1/(2.0*dx)*( -hjph**3 / dx  - hjmh**3 / dx )
        A[i,i+1] = - beta1/(2.0*dx)*(hjph**3 / dx )
    
        RHS[i] = G[i]
    
    i = n-1   
    A[i,i-1] = 0.0
    A[i,i] = 1.0
    
    RHS[i] = u1
    
    u = solve(A,RHS)
    
    return u



g = 9.81
a0 = 1.0
a1 = 0.7
sx = -25.0
ex = 25.0


expn = 25
L2uO = zeros(expn)
L2uN = zeros(expn)
dxs = zeros(expn)
for j in range(expn):
    n= 2.0**(j/2.0)
    dx = ((ex - sx))/n
    x = arange(sx ,ex+ dx,dx)
    
    h,u,G = Peakon(x,dx,a0) #Soliton(x,0.0,dx,g,a0,a1)
    
    un = FDHGforu(h,G,dx,2.0/3.0,0,1,0,0,1,0)
    un2 =FDHGforuAvg(h,G,dx,2.0/3.0,0,1,0,0,1,0)
    
    L2uO[j] = norm(un -u,ord=2)/norm(u,ord=2)
    L2uN[j] = norm(un2 -u,ord=2)/norm(u,ord=2)
    dxs[j] =dx
    
    
    
