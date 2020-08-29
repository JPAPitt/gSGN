"""
Python code to generate csv files that contain integrals of the basis functions
needed to generate the matrices of the Finite Element Method.

by Jordan Pitt 11/10/2018
"""

"""
############################# IMPORTS #########################################
"""


from numpy import linspace,sqrt,exp,sign
from matplotlib.pyplot import plot,show,title

a0 = 1
a1 = 0.5
g = 1
c = sqrt(g*(1 + a1/a0))
t = 0
x = linspace(-25,25,1000)
dx = 50.0 /1000.0
h = len(x)*[0]
u = len(x)*[0]
nx = len(x)*[0]
ux =len(x)*[0]
uxx =len(x)*[0]

anx = len(x)*[0]
anxx = len(x)*[0]
aux =len(x)*[0]
auxx =len(x)*[0]
G =len(x)*[0]
aG =len(x)*[0]
Gc = len(x)*[0]
for i in range(1,len(x)-1):
    nim1 = a1*exp(-sqrt(3)*abs(x[i-1]) / a0)
    nip1 = a1*exp(-sqrt(3)*abs(x[i+1]) / a0)
    
    ni = a1*exp(-sqrt(3)*abs(x[i]) / a0)
    h[i] = a0 +  ni
    u[i] = c*(ni / (1+ni))
    
    nx[i] = sqrt(3)*a1/a0*sign(-x[i]) *exp(-sqrt(3)*abs(x[i]) / a0)
    nxx =  -( a1 /a0**2) *(2*sqrt(3)*a0*int(-x[i] == 0) - 3*sign(-x)**2 ) *exp(-sqrt(3)*abs(x[i]) / a0)
    ux[i] = c*nx[i]/(ni + 1)**2
    uxx[i] = c*((1 + ni)*nxx - 2*nx[i]**2)/(ni + 1)**3
    G[i] = u[i]*h[i] - (h[i]**3)*uxx[i]/3 - h[i]**2*nx[i]*ux[i]
    
    Gc[i] = u[i]*h[i] + c/3*(6*a1**2*(a1*exp(-sqrt(3)*abs(c*t - x[i])) + 1)*exp(-2*sqrt(3)*abs(c*t - x[i]))*sign(c*t - x[i])**2 \
    - 2*sqrt(3)*a1*(a1*exp(-sqrt(3)*abs(c*t - x[i])) + 1)**2*exp(-sqrt(3)*abs(c*t - x[i]))*int(c*t - x[i] == 0) \
    + 3*a1*(a1*exp(-sqrt(3)*abs(c*t - x[i])) + 1)**2*exp(-sqrt(3)*abs(c*t - x[i]))*sign(c*t - x[i])**2)
    
    anx[i] = (nip1 - nim1) / 2*dx
    anxx[i] = (nip1- 2*ni+ nim1) / dx**2
    aux[i] = c*anx[i]/(ni + 1)**2
    auxx[i] = c*((1 + ni)*anxx[i] - 2*anx[i]**2)/(ni + 1)**3
    aG[i] = u[i]*h[i] - (h[i]**3)*auxx[i]/3 - h[i]**2*anx[i]*aux[i]

    

plot(x,h)
title('h')
show()
plot(x,G)
title('G')
show()
plot(x,Gc)
title('Gc')
show()
plot(x,aG)
title('aG')
show()