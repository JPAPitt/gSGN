
#Assume we have the expressions at all the necessary points then solved the matrix equation

from numpy import *
from numpy.linalg import solve
from matplotlib.pyplot import plot


def FEMfroGElem(hp,up,dx,j):
    
    hjmh = hp[3*j]
    hjms = hp[3*j+1]
    hjps = hp[3*j+2]
    hjph = hp[3*j+3] 
    
    Ge = zeros((4,4))
    
    
    Ge[0,0] = dx/2*(16.0/105)
    Ge[0,1] = dx/2*(33.0/280)
    Ge[0,2] = dx/2*(-3.0/70)
    Ge[0,3] = dx/2*(19.0/840)
    Ge[1,0] = dx/2*(33.0/280)
    Ge[1,1] = dx/2*(27.0/35)
    Ge[1,2] = dx/2*(-27.0/280)
    Ge[1,3] = dx/2*(-3.0/70)
    Ge[2,0] = dx/2*(-3.0/70)
    Ge[2,1] = dx/2*(-27.0/280)
    Ge[2,2] = dx/2*(27.0/35)
    Ge[2,3] = dx/2*(33.0/280)
    Ge[3,0] = dx/2*(19.0/840)
    Ge[3,1] = dx/2*(-3.0/70)
    Ge[3,2] = dx/2*(33.0/280)
    Ge[3,3] = dx/2*(16.0/105)


    
    uhe = zeros((4,4))
    
    uhe[0,0] = dx/2*(17*hjmh/160 + 9*hjms/140 + hjph/168 - 27*hjps/1120)
    uhe[0,1] = dx/2*(9*hjmh/140 + 27*hjms/280 + 3*hjph/560 - 27*hjps/560)
    uhe[0,2] = dx/2*(-27*hjmh/1120 - 27*hjms/560 + 3*hjph/560 + 27*hjps/1120)
    uhe[0,3] = dx/2*(hjmh/168 + 3*hjms/560 + hjph/168 + 3*hjps/560)
    uhe[1,0] = dx/2*(9*hjmh/140 + 27*hjms/280 + 3*hjph/560 - 27*hjps/560)
    uhe[1,1] = dx/2*(27*hjmh/280 + 729*hjms/1120 + 27*hjph/1120)
    uhe[1,2] = dx/2*(-27*hjmh/560 - 27*hjph/560)
    uhe[1,3] = dx/2*(3*hjmh/560 + 27*hjms/1120 - 27*hjph/1120 - 27*hjps/560)
    uhe[2,0] = dx/2*(-27*hjmh/1120 - 27*hjms/560 + 3*hjph/560 + 27*hjps/1120)
    uhe[2,1] = dx/2*(-27*hjmh/560 - 27*hjph/560)
    uhe[2,2] = dx/2*(27*hjmh/1120 + 27*hjph/280 + 729*hjps/1120)
    uhe[2,3] = dx/2*(3*hjmh/560 - 27*hjms/560 + 9*hjph/140 + 27*hjps/280)
    uhe[3,0] = dx/2*(hjmh/168 + 3*hjms/560 + hjph/168 + 3*hjps/560)
    uhe[3,1] = dx/2*(3*hjmh/560 + 27*hjms/1120 - 27*hjph/1120 - 27*hjps/560)
    uhe[3,2] = dx/2*(3*hjmh/560 - 27*hjms/560 + 9*hjph/140 + 27*hjps/280)
    uhe[3,3] = dx/2*(hjmh/168 - 27*hjms/1120 + 17*hjph/160 + 9*hjps/140)
    

    
    
    h3uxe = zeros((4,4))
    h3uxe[0,0] = (1.0/3.0)*2/dx*(5377*hjmh**3/8960 + 490239*hjmh**2*hjms/640640 + 45223*hjmh**2*hjph/640640 - 7695*hjmh**2*hjps/23296 + 19035*hjmh*hjms**2/23296 + 2601*hjmh*hjms*hjph/20020 - 100521*hjmh*hjms*hjps/160160 + 1163*hjmh*hjph**2/183040 - 8109*hjmh*hjph*hjps/160160 + 14499*hjmh*hjps**2/116480 + 235467*hjms**3/366080 + 297837*hjms**2*hjph/2562560 - 728271*hjms**2*hjps/1281280 + 9963*hjms*hjph**2/2562560 - 1377*hjms*hjph*hjps/16016 + 124659*hjms*hjps**2/640640 + 953*hjph**3/73216 + 8397*hjph**2*hjps/1281280 + 13527*hjph*hjps**2/640640 + 729*hjps**3/1281280)
    h3uxe[0,1] = (1.0/3.0)*2/dx*(-1164861*hjmh**3/1281280 - 2657367*hjmh**2*hjms/2562560 - 251253*hjmh**2*hjph/2562560 + 584091*hjmh**2*hjps/1281280 - 2374353*hjmh*hjms**2/2562560 - 39609*hjmh*hjms*hjph/256256 + 963009*hjmh*hjms*hjps/1281280 - 28737*hjmh*hjph**2/2562560 + 75087*hjmh*hjph*hjps/1281280 - 197559*hjmh*hjps**2/1281280 - 518319*hjms**3/1281280 - 21141*hjms**2*hjph/183040 + 269001*hjms**2*hjps/512512 + 1377*hjms*hjph**2/116480 + 143613*hjms*hjph*hjps/1281280 - 518319*hjms*hjps**2/2562560 - 6939*hjph**3/116480 - 8343*hjph**2*hjps/197120 - 114453*hjph*hjps**2/2562560 - 150903*hjps**3/1281280)
    h3uxe[0,2] = (1.0/3.0)*2/dx*(503469*hjmh**3/1281280 + 108783*hjmh**2*hjms/320320 + 171*hjmh**2*hjph/4928 - 201771*hjmh**2*hjps/1281280 + 165483*hjmh*hjms**2/1281280 + 2187*hjmh*hjms*hjph/80080 - 729*hjmh*hjms*hjps/4928 + 15507*hjmh*hjph**2/1281280 - 243*hjmh*hjph*hjps/45760 + 43011*hjmh*hjps**2/1281280 - 59049*hjms**3/197120 + 729*hjms**2*hjph/232960 + 6561*hjms**2*hjps/116480 - 11097*hjms*hjph**2/232960 - 729*hjms*hjph*hjps/14560 + 6561*hjms*hjps**2/320320 + 47763*hjph**3/366080 + 12069*hjph**2*hjps/116480 + 13851*hjph*hjps**2/320320 + 6561*hjps**3/116480)
    h3uxe[0,3] = (1.0/3.0)*2/dx*(-107519*hjmh**3/1281280 - 173853*hjmh**2*hjms/2562560 - 18559*hjmh**2*hjph/2562560 + 8181*hjmh**2*hjps/256256 - 7209*hjmh*hjms**2/366080 - 3411*hjmh*hjms*hjph/1281280 + 30699*hjmh*hjms*hjps/1281280 - 18559*hjmh*hjph**2/2562560 - 3411*hjmh*hjph*hjps/1281280 - 4941*hjmh*hjps**2/1281280 + 78003*hjms**3/1281280 - 4941*hjms**2*hjph/1281280 - 6561*hjms**2*hjps/512512 + 8181*hjms*hjph**2/256256 + 30699*hjms*hjph*hjps/1281280 - 6561*hjms*hjps**2/512512 - 107519*hjph**3/1281280 - 173853*hjph**2*hjps/2562560 - 7209*hjph*hjps**2/366080 + 78003*hjps**3/1281280)
    h3uxe[1,0] = (1.0/3.0)*2/dx*(-1164861*hjmh**3/1281280 - 2657367*hjmh**2*hjms/2562560 - 251253*hjmh**2*hjph/2562560 + 584091*hjmh**2*hjps/1281280 - 2374353*hjmh*hjms**2/2562560 - 39609*hjmh*hjms*hjph/256256 + 963009*hjmh*hjms*hjps/1281280 - 28737*hjmh*hjph**2/2562560 + 75087*hjmh*hjph*hjps/1281280 - 197559*hjmh*hjps**2/1281280 - 518319*hjms**3/1281280 - 21141*hjms**2*hjph/183040 + 269001*hjms**2*hjps/512512 + 1377*hjms*hjph**2/116480 + 143613*hjms*hjph*hjps/1281280 - 518319*hjms*hjps**2/2562560 - 6939*hjph**3/116480 - 8343*hjph**2*hjps/197120 - 114453*hjph*hjps**2/2562560 - 150903*hjps**3/1281280)
    h3uxe[1,1] = (1.0/3.0)*2/dx*(323217*hjmh**3/232960 + 1858221*hjmh**2*hjms/1281280 + 13689*hjmh**2*hjph/98560 - 1630773*hjmh**2*hjps/2562560 + 19683*hjmh*hjms**2/18304 + 729*hjmh*hjms*hjph/3640 - 6561*hjmh*hjms*hjps/6160 + 18873*hjmh*hjph**2/640640 - 729*hjmh*hjph*hjps/16016 + 400221*hjmh*hjps**2/2562560 + 85293*hjms**3/98560 + 111537*hjms**2*hjph/1281280 - 19683*hjms**2*hjps/640640 - 21141*hjms*hjph**2/256256 - 6561*hjms*hjph*hjps/20020 + 846369*hjms*hjps**2/1281280 + 50409*hjph**3/183040 + 140697*hjph**2*hjps/640640 + 19683*hjph*hjps**2/183040 + 478953*hjps**3/512512)
    h3uxe[1,2] = (1.0/3.0)*2/dx*(-111429*hjmh**3/183040 - 1324593*hjmh**2*hjms/2562560 - 27135*hjmh**2*hjph/512512 + 292329*hjmh**2*hjps/1281280 - 98415*hjmh*hjms**2/512512 - 51759*hjmh*hjms*hjph/1281280 + 465831*hjmh*hjms*hjps/1281280 - 27135*hjmh*hjph**2/512512 - 51759*hjmh*hjph*hjps/1281280 - 6561*hjmh*hjps**2/1281280 - 662661*hjms**3/1281280 - 6561*hjms**2*hjph/1281280 - 1318761*hjms**2*hjps/2562560 + 292329*hjms*hjph**2/1281280 + 465831*hjms*hjph*hjps/1281280 - 1318761*hjms*hjps**2/2562560 - 111429*hjph**3/183040 - 1324593*hjph**2*hjps/2562560 - 98415*hjph*hjps**2/512512 - 662661*hjps**3/1281280)
    h3uxe[1,3] = (1.0/3.0)*2/dx*(47763*hjmh**3/366080 + 12069*hjmh**2*hjms/116480 + 15507*hjmh**2*hjph/1281280 - 11097*hjmh**2*hjps/232960 + 13851*hjmh*hjms**2/320320 - 243*hjmh*hjms*hjph/45760 - 729*hjmh*hjms*hjps/14560 + 171*hjmh*hjph**2/4928 + 2187*hjmh*hjph*hjps/80080 + 729*hjmh*hjps**2/232960 + 6561*hjms**3/116480 + 43011*hjms**2*hjph/1281280 + 6561*hjms**2*hjps/320320 - 201771*hjms*hjph**2/1281280 - 729*hjms*hjph*hjps/4928 + 6561*hjms*hjps**2/116480 + 503469*hjph**3/1281280 + 108783*hjph**2*hjps/320320 + 165483*hjph*hjps**2/1281280 - 59049*hjps**3/197120)
    h3uxe[2,0] = (1.0/3.0)*2/dx*(503469*hjmh**3/1281280 + 108783*hjmh**2*hjms/320320 + 171*hjmh**2*hjph/4928 - 201771*hjmh**2*hjps/1281280 + 165483*hjmh*hjms**2/1281280 + 2187*hjmh*hjms*hjph/80080 - 729*hjmh*hjms*hjps/4928 + 15507*hjmh*hjph**2/1281280 - 243*hjmh*hjph*hjps/45760 + 43011*hjmh*hjps**2/1281280 - 59049*hjms**3/197120 + 729*hjms**2*hjph/232960 + 6561*hjms**2*hjps/116480 - 11097*hjms*hjph**2/232960 - 729*hjms*hjph*hjps/14560 + 6561*hjms*hjps**2/320320 + 47763*hjph**3/366080 + 12069*hjph**2*hjps/116480 + 13851*hjph*hjps**2/320320 + 6561*hjps**3/116480)
    h3uxe[2,1] = (1.0/3.0)*2/dx*(-111429*hjmh**3/183040 - 1324593*hjmh**2*hjms/2562560 - 27135*hjmh**2*hjph/512512 + 292329*hjmh**2*hjps/1281280 - 98415*hjmh*hjms**2/512512 - 51759*hjmh*hjms*hjph/1281280 + 465831*hjmh*hjms*hjps/1281280 - 27135*hjmh*hjph**2/512512 - 51759*hjmh*hjph*hjps/1281280 - 6561*hjmh*hjps**2/1281280 - 662661*hjms**3/1281280 - 6561*hjms**2*hjph/1281280 - 1318761*hjms**2*hjps/2562560 + 292329*hjms*hjph**2/1281280 + 465831*hjms*hjph*hjps/1281280 - 1318761*hjms*hjps**2/2562560 - 111429*hjph**3/183040 - 1324593*hjph**2*hjps/2562560 - 98415*hjph*hjps**2/512512 - 662661*hjps**3/1281280)
    h3uxe[2,2] = (1.0/3.0)*2/dx*(50409*hjmh**3/183040 + 140697*hjmh**2*hjms/640640 + 18873*hjmh**2*hjph/640640 - 21141*hjmh**2*hjps/256256 + 19683*hjmh*hjms**2/183040 - 729*hjmh*hjms*hjph/16016 - 6561*hjmh*hjms*hjps/20020 + 13689*hjmh*hjph**2/98560 + 729*hjmh*hjph*hjps/3640 + 111537*hjmh*hjps**2/1281280 + 478953*hjms**3/512512 + 400221*hjms**2*hjph/2562560 + 846369*hjms**2*hjps/1281280 - 1630773*hjms*hjph**2/2562560 - 6561*hjms*hjph*hjps/6160 - 19683*hjms*hjps**2/640640 + 323217*hjph**3/232960 + 1858221*hjph**2*hjps/1281280 + 19683*hjph*hjps**2/18304 + 85293*hjps**3/98560)
    h3uxe[2,3] = (1.0/3.0)*2/dx*(-6939*hjmh**3/116480 - 8343*hjmh**2*hjms/197120 - 28737*hjmh**2*hjph/2562560 + 1377*hjmh**2*hjps/116480 - 114453*hjmh*hjms**2/2562560 + 75087*hjmh*hjms*hjph/1281280 + 143613*hjmh*hjms*hjps/1281280 - 251253*hjmh*hjph**2/2562560 - 39609*hjmh*hjph*hjps/256256 - 21141*hjmh*hjps**2/183040 - 150903*hjms**3/1281280 - 197559*hjms**2*hjph/1281280 - 518319*hjms**2*hjps/2562560 + 584091*hjms*hjph**2/1281280 + 963009*hjms*hjph*hjps/1281280 + 269001*hjms*hjps**2/512512 - 1164861*hjph**3/1281280 - 2657367*hjph**2*hjps/2562560 - 2374353*hjph*hjps**2/2562560 - 518319*hjps**3/1281280)
    h3uxe[3,0] = (1.0/3.0)*2/dx*(-107519*hjmh**3/1281280 - 173853*hjmh**2*hjms/2562560 - 18559*hjmh**2*hjph/2562560 + 8181*hjmh**2*hjps/256256 - 7209*hjmh*hjms**2/366080 - 3411*hjmh*hjms*hjph/1281280 + 30699*hjmh*hjms*hjps/1281280 - 18559*hjmh*hjph**2/2562560 - 3411*hjmh*hjph*hjps/1281280 - 4941*hjmh*hjps**2/1281280 + 78003*hjms**3/1281280 - 4941*hjms**2*hjph/1281280 - 6561*hjms**2*hjps/512512 + 8181*hjms*hjph**2/256256 + 30699*hjms*hjph*hjps/1281280 - 6561*hjms*hjps**2/512512 - 107519*hjph**3/1281280 - 173853*hjph**2*hjps/2562560 - 7209*hjph*hjps**2/366080 + 78003*hjps**3/1281280)
    h3uxe[3,1] = (1.0/3.0)*2/dx*(47763*hjmh**3/366080 + 12069*hjmh**2*hjms/116480 + 15507*hjmh**2*hjph/1281280 - 11097*hjmh**2*hjps/232960 + 13851*hjmh*hjms**2/320320 - 243*hjmh*hjms*hjph/45760 - 729*hjmh*hjms*hjps/14560 + 171*hjmh*hjph**2/4928 + 2187*hjmh*hjph*hjps/80080 + 729*hjmh*hjps**2/232960 + 6561*hjms**3/116480 + 43011*hjms**2*hjph/1281280 + 6561*hjms**2*hjps/320320 - 201771*hjms*hjph**2/1281280 - 729*hjms*hjph*hjps/4928 + 6561*hjms*hjps**2/116480 + 503469*hjph**3/1281280 + 108783*hjph**2*hjps/320320 + 165483*hjph*hjps**2/1281280 - 59049*hjps**3/197120)
    h3uxe[3,2] = (1.0/3.0)*2/dx*(-6939*hjmh**3/116480 - 8343*hjmh**2*hjms/197120 - 28737*hjmh**2*hjph/2562560 + 1377*hjmh**2*hjps/116480 - 114453*hjmh*hjms**2/2562560 + 75087*hjmh*hjms*hjph/1281280 + 143613*hjmh*hjms*hjps/1281280 - 251253*hjmh*hjph**2/2562560 - 39609*hjmh*hjph*hjps/256256 - 21141*hjmh*hjps**2/183040 - 150903*hjms**3/1281280 - 197559*hjms**2*hjph/1281280 - 518319*hjms**2*hjps/2562560 + 584091*hjms*hjph**2/1281280 + 963009*hjms*hjph*hjps/1281280 + 269001*hjms*hjps**2/512512 - 1164861*hjph**3/1281280 - 2657367*hjph**2*hjps/2562560 - 2374353*hjph*hjps**2/2562560 - 518319*hjps**3/1281280)
    h3uxe[3,3] = (1.0/3.0)*2/dx*(953*hjmh**3/73216 + 8397*hjmh**2*hjms/1281280 + 1163*hjmh**2*hjph/183040 + 9963*hjmh**2*hjps/2562560 + 13527*hjmh*hjms**2/640640 - 8109*hjmh*hjms*hjph/160160 - 1377*hjmh*hjms*hjps/16016 + 45223*hjmh*hjph**2/640640 + 2601*hjmh*hjph*hjps/20020 + 297837*hjmh*hjps**2/2562560 + 729*hjms**3/1281280 + 14499*hjms**2*hjph/116480 + 124659*hjms**2*hjps/640640 - 7695*hjms*hjph**2/23296 - 100521*hjms*hjph*hjps/160160 - 728271*hjms*hjps**2/1281280 + 5377*hjph**3/8960 + 490239*hjph**2*hjps/640640 + 19035*hjph*hjps**2/23296 + 235467*hjps**3/366080)   

    
    # print(Ge)
    # print(uhe)
    # print()
    return Ge, uhe,h3uxe

def FEMforG(xn, hp,up,dx):
    
    n = len(xn)
    A = zeros((3*n+1,3*n+1))
    b = zeros((3*n+1,3*n+1))
            
    
    for i in range(n):
        #elementwisematrices
        #Ge, uhe, h3uxe = FEMfroGElem(hp,up,dx,i)
        
        Ge, uhe,h3uxe = FEMfroGElem(hp,up,dx,i)
        
        A[3*i : 3*i + 4, 3*i : 3*i + 4] =  A[3*i : 3*i + 4, 3*i : 3*i + 4] + Ge
        b[3*i : 3*i + 4, 3*i : 3*i + 4] =  b[3*i : 3*i + 4, 3*i : 3*i + 4] +  uhe + h3uxe

    A[0,0:4] = [1,0,0,0]
    b[0,0:4] = [0,0,0,0]
    
    A[-1,-4:] = [0,0,0,1]
    b[-1,-4:] = [0,0,0,0]
    

    return A,b


def initialhu(x,dx,a1):
    n = len(x)
    hp = []
    up = []
    xp = []

    c = sqrt(1 + a1)
    
    #i = 0
    i=0
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
    
    #u   
    ujmh = c*(hjmh - 1)/hjmh 
    ujms = c*(hjms - 1)/hjms 
    ujps = c*(hjps - 1)/hjps 
    ujph = c*(hjph - 1)/hjph 
    
    ucell = [ujmh,ujms,ujps,ujph]
    up = up + ucell
    
    for i in range(1,n):
        xjms = x[i] - dx/6
        xj = x[i]
        xjps = x[i] + dx/6
        xjph = x[i] + dx/2
        
        #h
        hjms = 1 + a1*exp(-sqrt(3)*abs(xjms))
        hjps = 1 + a1*exp(-sqrt(3)*abs(xjps))
        hjph = 1 + a1*exp(-sqrt(3)*abs(xjph))
        hcell = [hjms,hjps,hjph]
        xhcell = [xjms,xjps,xjph]
        hp = hp + hcell
        xp = xp + xhcell
        
        #u      
        ujms = c*(hjms - 1)/hjms 
        ujps = c*(hjps - 1)/hjps 
        ujph = c*(hjph - 1)/hjph 
        
        ucell = [ujms,ujps,ujph]
        up = up + ucell
        
        
    return hp,up,xp

def GetCellMidpoint(xn,q,dx):
    
    n = len(xn)
    qn = zeros(n)
    # cubic passes through y1, y2, y3, y4
    # pa = -9 Gjmh + 27 Gms - 27 Gps + 9 Gph / 2*dx3
    # pb = 9 Gjmh - 9 Gms - 9 Gps + 9 Gph / 4*dx2
    # pc = Gjmh - 27 Gms + 27 Gps - Gph / 8*dx
    # pd = -Gjmh + 9 Gms + 9 Gps - Gph / 16
    qn[0]=q[0]
    for i in range(1,n-1):
        qmh = q[3*i]
        qms = q[3*i + 1]
        qps = q[3*i + 2]
        qph = q[3*i + 3]
        qn[i] = (-qmh + 9*qms + 9*qps - qph) / 16.0
    
    qn[-1]=q[-1]
    return qn

n= 1000
dx = 10.0/n
xn = arange(-5,5+ dx,dx)
a1 = 0.5


h,u,x = initialhu(xn,dx,a1)

Guh = array(h)*array(u)

A,b = FEMforG(xn,h,u,dx)

Gc = solve(A, b.dot(u))

Gm = GetCellMidpoint(xn,Gc,dx)


