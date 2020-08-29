
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
    
    
    Ge[0,0] = dx/2.0*(16.0/105)
    Ge[0,1] = dx/2.0*(33.0/280 )
    Ge[0,2] = dx/2.0*(- 3.0/70)
    Ge[0,3] = dx/2.0*(19.0/840)
    
    Ge[1,0] = dx/2.0*(33.0/280)
    Ge[1,1] = dx/2.0*(27.0/35)
    Ge[1,2] = dx/2.0*(-27.0/280 )
    Ge[1,3] = dx/2.0*(-3.0/70)

    Ge[2,0] = dx/2.0*(-3.0/70)
    Ge[2,1] = dx/2.0*(- 27.0/280 )
    Ge[2,2] = dx/2.0*(27.0/35 )
    Ge[2,3] = dx/2.0*(33.0/280)    

    Ge[3,0] = dx/2.0*(19.0/840 )
    Ge[3,1] = dx/2.0*(-3.0/70)
    Ge[3,2] = dx/2.0*(33.0/280)
    Ge[3,3] = dx/2.0*(16.0/105)    


    
    uhe = zeros((4,4))
    
    uhe[0,0] = dx/2.0*((17.0/160)*hjmh+(9.0/140)*hjms+(-27.0/1120)*hjps+(1.0/168)*hjph)
    uhe[0,1] = dx/2.0*((9.0/140)*hjmh+(27.0/280)*hjms+(-27.0/560)*hjps+(3.0/560)*hjph)
    uhe[0,2] = dx/2.0*((-27.0/1120)*hjmh+(-27.0/560)*hjms+(27.0/1120)*hjps+(3.0/560)*hjph)
    uhe[0,3] = dx/2.0*((1.0/168)*hjmh+(3.0/560)*hjms+(3.0/560)*hjps+(1.0/168)*hjph)
    
    uhe[1,0] = dx/2.0*((9.0/140)*hjmh+(27.0/280)*hjms+(-27.0/560)*hjps+(3.0/560)*hjph)
    uhe[1,1] = dx/2.0*((27.0/280)*hjmh+(729.0/1120)*hjms+(729.0/1120)*hjps+(27.0/1120)*hjph)
    uhe[1,2] = dx/2.0*((-27.0/560)*hjmh+(-27.0/560)*hjms+(-27.0/560)*hjps+(-27.0/560)*hjph)
    uhe[1,3] = dx/2.0*((3.0/560)*hjmh+(27.0/1120)*hjms+(-27.0/560)*hjps+(-27.0/1120)*hjph)
    
    uhe[2,0] = dx/2.0*((-27.0/1120)*hjmh+(-27.0/560)*hjms+(27.0/1120)*hjps+(3.0/560)*hjph)
    uhe[2,1] = dx/2.0*((-27.0/560)*hjmh+(-27.0/560)*hjms+(-27.0/560)*hjps+(-27.0/560)*hjph)
    uhe[2,2] = dx/2.0*((27.0/1120)*hjmh+(27.0/1120)*hjms+(729.0/1120)*hjps+(27.0/280)*hjph)
    uhe[2,3] = dx/2.0*((3.0/560)*hjmh+(-27.0/560)*hjms+(27.0/280)*hjps+(9.0/140)*hjph)
    
    uhe[3,0] = dx/2.0*((1.0/168)*hjmh+(3.0/560)*hjms+(3.0/560)*hjps+(1.0/168)*hjph)
    uhe[3,1] = dx/2.0*((3.0/560)*hjmh+(27.0/1120)*hjms+(-27.0/560)*hjps+(-27.0/1120)*hjph)
    uhe[3,2] = dx/2.0*((3.0/560)*hjmh+(-27.0/560)*hjms+(27.0/280)*hjps+(9.0/140)*hjph)
    uhe[3,3] = dx/2.0*((1.0/168)*hjmh+(-27.0/1120)*hjms+(9.0/140)*hjps+(17.0/160)*hjph)
    

    
    
    # h3uxe = zeros((4,4))
    # h3uxe[0,0] = 2.0/dx*((5377.0/8960)*hjmh*hjmh*hjmh+(163413.0/640640)*hjms*hjmh*hjmh+(-2565.0/23296)*hjps*hjmh*hjmh+(45223.0/1921920)*hjph*hjmh*hjmh+(163413.0/640640)*hjmh*hjms*hjmh+(6345.0/23296)*hjms*hjms*hjmh+(-33507.0/320320)*hjps*hjms*hjmh+(867.0/40040)*hjph*hjms*hjmh+(-2565.0/23296)*hjmh*hjps*hjmh+(-33507.0/320320)*hjms*hjps*hjmh+(4833.0/116480)*hjps*hjps*hjmh+(-2703.0/320320)*hjph*hjps*hjmh+(45223.0/1921920)*hjmh*hjph*hjmh+(867.0/40040)*hjms*hjph*hjmh+(-2703.0/320320)*hjps*hjph*hjmh+(1163.0/549120)*hjph*hjph*hjmh+(163413.0/640640)*hjmh*hjmh*hjms+(6345.0/23296)*hjms*hjmh*hjms+(-33507.0/320320)*hjps*hjmh*hjms+(867.0/40040)*hjph*hjmh*hjms+(6345.0/23296)*hjmh*hjms*hjms+(235467.0/366080)*hjms*hjms*hjms+(-242757.0/1281280)*hjps*hjms*hjms+(99279.0/2562560)*hjph*hjms*hjms+(-33507.0/320320)*hjmh*hjps*hjms+(-242757.0/1281280)*hjms*hjps*hjms+(41553.0/640640)*hjps*hjps*hjms+(-459.0/32032)*hjph*hjps*hjms+(867.0/40040)*hjmh*hjph*hjms+(99279.0/2562560)*hjms*hjph*hjms+(-459.0/32032)*hjps*hjph*hjms+(3321.0/2562560)*hjph*hjph*hjms+(-2565.0/23296)*hjmh*hjmh*hjps+(-33507.0/320320)*hjms*hjmh*hjps+(4833.0/116480)*hjps*hjmh*hjps+(-2703.0/320320)*hjph*hjmh*hjps+(-33507.0/320320)*hjmh*hjms*hjps+(-242757.0/1281280)*hjms*hjms*hjps+(41553.0/640640)*hjps*hjms*hjps+(-459.0/32032)*hjph*hjms*hjps+(4833.0/116480)*hjmh*hjps*hjps+(41553.0/640640)*hjms*hjps*hjps+(729.0/1281280)*hjps*hjps*hjps+(4509.0/640640)*hjph*hjps*hjps+(-2703.0/320320)*hjmh*hjph*hjps+(-459.0/32032)*hjms*hjph*hjps+(4509.0/640640)*hjps*hjph*hjps+(2799.0/1281280)*hjph*hjph*hjps+(45223.0/1921920)*hjmh*hjmh*hjph+(867.0/40040)*hjms*hjmh*hjph+(-2703.0/320320)*hjps*hjmh*hjph+(1163.0/549120)*hjph*hjmh*hjph+(867.0/40040)*hjmh*hjms*hjph+(99279.0/2562560)*hjms*hjms*hjph+(-459.0/32032)*hjps*hjms*hjph+(3321.0/2562560)*hjph*hjms*hjph+(-2703.0/320320)*hjmh*hjps*hjph+(-459.0/32032)*hjms*hjps*hjph+(4509.0/640640)*hjps*hjps*hjph+(2799.0/1281280)*hjph*hjps*hjph+(1163.0/549120)*hjmh*hjph*hjph+(3321.0/2562560)*hjms*hjph*hjph+(2799.0/1281280)*hjps*hjph*hjph+(953.0/73216)*hjph*hjph*hjph)
    # h3uxe[0,1] = 2.0/dx*((-1164861.0/1281280)*hjmh*hjmh*hjmh+(-885789.0/2562560)*hjms*hjmh*hjmh+(194697.0/1281280)*hjps*hjmh*hjmh+(-83751.0/2562560)*hjph*hjmh*hjmh+(-885789.0/2562560)*hjmh*hjms*hjmh+(-791451.0/2562560)*hjms*hjms*hjmh+(321003.0/2562560)*hjps*hjms*hjmh+(-13203.0/512512)*hjph*hjms*hjmh+(194697.0/1281280)*hjmh*hjps*hjmh+(321003.0/2562560)*hjms*hjps*hjmh+(-65853.0/1281280)*hjps*hjps*hjmh+(25029.0/2562560)*hjph*hjps*hjmh+(-83751.0/2562560)*hjmh*hjph*hjmh+(-13203.0/512512)*hjms*hjph*hjmh+(25029.0/2562560)*hjps*hjph*hjmh+(-9579.0/2562560)*hjph*hjph*hjmh+(-885789.0/2562560)*hjmh*hjmh*hjms+(-791451.0/2562560)*hjms*hjmh*hjms+(321003.0/2562560)*hjps*hjmh*hjms+(-13203.0/512512)*hjph*hjmh*hjms+(-791451.0/2562560)*hjmh*hjms*hjms+(-518319.0/1281280)*hjms*hjms*hjms+(89667.0/512512)*hjps*hjms*hjms+(-7047.0/183040)*hjph*hjms*hjms+(321003.0/2562560)*hjmh*hjps*hjms+(89667.0/512512)*hjms*hjps*hjms+(-172773.0/2562560)*hjps*hjps*hjms+(47871.0/2562560)*hjph*hjps*hjms+(-13203.0/512512)*hjmh*hjph*hjms+(-7047.0/183040)*hjms*hjph*hjms+(47871.0/2562560)*hjps*hjph*hjms+(459.0/116480)*hjph*hjph*hjms+(194697.0/1281280)*hjmh*hjmh*hjps+(321003.0/2562560)*hjms*hjmh*hjps+(-65853.0/1281280)*hjps*hjmh*hjps+(25029.0/2562560)*hjph*hjmh*hjps+(321003.0/2562560)*hjmh*hjms*hjps+(89667.0/512512)*hjms*hjms*hjps+(-172773.0/2562560)*hjps*hjms*hjps+(47871.0/2562560)*hjph*hjms*hjps+(-65853.0/1281280)*hjmh*hjps*hjps+(-172773.0/2562560)*hjms*hjps*hjps+(-150903.0/1281280)*hjps*hjps*hjps+(-38151.0/2562560)*hjph*hjps*hjps+(25029.0/2562560)*hjmh*hjph*hjps+(47871.0/2562560)*hjms*hjph*hjps+(-38151.0/2562560)*hjps*hjph*hjps+(-2781.0/197120)*hjph*hjph*hjps+(-83751.0/2562560)*hjmh*hjmh*hjph+(-13203.0/512512)*hjms*hjmh*hjph+(25029.0/2562560)*hjps*hjmh*hjph+(-9579.0/2562560)*hjph*hjmh*hjph+(-13203.0/512512)*hjmh*hjms*hjph+(-7047.0/183040)*hjms*hjms*hjph+(47871.0/2562560)*hjps*hjms*hjph+(459.0/116480)*hjph*hjms*hjph+(25029.0/2562560)*hjmh*hjps*hjph+(47871.0/2562560)*hjms*hjps*hjph+(-38151.0/2562560)*hjps*hjps*hjph+(-2781.0/197120)*hjph*hjps*hjph+(-9579.0/2562560)*hjmh*hjph*hjph+(459.0/116480)*hjms*hjph*hjph+(-2781.0/197120)*hjps*hjph*hjph+(-6939.0/116480)*hjph*hjph*hjph)
    # h3uxe[0,2] = 2.0/dx*((503469.0/1281280)*hjmh*hjmh*hjmh+(36261.0/320320)*hjms*hjmh*hjmh+(-67257.0/1281280)*hjps*hjmh*hjmh+(57.0/4928)*hjph*hjmh*hjmh+(36261.0/320320)*hjmh*hjms*hjmh+(55161.0/1281280)*hjms*hjms*hjmh+(-243.0/9856)*hjps*hjms*hjmh+(729.0/160160)*hjph*hjms*hjmh+(-67257.0/1281280)*hjmh*hjps*hjmh+(-243.0/9856)*hjms*hjps*hjmh+(14337.0/1281280)*hjps*hjps*hjmh+(-81.0/91520)*hjph*hjps*hjmh+(57.0/4928)*hjmh*hjph*hjmh+(729.0/160160)*hjms*hjph*hjmh+(-81.0/91520)*hjps*hjph*hjmh+(5169.0/1281280)*hjph*hjph*hjmh+(36261.0/320320)*hjmh*hjmh*hjms+(55161.0/1281280)*hjms*hjmh*hjms+(-243.0/9856)*hjps*hjmh*hjms+(729.0/160160)*hjph*hjmh*hjms+(55161.0/1281280)*hjmh*hjms*hjms+(-59049.0/197120)*hjms*hjms*hjms+(2187.0/116480)*hjps*hjms*hjms+(243.0/232960)*hjph*hjms*hjms+(-243.0/9856)*hjmh*hjps*hjms+(2187.0/116480)*hjms*hjps*hjms+(2187.0/320320)*hjps*hjps*hjms+(-243.0/29120)*hjph*hjps*hjms+(729.0/160160)*hjmh*hjph*hjms+(243.0/232960)*hjms*hjph*hjms+(-243.0/29120)*hjps*hjph*hjms+(-3699.0/232960)*hjph*hjph*hjms+(-67257.0/1281280)*hjmh*hjmh*hjps+(-243.0/9856)*hjms*hjmh*hjps+(14337.0/1281280)*hjps*hjmh*hjps+(-81.0/91520)*hjph*hjmh*hjps+(-243.0/9856)*hjmh*hjms*hjps+(2187.0/116480)*hjms*hjms*hjps+(2187.0/320320)*hjps*hjms*hjps+(-243.0/29120)*hjph*hjms*hjps+(14337.0/1281280)*hjmh*hjps*hjps+(2187.0/320320)*hjms*hjps*hjps+(6561.0/116480)*hjps*hjps*hjps+(4617.0/320320)*hjph*hjps*hjps+(-81.0/91520)*hjmh*hjph*hjps+(-243.0/29120)*hjms*hjph*hjps+(4617.0/320320)*hjps*hjph*hjps+(4023.0/116480)*hjph*hjph*hjps+(57.0/4928)*hjmh*hjmh*hjph+(729.0/160160)*hjms*hjmh*hjph+(-81.0/91520)*hjps*hjmh*hjph+(5169.0/1281280)*hjph*hjmh*hjph+(729.0/160160)*hjmh*hjms*hjph+(243.0/232960)*hjms*hjms*hjph+(-243.0/29120)*hjps*hjms*hjph+(-3699.0/232960)*hjph*hjms*hjph+(-81.0/91520)*hjmh*hjps*hjph+(-243.0/29120)*hjms*hjps*hjph+(4617.0/320320)*hjps*hjps*hjph+(4023.0/116480)*hjph*hjps*hjph+(5169.0/1281280)*hjmh*hjph*hjph+(-3699.0/232960)*hjms*hjph*hjph+(4023.0/116480)*hjps*hjph*hjph+(47763.0/366080)*hjph*hjph*hjph)
    # h3uxe[0,3] = 2.0/dx*((-107519.0/1281280)*hjmh*hjmh*hjmh+(-57951.0/2562560)*hjms*hjmh*hjmh+(2727.0/256256)*hjps*hjmh*hjmh+(-18559.0/7687680)*hjph*hjmh*hjmh+(-57951.0/2562560)*hjmh*hjms*hjmh+(-2403.0/366080)*hjms*hjms*hjmh+(10233.0/2562560)*hjps*hjms*hjmh+(-1137.0/2562560)*hjph*hjms*hjmh+(2727.0/256256)*hjmh*hjps*hjmh+(10233.0/2562560)*hjms*hjps*hjmh+(-1647.0/1281280)*hjps*hjps*hjmh+(-1137.0/2562560)*hjph*hjps*hjmh+(-18559.0/7687680)*hjmh*hjph*hjmh+(-1137.0/2562560)*hjms*hjph*hjmh+(-1137.0/2562560)*hjps*hjph*hjmh+(-18559.0/7687680)*hjph*hjph*hjmh+(-57951.0/2562560)*hjmh*hjmh*hjms+(-2403.0/366080)*hjms*hjmh*hjms+(10233.0/2562560)*hjps*hjmh*hjms+(-1137.0/2562560)*hjph*hjmh*hjms+(-2403.0/366080)*hjmh*hjms*hjms+(78003.0/1281280)*hjms*hjms*hjms+(-2187.0/512512)*hjps*hjms*hjms+(-1647.0/1281280)*hjph*hjms*hjms+(10233.0/2562560)*hjmh*hjps*hjms+(-2187.0/512512)*hjms*hjps*hjms+(-2187.0/512512)*hjps*hjps*hjms+(10233.0/2562560)*hjph*hjps*hjms+(-1137.0/2562560)*hjmh*hjph*hjms+(-1647.0/1281280)*hjms*hjph*hjms+(10233.0/2562560)*hjps*hjph*hjms+(2727.0/256256)*hjph*hjph*hjms+(2727.0/256256)*hjmh*hjmh*hjps+(10233.0/2562560)*hjms*hjmh*hjps+(-1647.0/1281280)*hjps*hjmh*hjps+(-1137.0/2562560)*hjph*hjmh*hjps+(10233.0/2562560)*hjmh*hjms*hjps+(-2187.0/512512)*hjms*hjms*hjps+(-2187.0/512512)*hjps*hjms*hjps+(10233.0/2562560)*hjph*hjms*hjps+(-1647.0/1281280)*hjmh*hjps*hjps+(-2187.0/512512)*hjms*hjps*hjps+(78003.0/1281280)*hjps*hjps*hjps+(-2403.0/366080)*hjph*hjps*hjps+(-1137.0/2562560)*hjmh*hjph*hjps+(10233.0/2562560)*hjms*hjph*hjps+(-2403.0/366080)*hjps*hjph*hjps+(-57951.0/2562560)*hjph*hjph*hjps+(-18559.0/7687680)*hjmh*hjmh*hjph+(-1137.0/2562560)*hjms*hjmh*hjph+(-1137.0/2562560)*hjps*hjmh*hjph+(-18559.0/7687680)*hjph*hjmh*hjph+(-1137.0/2562560)*hjmh*hjms*hjph+(-1647.0/1281280)*hjms*hjms*hjph+(10233.0/2562560)*hjps*hjms*hjph+(2727.0/256256)*hjph*hjms*hjph+(-1137.0/2562560)*hjmh*hjps*hjph+(10233.0/2562560)*hjms*hjps*hjph+(-2403.0/366080)*hjps*hjps*hjph+(-57951.0/2562560)*hjph*hjps*hjph+(-18559.0/7687680)*hjmh*hjph*hjph+(2727.0/256256)*hjms*hjph*hjph+(-57951.0/2562560)*hjps*hjph*hjph+(-107519.0/1281280)*hjph*hjph*hjph)

    # h3uxe[1,0] = 2.0/dx*((-1164861.0/1281280)*hjmh*hjmh*hjmh+(-885789.0/2562560)*hjms*hjmh*hjmh+(194697.0/1281280)*hjps*hjmh*hjmh+(-83751.0/2562560)*hjph*hjmh*hjmh+(-885789.0/2562560)*hjmh*hjms*hjmh+(-791451.0/2562560)*hjms*hjms*hjmh+(321003.0/2562560)*hjps*hjms*hjmh+(-13203.0/512512)*hjph*hjms*hjmh+(194697.0/1281280)*hjmh*hjps*hjmh+(321003.0/2562560)*hjms*hjps*hjmh+(-65853.0/1281280)*hjps*hjps*hjmh+(25029.0/2562560)*hjph*hjps*hjmh+(-83751.0/2562560)*hjmh*hjph*hjmh+(-13203.0/512512)*hjms*hjph*hjmh+(25029.0/2562560)*hjps*hjph*hjmh+(-9579.0/2562560)*hjph*hjph*hjmh+(-885789.0/2562560)*hjmh*hjmh*hjms+(-791451.0/2562560)*hjms*hjmh*hjms+(321003.0/2562560)*hjps*hjmh*hjms+(-13203.0/512512)*hjph*hjmh*hjms+(-791451.0/2562560)*hjmh*hjms*hjms+(-518319.0/1281280)*hjms*hjms*hjms+(89667.0/512512)*hjps*hjms*hjms+(-7047.0/183040)*hjph*hjms*hjms+(321003.0/2562560)*hjmh*hjps*hjms+(89667.0/512512)*hjms*hjps*hjms+(-172773.0/2562560)*hjps*hjps*hjms+(47871.0/2562560)*hjph*hjps*hjms+(-13203.0/512512)*hjmh*hjph*hjms+(-7047.0/183040)*hjms*hjph*hjms+(47871.0/2562560)*hjps*hjph*hjms+(459.0/116480)*hjph*hjph*hjms+(194697.0/1281280)*hjmh*hjmh*hjps+(321003.0/2562560)*hjms*hjmh*hjps+(-65853.0/1281280)*hjps*hjmh*hjps+(25029.0/2562560)*hjph*hjmh*hjps+(321003.0/2562560)*hjmh*hjms*hjps+(89667.0/512512)*hjms*hjms*hjps+(-172773.0/2562560)*hjps*hjms*hjps+(47871.0/2562560)*hjph*hjms*hjps+(-65853.0/1281280)*hjmh*hjps*hjps+(-172773.0/2562560)*hjms*hjps*hjps+(-150903.0/1281280)*hjps*hjps*hjps+(-38151.0/2562560)*hjph*hjps*hjps+(25029.0/2562560)*hjmh*hjph*hjps+(47871.0/2562560)*hjms*hjph*hjps+(-38151.0/2562560)*hjps*hjph*hjps+(-2781.0/197120)*hjph*hjph*hjps+(-83751.0/2562560)*hjmh*hjmh*hjph+(-13203.0/512512)*hjms*hjmh*hjph+(25029.0/2562560)*hjps*hjmh*hjph+(-9579.0/2562560)*hjph*hjmh*hjph+(-13203.0/512512)*hjmh*hjms*hjph+(-7047.0/183040)*hjms*hjms*hjph+(47871.0/2562560)*hjps*hjms*hjph+(459.0/116480)*hjph*hjms*hjph+(25029.0/2562560)*hjmh*hjps*hjph+(47871.0/2562560)*hjms*hjps*hjph+(-38151.0/2562560)*hjps*hjps*hjph+(-2781.0/197120)*hjph*hjps*hjph+(-9579.0/2562560)*hjmh*hjph*hjph+(459.0/116480)*hjms*hjph*hjph+(-2781.0/197120)*hjps*hjph*hjph+(-6939.0/116480)*hjph*hjph*hjph)
    # h3uxe[1,1] = 2.0/dx*((323217.0/232960)*hjmh*hjmh*hjmh+(619407.0/1281280)*hjms*hjmh*hjmh+(-543591.0/2562560)*hjps*hjmh*hjmh+(4563.0/98560)*hjph*hjmh*hjmh+(619407.0/1281280)*hjmh*hjms*hjmh+(6561.0/18304)*hjms*hjms*hjmh+(-2187.0/12320)*hjps*hjms*hjmh+(243.0/7280)*hjph*hjms*hjmh+(-543591.0/2562560)*hjmh*hjps*hjmh+(-2187.0/12320)*hjms*hjps*hjmh+(133407.0/2562560)*hjps*hjps*hjmh+(-243.0/32032)*hjph*hjps*hjmh+(4563.0/98560)*hjmh*hjph*hjmh+(243.0/7280)*hjms*hjph*hjmh+(-243.0/32032)*hjps*hjph*hjmh+(6291.0/640640)*hjph*hjph*hjmh+(619407.0/1281280)*hjmh*hjmh*hjms+(6561.0/18304)*hjms*hjmh*hjms+(-2187.0/12320)*hjps*hjmh*hjms+(243.0/7280)*hjph*hjmh*hjms+(6561.0/18304)*hjmh*hjms*hjms+(85293.0/98560)*hjms*hjms*hjms+(-6561.0/640640)*hjps*hjms*hjms+(37179.0/1281280)*hjph*hjms*hjms+(-2187.0/12320)*hjmh*hjps*hjms+(-6561.0/640640)*hjms*hjps*hjms+(282123.0/1281280)*hjps*hjps*hjms+(-2187.0/40040)*hjph*hjps*hjms+(243.0/7280)*hjmh*hjph*hjms+(37179.0/1281280)*hjms*hjph*hjms+(-2187.0/40040)*hjps*hjph*hjms+(-7047.0/256256)*hjph*hjph*hjms+(-543591.0/2562560)*hjmh*hjmh*hjps+(-2187.0/12320)*hjms*hjmh*hjps+(133407.0/2562560)*hjps*hjmh*hjps+(-243.0/32032)*hjph*hjmh*hjps+(-2187.0/12320)*hjmh*hjms*hjps+(-6561.0/640640)*hjms*hjms*hjps+(282123.0/1281280)*hjps*hjms*hjps+(-2187.0/40040)*hjph*hjms*hjps+(133407.0/2562560)*hjmh*hjps*hjps+(282123.0/1281280)*hjms*hjps*hjps+(478953.0/512512)*hjps*hjps*hjps+(6561.0/183040)*hjph*hjps*hjps+(-243.0/32032)*hjmh*hjph*hjps+(-2187.0/40040)*hjms*hjph*hjps+(6561.0/183040)*hjps*hjph*hjps+(46899.0/640640)*hjph*hjph*hjps+(4563.0/98560)*hjmh*hjmh*hjph+(243.0/7280)*hjms*hjmh*hjph+(-243.0/32032)*hjps*hjmh*hjph+(6291.0/640640)*hjph*hjmh*hjph+(243.0/7280)*hjmh*hjms*hjph+(37179.0/1281280)*hjms*hjms*hjph+(-2187.0/40040)*hjps*hjms*hjph+(-7047.0/256256)*hjph*hjms*hjph+(-243.0/32032)*hjmh*hjps*hjph+(-2187.0/40040)*hjms*hjps*hjph+(6561.0/183040)*hjps*hjps*hjph+(46899.0/640640)*hjph*hjps*hjph+(6291.0/640640)*hjmh*hjph*hjph+(-7047.0/256256)*hjms*hjph*hjph+(46899.0/640640)*hjps*hjph*hjph+(50409.0/183040)*hjph*hjph*hjph)
    # h3uxe[1,2] = 2.0/dx*((-111429.0/183040)*hjmh*hjmh*hjmh+(-441531.0/2562560)*hjms*hjmh*hjmh+(97443.0/1281280)*hjps*hjmh*hjmh+(-9045.0/512512)*hjph*hjmh*hjmh+(-441531.0/2562560)*hjmh*hjms*hjmh+(-32805.0/512512)*hjms*hjms*hjmh+(155277.0/2562560)*hjps*hjms*hjmh+(-17253.0/2562560)*hjph*hjms*hjmh+(97443.0/1281280)*hjmh*hjps*hjmh+(155277.0/2562560)*hjms*hjps*hjmh+(-2187.0/1281280)*hjps*hjps*hjmh+(-17253.0/2562560)*hjph*hjps*hjmh+(-9045.0/512512)*hjmh*hjph*hjmh+(-17253.0/2562560)*hjms*hjph*hjmh+(-17253.0/2562560)*hjps*hjph*hjmh+(-9045.0/512512)*hjph*hjph*hjmh+(-441531.0/2562560)*hjmh*hjmh*hjms+(-32805.0/512512)*hjms*hjmh*hjms+(155277.0/2562560)*hjps*hjmh*hjms+(-17253.0/2562560)*hjph*hjmh*hjms+(-32805.0/512512)*hjmh*hjms*hjms+(-662661.0/1281280)*hjms*hjms*hjms+(-439587.0/2562560)*hjps*hjms*hjms+(-2187.0/1281280)*hjph*hjms*hjms+(155277.0/2562560)*hjmh*hjps*hjms+(-439587.0/2562560)*hjms*hjps*hjms+(-439587.0/2562560)*hjps*hjps*hjms+(155277.0/2562560)*hjph*hjps*hjms+(-17253.0/2562560)*hjmh*hjph*hjms+(-2187.0/1281280)*hjms*hjph*hjms+(155277.0/2562560)*hjps*hjph*hjms+(97443.0/1281280)*hjph*hjph*hjms+(97443.0/1281280)*hjmh*hjmh*hjps+(155277.0/2562560)*hjms*hjmh*hjps+(-2187.0/1281280)*hjps*hjmh*hjps+(-17253.0/2562560)*hjph*hjmh*hjps+(155277.0/2562560)*hjmh*hjms*hjps+(-439587.0/2562560)*hjms*hjms*hjps+(-439587.0/2562560)*hjps*hjms*hjps+(155277.0/2562560)*hjph*hjms*hjps+(-2187.0/1281280)*hjmh*hjps*hjps+(-439587.0/2562560)*hjms*hjps*hjps+(-662661.0/1281280)*hjps*hjps*hjps+(-32805.0/512512)*hjph*hjps*hjps+(-17253.0/2562560)*hjmh*hjph*hjps+(155277.0/2562560)*hjms*hjph*hjps+(-32805.0/512512)*hjps*hjph*hjps+(-441531.0/2562560)*hjph*hjph*hjps+(-9045.0/512512)*hjmh*hjmh*hjph+(-17253.0/2562560)*hjms*hjmh*hjph+(-17253.0/2562560)*hjps*hjmh*hjph+(-9045.0/512512)*hjph*hjmh*hjph+(-17253.0/2562560)*hjmh*hjms*hjph+(-2187.0/1281280)*hjms*hjms*hjph+(155277.0/2562560)*hjps*hjms*hjph+(97443.0/1281280)*hjph*hjms*hjph+(-17253.0/2562560)*hjmh*hjps*hjph+(155277.0/2562560)*hjms*hjps*hjph+(-32805.0/512512)*hjps*hjps*hjph+(-441531.0/2562560)*hjph*hjps*hjph+(-9045.0/512512)*hjmh*hjph*hjph+(97443.0/1281280)*hjms*hjph*hjph+(-441531.0/2562560)*hjps*hjph*hjph+(-111429.0/183040)*hjph*hjph*hjph)
    # h3uxe[1,3] = 2.0/dx*((47763.0/366080)*hjmh*hjmh*hjmh+(4023.0/116480)*hjms*hjmh*hjmh+(-3699.0/232960)*hjps*hjmh*hjmh+(5169.0/1281280)*hjph*hjmh*hjmh+(4023.0/116480)*hjmh*hjms*hjmh+(4617.0/320320)*hjms*hjms*hjmh+(-243.0/29120)*hjps*hjms*hjmh+(-81.0/91520)*hjph*hjms*hjmh+(-3699.0/232960)*hjmh*hjps*hjmh+(-243.0/29120)*hjms*hjps*hjmh+(243.0/232960)*hjps*hjps*hjmh+(729.0/160160)*hjph*hjps*hjmh+(5169.0/1281280)*hjmh*hjph*hjmh+(-81.0/91520)*hjms*hjph*hjmh+(729.0/160160)*hjps*hjph*hjmh+(57.0/4928)*hjph*hjph*hjmh+(4023.0/116480)*hjmh*hjmh*hjms+(4617.0/320320)*hjms*hjmh*hjms+(-243.0/29120)*hjps*hjmh*hjms+(-81.0/91520)*hjph*hjmh*hjms+(4617.0/320320)*hjmh*hjms*hjms+(6561.0/116480)*hjms*hjms*hjms+(2187.0/320320)*hjps*hjms*hjms+(14337.0/1281280)*hjph*hjms*hjms+(-243.0/29120)*hjmh*hjps*hjms+(2187.0/320320)*hjms*hjps*hjms+(2187.0/116480)*hjps*hjps*hjms+(-243.0/9856)*hjph*hjps*hjms+(-81.0/91520)*hjmh*hjph*hjms+(14337.0/1281280)*hjms*hjph*hjms+(-243.0/9856)*hjps*hjph*hjms+(-67257.0/1281280)*hjph*hjph*hjms+(-3699.0/232960)*hjmh*hjmh*hjps+(-243.0/29120)*hjms*hjmh*hjps+(243.0/232960)*hjps*hjmh*hjps+(729.0/160160)*hjph*hjmh*hjps+(-243.0/29120)*hjmh*hjms*hjps+(2187.0/320320)*hjms*hjms*hjps+(2187.0/116480)*hjps*hjms*hjps+(-243.0/9856)*hjph*hjms*hjps+(243.0/232960)*hjmh*hjps*hjps+(2187.0/116480)*hjms*hjps*hjps+(-59049.0/197120)*hjps*hjps*hjps+(55161.0/1281280)*hjph*hjps*hjps+(729.0/160160)*hjmh*hjph*hjps+(-243.0/9856)*hjms*hjph*hjps+(55161.0/1281280)*hjps*hjph*hjps+(36261.0/320320)*hjph*hjph*hjps+(5169.0/1281280)*hjmh*hjmh*hjph+(-81.0/91520)*hjms*hjmh*hjph+(729.0/160160)*hjps*hjmh*hjph+(57.0/4928)*hjph*hjmh*hjph+(-81.0/91520)*hjmh*hjms*hjph+(14337.0/1281280)*hjms*hjms*hjph+(-243.0/9856)*hjps*hjms*hjph+(-67257.0/1281280)*hjph*hjms*hjph+(729.0/160160)*hjmh*hjps*hjph+(-243.0/9856)*hjms*hjps*hjph+(55161.0/1281280)*hjps*hjps*hjph+(36261.0/320320)*hjph*hjps*hjph+(57.0/4928)*hjmh*hjph*hjph+(-67257.0/1281280)*hjms*hjph*hjph+(36261.0/320320)*hjps*hjph*hjph+(503469.0/1281280)*hjph*hjph*hjph)

    # h3uxe[2,0] = 2.0/dx*((503469.0/1281280)*hjmh*hjmh*hjmh+(36261.0/320320)*hjms*hjmh*hjmh+(-67257.0/1281280)*hjps*hjmh*hjmh+(57.0/4928)*hjph*hjmh*hjmh+(36261.0/320320)*hjmh*hjms*hjmh+(55161.0/1281280)*hjms*hjms*hjmh+(-243.0/9856)*hjps*hjms*hjmh+(729.0/160160)*hjph*hjms*hjmh+(-67257.0/1281280)*hjmh*hjps*hjmh+(-243.0/9856)*hjms*hjps*hjmh+(14337.0/1281280)*hjps*hjps*hjmh+(-81.0/91520)*hjph*hjps*hjmh+(57.0/4928)*hjmh*hjph*hjmh+(729.0/160160)*hjms*hjph*hjmh+(-81.0/91520)*hjps*hjph*hjmh+(5169.0/1281280)*hjph*hjph*hjmh+(36261.0/320320)*hjmh*hjmh*hjms+(55161.0/1281280)*hjms*hjmh*hjms+(-243.0/9856)*hjps*hjmh*hjms+(729.0/160160)*hjph*hjmh*hjms+(55161.0/1281280)*hjmh*hjms*hjms+(-59049.0/197120)*hjms*hjms*hjms+(2187.0/116480)*hjps*hjms*hjms+(243.0/232960)*hjph*hjms*hjms+(-243.0/9856)*hjmh*hjps*hjms+(2187.0/116480)*hjms*hjps*hjms+(2187.0/320320)*hjps*hjps*hjms+(-243.0/29120)*hjph*hjps*hjms+(729.0/160160)*hjmh*hjph*hjms+(243.0/232960)*hjms*hjph*hjms+(-243.0/29120)*hjps*hjph*hjms+(-3699.0/232960)*hjph*hjph*hjms+(-67257.0/1281280)*hjmh*hjmh*hjps+(-243.0/9856)*hjms*hjmh*hjps+(14337.0/1281280)*hjps*hjmh*hjps+(-81.0/91520)*hjph*hjmh*hjps+(-243.0/9856)*hjmh*hjms*hjps+(2187.0/116480)*hjms*hjms*hjps+(2187.0/320320)*hjps*hjms*hjps+(-243.0/29120)*hjph*hjms*hjps+(14337.0/1281280)*hjmh*hjps*hjps+(2187.0/320320)*hjms*hjps*hjps+(6561.0/116480)*hjps*hjps*hjps+(4617.0/320320)*hjph*hjps*hjps+(-81.0/91520)*hjmh*hjph*hjps+(-243.0/29120)*hjms*hjph*hjps+(4617.0/320320)*hjps*hjph*hjps+(4023.0/116480)*hjph*hjph*hjps+(57.0/4928)*hjmh*hjmh*hjph+(729.0/160160)*hjms*hjmh*hjph+(-81.0/91520)*hjps*hjmh*hjph+(5169.0/1281280)*hjph*hjmh*hjph+(729.0/160160)*hjmh*hjms*hjph+(243.0/232960)*hjms*hjms*hjph+(-243.0/29120)*hjps*hjms*hjph+(-3699.0/232960)*hjph*hjms*hjph+(-81.0/91520)*hjmh*hjps*hjph+(-243.0/29120)*hjms*hjps*hjph+(4617.0/320320)*hjps*hjps*hjph+(4023.0/116480)*hjph*hjps*hjph+(5169.0/1281280)*hjmh*hjph*hjph+(-3699.0/232960)*hjms*hjph*hjph+(4023.0/116480)*hjps*hjph*hjph+(47763.0/366080)*hjph*hjph*hjph)
    # h3uxe[2,1] = 2.0/dx*((-111429.0/183040)*hjmh*hjmh*hjmh+(-441531.0/2562560)*hjms*hjmh*hjmh+(97443.0/1281280)*hjps*hjmh*hjmh+(-9045.0/512512)*hjph*hjmh*hjmh+(-441531.0/2562560)*hjmh*hjms*hjmh+(-32805.0/512512)*hjms*hjms*hjmh+(155277.0/2562560)*hjps*hjms*hjmh+(-17253.0/2562560)*hjph*hjms*hjmh+(97443.0/1281280)*hjmh*hjps*hjmh+(155277.0/2562560)*hjms*hjps*hjmh+(-2187.0/1281280)*hjps*hjps*hjmh+(-17253.0/2562560)*hjph*hjps*hjmh+(-9045.0/512512)*hjmh*hjph*hjmh+(-17253.0/2562560)*hjms*hjph*hjmh+(-17253.0/2562560)*hjps*hjph*hjmh+(-9045.0/512512)*hjph*hjph*hjmh+(-441531.0/2562560)*hjmh*hjmh*hjms+(-32805.0/512512)*hjms*hjmh*hjms+(155277.0/2562560)*hjps*hjmh*hjms+(-17253.0/2562560)*hjph*hjmh*hjms+(-32805.0/512512)*hjmh*hjms*hjms+(-662661.0/1281280)*hjms*hjms*hjms+(-439587.0/2562560)*hjps*hjms*hjms+(-2187.0/1281280)*hjph*hjms*hjms+(155277.0/2562560)*hjmh*hjps*hjms+(-439587.0/2562560)*hjms*hjps*hjms+(-439587.0/2562560)*hjps*hjps*hjms+(155277.0/2562560)*hjph*hjps*hjms+(-17253.0/2562560)*hjmh*hjph*hjms+(-2187.0/1281280)*hjms*hjph*hjms+(155277.0/2562560)*hjps*hjph*hjms+(97443.0/1281280)*hjph*hjph*hjms+(97443.0/1281280)*hjmh*hjmh*hjps+(155277.0/2562560)*hjms*hjmh*hjps+(-2187.0/1281280)*hjps*hjmh*hjps+(-17253.0/2562560)*hjph*hjmh*hjps+(155277.0/2562560)*hjmh*hjms*hjps+(-439587.0/2562560)*hjms*hjms*hjps+(-439587.0/2562560)*hjps*hjms*hjps+(155277.0/2562560)*hjph*hjms*hjps+(-2187.0/1281280)*hjmh*hjps*hjps+(-439587.0/2562560)*hjms*hjps*hjps+(-662661.0/1281280)*hjps*hjps*hjps+(-32805.0/512512)*hjph*hjps*hjps+(-17253.0/2562560)*hjmh*hjph*hjps+(155277.0/2562560)*hjms*hjph*hjps+(-32805.0/512512)*hjps*hjph*hjps+(-441531.0/2562560)*hjph*hjph*hjps+(-9045.0/512512)*hjmh*hjmh*hjph+(-17253.0/2562560)*hjms*hjmh*hjph+(-17253.0/2562560)*hjps*hjmh*hjph+(-9045.0/512512)*hjph*hjmh*hjph+(-17253.0/2562560)*hjmh*hjms*hjph+(-2187.0/1281280)*hjms*hjms*hjph+(155277.0/2562560)*hjps*hjms*hjph+(97443.0/1281280)*hjph*hjms*hjph+(-17253.0/2562560)*hjmh*hjps*hjph+(155277.0/2562560)*hjms*hjps*hjph+(-32805.0/512512)*hjps*hjps*hjph+(-441531.0/2562560)*hjph*hjps*hjph+(-9045.0/512512)*hjmh*hjph*hjph+(97443.0/1281280)*hjms*hjph*hjph+(-441531.0/2562560)*hjps*hjph*hjph+(-111429.0/183040)*hjph*hjph*hjph)
    # h3uxe[2,2] = 2.0/dx*((50409.0/183040)*hjmh*hjmh*hjmh+(46899.0/640640)*hjms*hjmh*hjmh+(-7047.0/256256)*hjps*hjmh*hjmh+(6291.0/640640)*hjph*hjmh*hjmh+(46899.0/640640)*hjmh*hjms*hjmh+(6561.0/183040)*hjms*hjms*hjmh+(-2187.0/40040)*hjps*hjms*hjmh+(-243.0/32032)*hjph*hjms*hjmh+(-7047.0/256256)*hjmh*hjps*hjmh+(-2187.0/40040)*hjms*hjps*hjmh+(37179.0/1281280)*hjps*hjps*hjmh+(243.0/7280)*hjph*hjps*hjmh+(6291.0/640640)*hjmh*hjph*hjmh+(-243.0/32032)*hjms*hjph*hjmh+(243.0/7280)*hjps*hjph*hjmh+(4563.0/98560)*hjph*hjph*hjmh+(46899.0/640640)*hjmh*hjmh*hjms+(6561.0/183040)*hjms*hjmh*hjms+(-2187.0/40040)*hjps*hjmh*hjms+(-243.0/32032)*hjph*hjmh*hjms+(6561.0/183040)*hjmh*hjms*hjms+(478953.0/512512)*hjms*hjms*hjms+(282123.0/1281280)*hjps*hjms*hjms+(133407.0/2562560)*hjph*hjms*hjms+(-2187.0/40040)*hjmh*hjps*hjms+(282123.0/1281280)*hjms*hjps*hjms+(-6561.0/640640)*hjps*hjps*hjms+(-2187.0/12320)*hjph*hjps*hjms+(-243.0/32032)*hjmh*hjph*hjms+(133407.0/2562560)*hjms*hjph*hjms+(-2187.0/12320)*hjps*hjph*hjms+(-543591.0/2562560)*hjph*hjph*hjms+(-7047.0/256256)*hjmh*hjmh*hjps+(-2187.0/40040)*hjms*hjmh*hjps+(37179.0/1281280)*hjps*hjmh*hjps+(243.0/7280)*hjph*hjmh*hjps+(-2187.0/40040)*hjmh*hjms*hjps+(282123.0/1281280)*hjms*hjms*hjps+(-6561.0/640640)*hjps*hjms*hjps+(-2187.0/12320)*hjph*hjms*hjps+(37179.0/1281280)*hjmh*hjps*hjps+(-6561.0/640640)*hjms*hjps*hjps+(85293.0/98560)*hjps*hjps*hjps+(6561.0/18304)*hjph*hjps*hjps+(243.0/7280)*hjmh*hjph*hjps+(-2187.0/12320)*hjms*hjph*hjps+(6561.0/18304)*hjps*hjph*hjps+(619407.0/1281280)*hjph*hjph*hjps+(6291.0/640640)*hjmh*hjmh*hjph+(-243.0/32032)*hjms*hjmh*hjph+(243.0/7280)*hjps*hjmh*hjph+(4563.0/98560)*hjph*hjmh*hjph+(-243.0/32032)*hjmh*hjms*hjph+(133407.0/2562560)*hjms*hjms*hjph+(-2187.0/12320)*hjps*hjms*hjph+(-543591.0/2562560)*hjph*hjms*hjph+(243.0/7280)*hjmh*hjps*hjph+(-2187.0/12320)*hjms*hjps*hjph+(6561.0/18304)*hjps*hjps*hjph+(619407.0/1281280)*hjph*hjps*hjph+(4563.0/98560)*hjmh*hjph*hjph+(-543591.0/2562560)*hjms*hjph*hjph+(619407.0/1281280)*hjps*hjph*hjph+(323217.0/232960)*hjph*hjph*hjph)
    # h3uxe[2,3] = 2.0/dx*((-6939.0/116480)*hjmh*hjmh*hjmh+(-2781.0/197120)*hjms*hjmh*hjmh+(459.0/116480)*hjps*hjmh*hjmh+(-9579.0/2562560)*hjph*hjmh*hjmh+(-2781.0/197120)*hjmh*hjms*hjmh+(-38151.0/2562560)*hjms*hjms*hjmh+(47871.0/2562560)*hjps*hjms*hjmh+(25029.0/2562560)*hjph*hjms*hjmh+(459.0/116480)*hjmh*hjps*hjmh+(47871.0/2562560)*hjms*hjps*hjmh+(-7047.0/183040)*hjps*hjps*hjmh+(-13203.0/512512)*hjph*hjps*hjmh+(-9579.0/2562560)*hjmh*hjph*hjmh+(25029.0/2562560)*hjms*hjph*hjmh+(-13203.0/512512)*hjps*hjph*hjmh+(-83751.0/2562560)*hjph*hjph*hjmh+(-2781.0/197120)*hjmh*hjmh*hjms+(-38151.0/2562560)*hjms*hjmh*hjms+(47871.0/2562560)*hjps*hjmh*hjms+(25029.0/2562560)*hjph*hjmh*hjms+(-38151.0/2562560)*hjmh*hjms*hjms+(-150903.0/1281280)*hjms*hjms*hjms+(-172773.0/2562560)*hjps*hjms*hjms+(-65853.0/1281280)*hjph*hjms*hjms+(47871.0/2562560)*hjmh*hjps*hjms+(-172773.0/2562560)*hjms*hjps*hjms+(89667.0/512512)*hjps*hjps*hjms+(321003.0/2562560)*hjph*hjps*hjms+(25029.0/2562560)*hjmh*hjph*hjms+(-65853.0/1281280)*hjms*hjph*hjms+(321003.0/2562560)*hjps*hjph*hjms+(194697.0/1281280)*hjph*hjph*hjms+(459.0/116480)*hjmh*hjmh*hjps+(47871.0/2562560)*hjms*hjmh*hjps+(-7047.0/183040)*hjps*hjmh*hjps+(-13203.0/512512)*hjph*hjmh*hjps+(47871.0/2562560)*hjmh*hjms*hjps+(-172773.0/2562560)*hjms*hjms*hjps+(89667.0/512512)*hjps*hjms*hjps+(321003.0/2562560)*hjph*hjms*hjps+(-7047.0/183040)*hjmh*hjps*hjps+(89667.0/512512)*hjms*hjps*hjps+(-518319.0/1281280)*hjps*hjps*hjps+(-791451.0/2562560)*hjph*hjps*hjps+(-13203.0/512512)*hjmh*hjph*hjps+(321003.0/2562560)*hjms*hjph*hjps+(-791451.0/2562560)*hjps*hjph*hjps+(-885789.0/2562560)*hjph*hjph*hjps+(-9579.0/2562560)*hjmh*hjmh*hjph+(25029.0/2562560)*hjms*hjmh*hjph+(-13203.0/512512)*hjps*hjmh*hjph+(-83751.0/2562560)*hjph*hjmh*hjph+(25029.0/2562560)*hjmh*hjms*hjph+(-65853.0/1281280)*hjms*hjms*hjph+(321003.0/2562560)*hjps*hjms*hjph+(194697.0/1281280)*hjph*hjms*hjph+(-13203.0/512512)*hjmh*hjps*hjph+(321003.0/2562560)*hjms*hjps*hjph+(-791451.0/2562560)*hjps*hjps*hjph+(-885789.0/2562560)*hjph*hjps*hjph+(-83751.0/2562560)*hjmh*hjph*hjph+(194697.0/1281280)*hjms*hjph*hjph+(-885789.0/2562560)*hjps*hjph*hjph+(-1164861.0/1281280)*hjph*hjph*hjph)
    
    # h3uxe[3,0] = 2.0/dx*((-107519.0/1281280)*hjmh*hjmh*hjmh+(-57951.0/2562560)*hjms*hjmh*hjmh+(2727.0/256256)*hjps*hjmh*hjmh+(-18559.0/7687680)*hjph*hjmh*hjmh+(-57951.0/2562560)*hjmh*hjms*hjmh+(-2403.0/366080)*hjms*hjms*hjmh+(10233.0/2562560)*hjps*hjms*hjmh+(-1137.0/2562560)*hjph*hjms*hjmh+(2727.0/256256)*hjmh*hjps*hjmh+(10233.0/2562560)*hjms*hjps*hjmh+(-1647.0/1281280)*hjps*hjps*hjmh+(-1137.0/2562560)*hjph*hjps*hjmh+(-18559.0/7687680)*hjmh*hjph*hjmh+(-1137.0/2562560)*hjms*hjph*hjmh+(-1137.0/2562560)*hjps*hjph*hjmh+(-18559.0/7687680)*hjph*hjph*hjmh+(-57951.0/2562560)*hjmh*hjmh*hjms+(-2403.0/366080)*hjms*hjmh*hjms+(10233.0/2562560)*hjps*hjmh*hjms+(-1137.0/2562560)*hjph*hjmh*hjms+(-2403.0/366080)*hjmh*hjms*hjms+(78003.0/1281280)*hjms*hjms*hjms+(-2187.0/512512)*hjps*hjms*hjms+(-1647.0/1281280)*hjph*hjms*hjms+(10233.0/2562560)*hjmh*hjps*hjms+(-2187.0/512512)*hjms*hjps*hjms+(-2187.0/512512)*hjps*hjps*hjms+(10233.0/2562560)*hjph*hjps*hjms+(-1137.0/2562560)*hjmh*hjph*hjms+(-1647.0/1281280)*hjms*hjph*hjms+(10233.0/2562560)*hjps*hjph*hjms+(2727.0/256256)*hjph*hjph*hjms+(2727.0/256256)*hjmh*hjmh*hjps+(10233.0/2562560)*hjms*hjmh*hjps+(-1647.0/1281280)*hjps*hjmh*hjps+(-1137.0/2562560)*hjph*hjmh*hjps+(10233.0/2562560)*hjmh*hjms*hjps+(-2187.0/512512)*hjms*hjms*hjps+(-2187.0/512512)*hjps*hjms*hjps+(10233.0/2562560)*hjph*hjms*hjps+(-1647.0/1281280)*hjmh*hjps*hjps+(-2187.0/512512)*hjms*hjps*hjps+(78003.0/1281280)*hjps*hjps*hjps+(-2403.0/366080)*hjph*hjps*hjps+(-1137.0/2562560)*hjmh*hjph*hjps+(10233.0/2562560)*hjms*hjph*hjps+(-2403.0/366080)*hjps*hjph*hjps+(-57951.0/2562560)*hjph*hjph*hjps+(-18559.0/7687680)*hjmh*hjmh*hjph+(-1137.0/2562560)*hjms*hjmh*hjph+(-1137.0/2562560)*hjps*hjmh*hjph+(-18559.0/7687680)*hjph*hjmh*hjph+(-1137.0/2562560)*hjmh*hjms*hjph+(-1647.0/1281280)*hjms*hjms*hjph+(10233.0/2562560)*hjps*hjms*hjph+(2727.0/256256)*hjph*hjms*hjph+(-1137.0/2562560)*hjmh*hjps*hjph+(10233.0/2562560)*hjms*hjps*hjph+(-2403.0/366080)*hjps*hjps*hjph+(-57951.0/2562560)*hjph*hjps*hjph+(-18559.0/7687680)*hjmh*hjph*hjph+(2727.0/256256)*hjms*hjph*hjph+(-57951.0/2562560)*hjps*hjph*hjph+(-107519.0/1281280)*hjph*hjph*hjph)
    # h3uxe[3,1] = 2.0/dx*((47763.0/366080)*hjmh*hjmh*hjmh+(4023.0/116480)*hjms*hjmh*hjmh+(-3699.0/232960)*hjps*hjmh*hjmh+(5169.0/1281280)*hjph*hjmh*hjmh+(4023.0/116480)*hjmh*hjms*hjmh+(4617.0/320320)*hjms*hjms*hjmh+(-243.0/29120)*hjps*hjms*hjmh+(-81.0/91520)*hjph*hjms*hjmh+(-3699.0/232960)*hjmh*hjps*hjmh+(-243.0/29120)*hjms*hjps*hjmh+(243.0/232960)*hjps*hjps*hjmh+(729.0/160160)*hjph*hjps*hjmh+(5169.0/1281280)*hjmh*hjph*hjmh+(-81.0/91520)*hjms*hjph*hjmh+(729.0/160160)*hjps*hjph*hjmh+(57.0/4928)*hjph*hjph*hjmh+(4023.0/116480)*hjmh*hjmh*hjms+(4617.0/320320)*hjms*hjmh*hjms+(-243.0/29120)*hjps*hjmh*hjms+(-81.0/91520)*hjph*hjmh*hjms+(4617.0/320320)*hjmh*hjms*hjms+(6561.0/116480)*hjms*hjms*hjms+(2187.0/320320)*hjps*hjms*hjms+(14337.0/1281280)*hjph*hjms*hjms+(-243.0/29120)*hjmh*hjps*hjms+(2187.0/320320)*hjms*hjps*hjms+(2187.0/116480)*hjps*hjps*hjms+(-243.0/9856)*hjph*hjps*hjms+(-81.0/91520)*hjmh*hjph*hjms+(14337.0/1281280)*hjms*hjph*hjms+(-243.0/9856)*hjps*hjph*hjms+(-67257.0/1281280)*hjph*hjph*hjms+(-3699.0/232960)*hjmh*hjmh*hjps+(-243.0/29120)*hjms*hjmh*hjps+(243.0/232960)*hjps*hjmh*hjps+(729.0/160160)*hjph*hjmh*hjps+(-243.0/29120)*hjmh*hjms*hjps+(2187.0/320320)*hjms*hjms*hjps+(2187.0/116480)*hjps*hjms*hjps+(-243.0/9856)*hjph*hjms*hjps+(243.0/232960)*hjmh*hjps*hjps+(2187.0/116480)*hjms*hjps*hjps+(-59049.0/197120)*hjps*hjps*hjps+(55161.0/1281280)*hjph*hjps*hjps+(729.0/160160)*hjmh*hjph*hjps+(-243.0/9856)*hjms*hjph*hjps+(55161.0/1281280)*hjps*hjph*hjps+(36261.0/320320)*hjph*hjph*hjps+(5169.0/1281280)*hjmh*hjmh*hjph+(-81.0/91520)*hjms*hjmh*hjph+(729.0/160160)*hjps*hjmh*hjph+(57.0/4928)*hjph*hjmh*hjph+(-81.0/91520)*hjmh*hjms*hjph+(14337.0/1281280)*hjms*hjms*hjph+(-243.0/9856)*hjps*hjms*hjph+(-67257.0/1281280)*hjph*hjms*hjph+(729.0/160160)*hjmh*hjps*hjph+(-243.0/9856)*hjms*hjps*hjph+(55161.0/1281280)*hjps*hjps*hjph+(36261.0/320320)*hjph*hjps*hjph+(57.0/4928)*hjmh*hjph*hjph+(-67257.0/1281280)*hjms*hjph*hjph+(36261.0/320320)*hjps*hjph*hjph+(503469.0/1281280)*hjph*hjph*hjph)
    # h3uxe[3,2] = 2.0/dx*((-6939.0/116480)*hjmh*hjmh*hjmh+(-2781.0/197120)*hjms*hjmh*hjmh+(459.0/116480)*hjps*hjmh*hjmh+(-9579.0/2562560)*hjph*hjmh*hjmh+(-2781.0/197120)*hjmh*hjms*hjmh+(-38151.0/2562560)*hjms*hjms*hjmh+(47871.0/2562560)*hjps*hjms*hjmh+(25029.0/2562560)*hjph*hjms*hjmh+(459.0/116480)*hjmh*hjps*hjmh+(47871.0/2562560)*hjms*hjps*hjmh+(-7047.0/183040)*hjps*hjps*hjmh+(-13203.0/512512)*hjph*hjps*hjmh+(-9579.0/2562560)*hjmh*hjph*hjmh+(25029.0/2562560)*hjms*hjph*hjmh+(-13203.0/512512)*hjps*hjph*hjmh+(-83751.0/2562560)*hjph*hjph*hjmh+(-2781.0/197120)*hjmh*hjmh*hjms+(-38151.0/2562560)*hjms*hjmh*hjms+(47871.0/2562560)*hjps*hjmh*hjms+(25029.0/2562560)*hjph*hjmh*hjms+(-38151.0/2562560)*hjmh*hjms*hjms+(-150903.0/1281280)*hjms*hjms*hjms+(-172773.0/2562560)*hjps*hjms*hjms+(-65853.0/1281280)*hjph*hjms*hjms+(47871.0/2562560)*hjmh*hjps*hjms+(-172773.0/2562560)*hjms*hjps*hjms+(89667.0/512512)*hjps*hjps*hjms+(321003.0/2562560)*hjph*hjps*hjms+(25029.0/2562560)*hjmh*hjph*hjms+(-65853.0/1281280)*hjms*hjph*hjms+(321003.0/2562560)*hjps*hjph*hjms+(194697.0/1281280)*hjph*hjph*hjms+(459.0/116480)*hjmh*hjmh*hjps+(47871.0/2562560)*hjms*hjmh*hjps+(-7047.0/183040)*hjps*hjmh*hjps+(-13203.0/512512)*hjph*hjmh*hjps+(47871.0/2562560)*hjmh*hjms*hjps+(-172773.0/2562560)*hjms*hjms*hjps+(89667.0/512512)*hjps*hjms*hjps+(321003.0/2562560)*hjph*hjms*hjps+(-7047.0/183040)*hjmh*hjps*hjps+(89667.0/512512)*hjms*hjps*hjps+(-518319.0/1281280)*hjps*hjps*hjps+(-791451.0/2562560)*hjph*hjps*hjps+(-13203.0/512512)*hjmh*hjph*hjps+(321003.0/2562560)*hjms*hjph*hjps+(-791451.0/2562560)*hjps*hjph*hjps+(-885789.0/2562560)*hjph*hjph*hjps+(-9579.0/2562560)*hjmh*hjmh*hjph+(25029.0/2562560)*hjms*hjmh*hjph+(-13203.0/512512)*hjps*hjmh*hjph+(-83751.0/2562560)*hjph*hjmh*hjph+(25029.0/2562560)*hjmh*hjms*hjph+(-65853.0/1281280)*hjms*hjms*hjph+(321003.0/2562560)*hjps*hjms*hjph+(194697.0/1281280)*hjph*hjms*hjph+(-13203.0/512512)*hjmh*hjps*hjph+(321003.0/2562560)*hjms*hjps*hjph+(-791451.0/2562560)*hjps*hjps*hjph+(-885789.0/2562560)*hjph*hjps*hjph+(-83751.0/2562560)*hjmh*hjph*hjph+(194697.0/1281280)*hjms*hjph*hjph+(-885789.0/2562560)*hjps*hjph*hjph+(-1164861.0/1281280)*hjph*hjph*hjph)
    # h3uxe[3,3] = 2.0/dx*((953.0/73216)*hjmh*hjmh*hjmh+(2799.0/1281280)*hjms*hjmh*hjmh+(3321.0/2562560)*hjps*hjmh*hjmh+(1163.0/549120)*hjph*hjmh*hjmh+(2799.0/1281280)*hjmh*hjms*hjmh+(4509.0/640640)*hjms*hjms*hjmh+(-459.0/32032)*hjps*hjms*hjmh+(-2703.0/320320)*hjph*hjms*hjmh+(3321.0/2562560)*hjmh*hjps*hjmh+(-459.0/32032)*hjms*hjps*hjmh+(99279.0/2562560)*hjps*hjps*hjmh+(867.0/40040)*hjph*hjps*hjmh+(1163.0/549120)*hjmh*hjph*hjmh+(-2703.0/320320)*hjms*hjph*hjmh+(867.0/40040)*hjps*hjph*hjmh+(45223.0/1921920)*hjph*hjph*hjmh+(2799.0/1281280)*hjmh*hjmh*hjms+(4509.0/640640)*hjms*hjmh*hjms+(-459.0/32032)*hjps*hjmh*hjms+(-2703.0/320320)*hjph*hjmh*hjms+(4509.0/640640)*hjmh*hjms*hjms+(729.0/1281280)*hjms*hjms*hjms+(41553.0/640640)*hjps*hjms*hjms+(4833.0/116480)*hjph*hjms*hjms+(-459.0/32032)*hjmh*hjps*hjms+(41553.0/640640)*hjms*hjps*hjms+(-242757.0/1281280)*hjps*hjps*hjms+(-33507.0/320320)*hjph*hjps*hjms+(-2703.0/320320)*hjmh*hjph*hjms+(4833.0/116480)*hjms*hjph*hjms+(-33507.0/320320)*hjps*hjph*hjms+(-2565.0/23296)*hjph*hjph*hjms+(3321.0/2562560)*hjmh*hjmh*hjps+(-459.0/32032)*hjms*hjmh*hjps+(99279.0/2562560)*hjps*hjmh*hjps+(867.0/40040)*hjph*hjmh*hjps+(-459.0/32032)*hjmh*hjms*hjps+(41553.0/640640)*hjms*hjms*hjps+(-242757.0/1281280)*hjps*hjms*hjps+(-33507.0/320320)*hjph*hjms*hjps+(99279.0/2562560)*hjmh*hjps*hjps+(-242757.0/1281280)*hjms*hjps*hjps+(235467.0/366080)*hjps*hjps*hjps+(6345.0/23296)*hjph*hjps*hjps+(867.0/40040)*hjmh*hjph*hjps+(-33507.0/320320)*hjms*hjph*hjps+(6345.0/23296)*hjps*hjph*hjps+(163413.0/640640)*hjph*hjph*hjps+(1163.0/549120)*hjmh*hjmh*hjph+(-2703.0/320320)*hjms*hjmh*hjph+(867.0/40040)*hjps*hjmh*hjph+(45223.0/1921920)*hjph*hjmh*hjph+(-2703.0/320320)*hjmh*hjms*hjph+(4833.0/116480)*hjms*hjms*hjph+(-33507.0/320320)*hjps*hjms*hjph+(-2565.0/23296)*hjph*hjms*hjph+(867.0/40040)*hjmh*hjps*hjph+(-33507.0/320320)*hjms*hjps*hjph+(6345.0/23296)*hjps*hjps*hjph+(163413.0/640640)*hjph*hjps*hjph+(45223.0/1921920)*hjmh*hjph*hjph+(-2565.0/23296)*hjms*hjph*hjph+(163413.0/640640)*hjps*hjph*hjph+(5377.0/8960)*hjph*hjph*hjph)
        
    # h3uxe = 1.0/3.0*dot(h3uxe,ue)
    
    print(Ge)
    print(uhe)
    print()
    return Ge, uhe

def FEMforG(xn, hp,up,dx):
    
    n = len(xn)
    A = zeros((3*n+1,3*n+1))
    b = zeros((3*n+1,3*n+1))
    
    # i = 0
    i=0
    A[3*i : 3*i + 4, 3*i : 3*i + 4] =  eye(4,4)
    b[3*i : 3*i + 4, 3*i : 3*i + 4] =  hp[0]*up[0]*eye(4,4)
 
    # i = 0
    i=n-1
    A[3*i : 3*i + 4, 3*i : 3*i + 4] =  eye(4,4)
    b[3*i : 3*i + 4, 3*i : 3*i + 4] =  hp[-1]*up[-1]*eye(4,4)
    
    for i in range(1,n-1):
        #elementwisematrices
        #Ge, uhe, h3uxe = FEMfroGElem(hp,up,dx,i)
        
        Ge, uhe = FEMfroGElem(hp,up,dx,i)
        
        A[3*i : 3*i + 4, 3*i : 3*i + 4] =  A[3*i : 3*i + 4, 3*i : 3*i + 4] + Ge
        b[3*i : 3*i + 4, 3*i : 3*i + 4] =  b[3*i : 3*i + 4, 3*i : 3*i + 4]  +  uhe #+ h3uxe


    
    #print(A)
    #print(b)
    
    #print(Ge)
    #print(uhe)

    #LHS = dot(b,up)
    #G = solve(A,LHS)
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
    hjmh = 2 #+ a1*exp(-sqrt(3)*abs(xjmh))
    hjms = 2 #+ a1*exp(-sqrt(3)*abs(xjms))
    hjps = 2 #+ a1*exp(-sqrt(3)*abs(xjps))
    hjph = 2 #+ a1*exp(-sqrt(3)*abs(xjph))
    hcell = [hjmh,hjms,hjps,hjph]
    xhcell = [xjmh,xjms,xjps,xjph]
    hp = hp + hcell
    xp = xp + xhcell
    
    #u      
    ujmh = 1
    ujms = 1
    ujps = 1
    ujph = 1
    
    ucell = [ujmh,ujms,ujps,ujph]
    up = up + ucell
    
    for i in range(1,n):
        xjms = x[i] - dx/6
        xj = x[i]
        xjps = x[i] + dx/6
        xjph = x[i] + dx/2
        
        #h
        hjms = 2# + a1*exp(-sqrt(3)*abs(xjms))
        hjps = 2# + a1*exp(-sqrt(3)*abs(xjps))
        hjph = 2# + a1*exp(-sqrt(3)*abs(xjph))
        hcell = [hjms,hjps,hjph]
        xhcell = [xjms,xjps,xjph]
        hp = hp + hcell
        xp = xp + xhcell
        
        #u      
        ujms = 1
        ujps =1
        ujph = 1
        #ujph = c*(hjph - 1)/hjph 
        
        ucell = [ujms,ujps,ujph]
        up = up + ucell
        
        
    return hp,up,xp

n= 2
dx = 10.0/n
xn = arange(-5,5+ dx,dx)
a1 = 0.5


h,u,x = initialhu(xn,dx,a1)

Guh = array(h)*array(u)

A,b = FEMforG(xn,h,u,dx)

Gc = solve(A, b.dot(u))