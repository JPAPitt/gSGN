"""
Python code to generate csv files that contain integrals of the basis functions
needed to generate the matrices of the Finite Element Method.

by Jordan Pitt 11/10/2018
"""

"""
############################# IMPORTS #########################################
"""


from sympy import * #Related website: https://www.sympy.org/en/index.html
from IPython.display import display

# symbols
x,t,dx,dt,ga,b1 = symbols('x,t,dx,dt,g,beta1', positive = True, nonzero = True,real=True)

def PolyFromPoints(yjmh,yjms,yjps,yjph,dx):
    a3 = (-9*yjmh + 27*yjms - 27*yjps + 9*yjph)/ (2*dx**3)
    a2 = (9*yjmh - 9*yjms - 9*yjps + 9*yjph  )/ (4*dx**2)
    a1 = (yjmh  - 27*yjms  + 27*yjps  - yjph  )/ (8*dx)
    a0 = (-yjmh+ 9*yjms + 9*yjps  - yjph)/ 16
    return a0,a1,a2,a3

def SolveForuEdges(beta1,ha0,ha1,ha2,ha3,Ga0,Ga1,Ga2,Ga3 ,dx):
    print(beta1,ha0,ha1,ha2,ha3,Ga0,Ga1,Ga2,Ga3 ,dx)
    if (abs(beta1) > 10.0**(-10)):
        if (abs(ha3) < 10.0**(-15)):
            ua3 =(108*Ga0*beta1**2*ha0**4*ha1*ha3**2 - 36*Ga0*beta1**2*ha0**4*ha2**2*ha3 + 576*Ga0*beta1**2*ha0**3*ha1**2*ha2*ha3 - 360*Ga0*beta1**2*ha0**3*ha1*ha2**3 - 216*Ga0*beta1**2*ha0**2*ha1**4*ha3 + 240*Ga0*beta1**2*ha0**2*ha1**3*ha2**2 + 24*Ga0*beta1**2*ha0*ha1**5*ha2 + 24*Ga0*beta1**2*ha1**7 - 21*Ga0*beta1*ha0**2*ha1**2*ha3 + 72*Ga0*beta1*ha0**2*ha1*ha2**2 - 22*Ga0*beta1*ha0*ha1**3*ha2 - 11*Ga0*beta1*ha1**5 + Ga0*ha0**2*ha3 - 2*Ga0*ha0*ha1*ha2 + Ga0*ha1**3 - 198*Ga1*beta1**2*ha0**5*ha3**2 - 426*Ga1*beta1**2*ha0**4*ha1*ha2*ha3 + 192*Ga1*beta1**2*ha0**4*ha2**3 + 240*Ga1*beta1**2*ha0**3*ha1**3*ha3 - 312*Ga1*beta1**2*ha0**3*ha1**2*ha2**2 - 48*Ga1*beta1**2*ha0**2*ha1**4*ha2 - 24*Ga1*beta1**2*ha0*ha1**6 + 10*Ga1*beta1*ha0**3*ha1*ha3 - 28*Ga1*beta1*ha0**3*ha2**2 + 33*Ga1*beta1*ha0**2*ha1**2*ha2 + 11*Ga1*beta1*ha0*ha1**4 + Ga1*ha0**2*ha2 - Ga1*ha0*ha1**2 + 132*Ga2*beta1**2*ha0**5*ha2*ha3 - 144*Ga2*beta1**2*ha0**4*ha1**2*ha3 + 168*Ga2*beta1**2*ha0**4*ha1*ha2**2 + 72*Ga2*beta1**2*ha0**3*ha1**3*ha2 + 24*Ga2*beta1**2*ha0**2*ha1**5 - 22*Ga2*beta1*ha0**4*ha3 - 44*Ga2*beta1*ha0**3*ha1*ha2 - 11*Ga2*beta1*ha0**2*ha1**3 + Ga2*ha0**2*ha1 + 90*Ga3*beta1**2*ha0**5*ha1*ha3 - 96*Ga3*beta1**2*ha0**5*ha2**2 - 6*Ga3*beta1**2*ha0**4*ha1**2*ha2 - 24*Ga3*beta1**2*ha0**3*ha1**4 + 22*Ga3*beta1*ha0**4*ha2 + 11*Ga3*beta1*ha0**3*ha1**2 - Ga3*ha0**3)/(ha0**4*(1188*beta1**3*ha0**4*ha3**2 - 4248*beta1**3*ha0**3*ha1*ha2*ha3 + 2304*beta1**3*ha0**3*ha2**3 - 792*beta1**3*ha0**2*ha1**3*ha3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 + 624*beta1**2*ha0**2*ha1*ha3 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
            ua2 = (54*Ga0*beta1**2*ha0**4*ha3**2 - 522*Ga0*beta1**2*ha0**3*ha1*ha2*ha3 + 144*Ga0*beta1**2*ha0**3*ha2**3 + 180*Ga0*beta1**2*ha0**2*ha1**3*ha3 + 36*Ga0*beta1**2*ha0**2*ha1**2*ha2**2 - 216*Ga0*beta1**2*ha0*ha1**4*ha2 - 108*Ga0*beta1**2*ha1**6 + 30*Ga0*beta1*ha0**2*ha1*ha3 - 36*Ga0*beta1*ha0**2*ha2**2 - 15*Ga0*beta1*ha0*ha1**2*ha2 + 39*Ga0*beta1*ha1**4 + Ga0*ha0*ha2 - Ga0*ha1**2 + 324*Ga1*beta1**2*ha0**4*ha2*ha3 - 243*Ga1*beta1**2*ha0**3*ha1**2*ha3 + 216*Ga1*beta1**2*ha0**3*ha1*ha2**2 + 324*Ga1*beta1**2*ha0**2*ha1**3*ha2 + 108*Ga1*beta1**2*ha0*ha1**5 - 9*Ga1*beta1*ha0**3*ha3 - 24*Ga1*beta1*ha0**2*ha1*ha2 - 39*Ga1*beta1*ha0*ha1**3 + Ga1*ha0*ha1 + 126*Ga2*beta1**2*ha0**4*ha1*ha3 - 144*Ga2*beta1**2*ha0**4*ha2**2 - 252*Ga2*beta1**2*ha0**3*ha1**2*ha2 - 108*Ga2*beta1**2*ha0**2*ha1**4 + 36*Ga2*beta1*ha0**3*ha2 + 39*Ga2*beta1*ha0**2*ha1**2 - Ga2*ha0**2 - 54*Ga3*beta1**2*ha0**5*ha3 + 72*Ga3*beta1**2*ha0**4*ha1*ha2 + 63*Ga3*beta1**2*ha0**3*ha1**3 - 21*Ga3*beta1*ha0**3*ha1)/(ha0**3*(1188*beta1**3*ha0**4*ha3**2 - 4248*beta1**3*ha0**3*ha1*ha2*ha3 + 2304*beta1**3*ha0**3*ha2**3 - 792*beta1**3*ha0**2*ha1**3*ha3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 + 624*beta1**2*ha0**2*ha1*ha3 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
            ua1 = (36*Ga0*beta1**2*ha0**3*ha2*ha3 - 432*Ga0*beta1**2*ha0**2*ha1**2*ha3 + 504*Ga0*beta1**2*ha0**2*ha1*ha2**2 - 384*Ga0*beta1**2*ha0*ha1**3*ha2 + 432*Ga0*beta1**2*ha1**5 + 6*Ga0*beta1*ha0**2*ha3 - 48*Ga0*beta1*ha0*ha1*ha2 - 27*Ga0*beta1*ha1**3 + Ga0*ha1 + 594*Ga1*beta1**2*ha0**3*ha1*ha3 - 576*Ga1*beta1**2*ha0**3*ha2**2 + 72*Ga1*beta1**2*ha0**2*ha1**2*ha2 - 432*Ga1*beta1**2*ha0*ha1**4 + 52*Ga1*beta1*ha0**2*ha2 + 27*Ga1*beta1*ha0*ha1**2 - Ga1*ha0 - 132*Ga2*beta1**2*ha0**4*ha3 + 72*Ga2*beta1**2*ha0**3*ha1*ha2 + 312*Ga2*beta1**2*ha0**2*ha1**3 - 4*Ga2*beta1*ha0**2*ha1 + 96*Ga3*beta1**2*ha0**4*ha2 - 162*Ga3*beta1**2*ha0**3*ha1**2 - 6*Ga3*beta1*ha0**3)/(ha0**2*(1188*beta1**3*ha0**4*ha3**2 - 4248*beta1**3*ha0**3*ha1*ha2*ha3 + 2304*beta1**3*ha0**3*ha2**3 - 792*beta1**3*ha0**2*ha1**3*ha3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 + 624*beta1**2*ha0**2*ha1*ha3 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
            ua0 = (1296*Ga0*beta1**3*ha0**4*ha3**2 - 5184*Ga0*beta1**3*ha0**3*ha1*ha2*ha3 + 2592*Ga0*beta1**3*ha0**3*ha2**3 - 1728*Ga0*beta1**3*ha0**2*ha1**3*ha3 + 2592*Ga0*beta1**3*ha0**2*ha1**2*ha2**2 - 2160*Ga0*beta1**3*ha0*ha1**4*ha2 + 1440*Ga0*beta1**3*ha1**6 + 702*Ga0*beta1**2*ha0**2*ha1*ha3 - 792*Ga0*beta1**2*ha0**2*ha2**2 - 378*Ga0*beta1**2*ha0*ha1**2*ha2 - 192*Ga0*beta1**2*ha1**4 + 54*Ga0*beta1*ha0*ha2 + 27*Ga0*beta1*ha1**2 - Ga0 + 648*Ga1*beta1**3*ha0**4*ha2*ha3 + 1296*Ga1*beta1**3*ha0**3*ha1**2*ha3 - 1296*Ga1*beta1**3*ha0**3*ha1*ha2**2 + 864*Ga1*beta1**3*ha0**2*ha1**3*ha2 - 1080*Ga1*beta1**3*ha0*ha1**5 - 18*Ga1*beta1**2*ha0**3*ha3 + 108*Ga1*beta1**2*ha0**2*ha1*ha2 + 3*Ga1*beta1**2*ha0*ha1**3 - Ga1*beta1*ha0*ha1 - 144*Ga2*beta1**3*ha0**4*ha1*ha3 - 288*Ga2*beta1**3*ha0**4*ha2**2 - 288*Ga2*beta1**3*ha0**3*ha1**2*ha2 + 720*Ga2*beta1**3*ha0**2*ha1**4 + 72*Ga2*beta1**2*ha0**3*ha2 + 66*Ga2*beta1**2*ha0**2*ha1**2 - 2*Ga2*beta1*ha0**2 - 108*Ga3*beta1**3*ha0**5*ha3 + 432*Ga3*beta1**3*ha0**4*ha1*ha2 - 360*Ga3*beta1**3*ha0**3*ha1**3 - 60*Ga3*beta1**2*ha0**3*ha1)/(ha0*(1188*beta1**3*ha0**4*ha3**2 - 4248*beta1**3*ha0**3*ha1*ha2*ha3 + 2304*beta1**3*ha0**3*ha2**3 - 792*beta1**3*ha0**2*ha1**3*ha3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 + 624*beta1**2*ha0**2*ha1*ha3 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
        elif (abs(ha2) < 10.0**(-15)):
            ua3 = (-360*Ga0*beta1**2*ha0**3*ha1*ha2**3 + 240*Ga0*beta1**2*ha0**2*ha1**3*ha2**2 + 24*Ga0*beta1**2*ha0*ha1**5*ha2 + 24*Ga0*beta1**2*ha1**7 + 72*Ga0*beta1*ha0**2*ha1*ha2**2 - 22*Ga0*beta1*ha0*ha1**3*ha2 - 11*Ga0*beta1*ha1**5 - 2*Ga0*ha0*ha1*ha2 + Ga0*ha1**3 + 192*Ga1*beta1**2*ha0**4*ha2**3 - 312*Ga1*beta1**2*ha0**3*ha1**2*ha2**2 - 48*Ga1*beta1**2*ha0**2*ha1**4*ha2 - 24*Ga1*beta1**2*ha0*ha1**6 - 28*Ga1*beta1*ha0**3*ha2**2 + 33*Ga1*beta1*ha0**2*ha1**2*ha2 + 11*Ga1*beta1*ha0*ha1**4 + Ga1*ha0**2*ha2 - Ga1*ha0*ha1**2 + 168*Ga2*beta1**2*ha0**4*ha1*ha2**2 + 72*Ga2*beta1**2*ha0**3*ha1**3*ha2 + 24*Ga2*beta1**2*ha0**2*ha1**5 - 44*Ga2*beta1*ha0**3*ha1*ha2 - 11*Ga2*beta1*ha0**2*ha1**3 + Ga2*ha0**2*ha1 - 96*Ga3*beta1**2*ha0**5*ha2**2 - 6*Ga3*beta1**2*ha0**4*ha1**2*ha2 - 24*Ga3*beta1**2*ha0**3*ha1**4 + 22*Ga3*beta1*ha0**4*ha2 + 11*Ga3*beta1*ha0**3*ha1**2 - Ga3*ha0**3)/(ha0**4*(2304*beta1**3*ha0**3*ha2**3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
            ua2 = (144*Ga0*beta1**2*ha0**3*ha2**3 + 36*Ga0*beta1**2*ha0**2*ha1**2*ha2**2 - 216*Ga0*beta1**2*ha0*ha1**4*ha2 - 108*Ga0*beta1**2*ha1**6 - 36*Ga0*beta1*ha0**2*ha2**2 - 15*Ga0*beta1*ha0*ha1**2*ha2 + 39*Ga0*beta1*ha1**4 + Ga0*ha0*ha2 - Ga0*ha1**2 + 216*Ga1*beta1**2*ha0**3*ha1*ha2**2 + 324*Ga1*beta1**2*ha0**2*ha1**3*ha2 + 108*Ga1*beta1**2*ha0*ha1**5 - 24*Ga1*beta1*ha0**2*ha1*ha2 - 39*Ga1*beta1*ha0*ha1**3 + Ga1*ha0*ha1 - 144*Ga2*beta1**2*ha0**4*ha2**2 - 252*Ga2*beta1**2*ha0**3*ha1**2*ha2 - 108*Ga2*beta1**2*ha0**2*ha1**4 + 36*Ga2*beta1*ha0**3*ha2 + 39*Ga2*beta1*ha0**2*ha1**2 - Ga2*ha0**2 + 72*Ga3*beta1**2*ha0**4*ha1*ha2 + 63*Ga3*beta1**2*ha0**3*ha1**3 - 21*Ga3*beta1*ha0**3*ha1)/(ha0**3*(2304*beta1**3*ha0**3*ha2**3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
            ua1 = (504*Ga0*beta1**2*ha0**2*ha1*ha2**2 - 384*Ga0*beta1**2*ha0*ha1**3*ha2 + 432*Ga0*beta1**2*ha1**5 - 48*Ga0*beta1*ha0*ha1*ha2 - 27*Ga0*beta1*ha1**3 + Ga0*ha1 - 576*Ga1*beta1**2*ha0**3*ha2**2 + 72*Ga1*beta1**2*ha0**2*ha1**2*ha2 - 432*Ga1*beta1**2*ha0*ha1**4 + 52*Ga1*beta1*ha0**2*ha2 + 27*Ga1*beta1*ha0*ha1**2 - Ga1*ha0 + 72*Ga2*beta1**2*ha0**3*ha1*ha2 + 312*Ga2*beta1**2*ha0**2*ha1**3 - 4*Ga2*beta1*ha0**2*ha1 + 96*Ga3*beta1**2*ha0**4*ha2 - 162*Ga3*beta1**2*ha0**3*ha1**2 - 6*Ga3*beta1*ha0**3)/(ha0**2*(2304*beta1**3*ha0**3*ha2**3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
            ua0 =(2592*Ga0*beta1**3*ha0**3*ha2**3 + 2592*Ga0*beta1**3*ha0**2*ha1**2*ha2**2 - 2160*Ga0*beta1**3*ha0*ha1**4*ha2 + 1440*Ga0*beta1**3*ha1**6 - 792*Ga0*beta1**2*ha0**2*ha2**2 - 378*Ga0*beta1**2*ha0*ha1**2*ha2 - 192*Ga0*beta1**2*ha1**4 + 54*Ga0*beta1*ha0*ha2 + 27*Ga0*beta1*ha1**2 - Ga0 - 1296*Ga1*beta1**3*ha0**3*ha1*ha2**2 + 864*Ga1*beta1**3*ha0**2*ha1**3*ha2 - 1080*Ga1*beta1**3*ha0*ha1**5 + 108*Ga1*beta1**2*ha0**2*ha1*ha2 + 3*Ga1*beta1**2*ha0*ha1**3 - Ga1*beta1*ha0*ha1 - 288*Ga2*beta1**3*ha0**4*ha2**2 - 288*Ga2*beta1**3*ha0**3*ha1**2*ha2 + 720*Ga2*beta1**3*ha0**2*ha1**4 + 72*Ga2*beta1**2*ha0**3*ha2 + 66*Ga2*beta1**2*ha0**2*ha1**2 - 2*Ga2*beta1*ha0**2 + 432*Ga3*beta1**3*ha0**4*ha1*ha2 - 360*Ga3*beta1**3*ha0**3*ha1**3 - 60*Ga3*beta1**2*ha0**3*ha1)/(ha0*(2304*beta1**3*ha0**3*ha2**3 + 1008*beta1**3*ha0**2*ha1**2*ha2**2 - 576*beta1**3*ha0*ha1**4*ha2 + 360*beta1**3*ha1**6 - 720*beta1**2*ha0**2*ha2**2 - 204*beta1**2*ha0*ha1**2*ha2 - 189*beta1**2*ha1**4 + 52*beta1*ha0*ha2 + 26*beta1*ha1**2 - 1))
        elif (abs(ha1) < 10.0**(-15)):
            ua3 = (Ga0*ha1**3 - Ga1*ha0*ha1**2 + Ga2*ha0**2*ha1 - Ga3*ha0**3)/(ha0**4*(15*beta1*ha1**2 - 1))
            ua2 = (-36*Ga0*beta1*ha1**4 + Ga0*ha1**2 + 36*Ga1*beta1*ha0*ha1**3 - Ga1*ha0*ha1 - 36*Ga2*beta1*ha0**2*ha1**2 + Ga2*ha0**2 + 21*Ga3*beta1*ha0**3*ha1)/(ha0**3*(120*beta1**2*ha1**4 - 23*beta1*ha1**2 + 1))
            ua1 = (432*Ga0*beta1**2*ha1**5 - 27*Ga0*beta1*ha1**3 + Ga0*ha1 - 432*Ga1*beta1**2*ha0*ha1**4 + 27*Ga1*beta1*ha0*ha1**2 - Ga1*ha0 + 312*Ga2*beta1**2*ha0**2*ha1**3 - 4*Ga2*beta1*ha0**2*ha1 - 162*Ga3*beta1**2*ha0**3*ha1**2 - 6*Ga3*beta1*ha0**3)/(ha0**2*(360*beta1**3*ha1**6 - 189*beta1**2*ha1**4 + 26*beta1*ha1**2 - 1))
            ua0 = (1440*Ga0*beta1**3*ha1**6 - 192*Ga0*beta1**2*ha1**4 + 27*Ga0*beta1*ha1**2 - Ga0 - 1080*Ga1*beta1**3*ha0*ha1**5 + 3*Ga1*beta1**2*ha0*ha1**3 - Ga1*beta1*ha0*ha1 + 720*Ga2*beta1**3*ha0**2*ha1**4 + 66*Ga2*beta1**2*ha0**2*ha1**2 - 2*Ga2*beta1*ha0**2 - 360*Ga3*beta1**3*ha0**3*ha1**3 - 60*Ga3*beta1**2*ha0**3*ha1)/(ha0*(360*beta1**3*ha1**6 - 189*beta1**2*ha1**4 + 26*beta1*ha1**2 - 1))
        else:
            ua3 =Ga3/ha0
            ua2 =Ga2/ha0
            ua1 =Ga1/ha0 + 6*Ga3*beta1*ha0
            ua0 =Ga0/ha0 + 2*Ga2*beta1*ha0
            
    else:
            ua3 =Ga3/ha3
            ua2 =Ga2/ha2
            ua1 =Ga1/ha1
            ua0 =Ga0/ha0
                                                  
    return ua0,ua1,ua3,ua2



dx = 0.1
h0 = 2
G0 = 1

print('A')
A = SolveForuEdges(2.0/3.0,2,2,2,2,1,1,1,1,dx)   
print(A)
print('A')
A = SolveForuEdges(2.0/3.0,2,2,2,0,1,1,1,1,dx)  
print(A)
print('A')
A = SolveForuEdges(2.0/3.0,2,2,0,0,1,1,1,1,dx)  
print(A)
print('A')
A = SolveForuEdges(2.0/3.0,2,0,0,0,1,1,1,1,dx)  
print(A)
# print('A')
# A = SolveForuEdges(2.0/3.0,1,2,3,4,0,0,0,0,dx)   
# print('B')
# B = SolveForuEdges(2.0/3.0,1,2,3,4,1,0,0,0,dx)   
# print('C')
# C = SolveForuEdges(2.0/3.0,1,2,3,4,1,1,0,0,dx)   
# print('D')
# D = SolveForuEdges(2.0/3.0,1,2,3,4,1,1,1,0,dx)   
# print('E')
# E = SolveForuEdges(2.0/3.0,1,2,3,4,1,1,1,1,dx)   
#jpLets Check                                                     