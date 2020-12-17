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
    
h = Function('h')(x,t)
u = Function('u')(x,t)

# #Symbols

Gaj = symbols('Gaj',real=True)

Gjph,Gjps,Gjms,Gjmh,SGj,hjph,hjps,hjms,hjmh,Shj = symbols('Gjph Gjps Gjms Gjmh Gj hjph hjps hjms hjmh hj',real=True)
ujph,ujps,ujms,ujmh,Suj = symbols('ujph ujps ujms ujmh uj',real=True)

Sdhjph,Sdhjps,Sdhjms,Sdhjmh,Sdhj =  symbols('dhjph dhjps dhjms dhjmh dhj',real =True)
Sdujph,Sdujps,Sdujms,Sdujmh,Sduj  =  symbols('u^d_{j+1/2} u^d_{j+1/6} u^d_{j-1/6} u^d_{j-1/2} u^{d}_j',real =True)
Sddujph,Sddujps,Sddujms,Sddujmh,Sdduj  =  symbols('u^{dd}_{j+1/2} u^{dd}_{j+1/6} u^{dd}_{j-1/6} u^{dd}_{j-1/2} u^{dd}_j',real =True)

ha0,ha1,ha2,ha3 = PolyFromPoints(hjmh,hjms,hjps,hjph,dx)
ua0,ua1,ua2,ua3 = PolyFromPoints(ujmh,ujms,ujps,ujph,dx) 
Ga0,Ga1,Ga2,Ga3 = PolyFromPoints(Gjmh,Gjms,Gjps,Gjph,dx)

hPoly = ha3*(x)**3 + ha2*(x)**2 + ha1*(x) + ha0
uPoly = ua3*(x)**3 + ua2*(x)**2 + ua1*(x) + ua0
GPoly = Ga3*(x)**3 + Ga2*(x)**2 + Ga1*(x) + Ga0

uj = ua0
hj = ha0
Gj = Ga0

duj = ua1
dhj = ha1
dGj = Ga1

dduj = 2*ua2
ddhj = 2*ha2
ddGj = 2*Ga2

Gajint = integrate(GPoly,(x,-dx/2,dx/2)).simplify()

dujph= diff(uPoly,x).doit().subs(x,dx/2).simplify()
dujmh= diff(uPoly,x).doit().subs(x,-dx/2).simplify()

dhjph= diff(hPoly,x).doit().subs(x,dx/2).simplify()
dhjmh= diff(hPoly,x).doit().subs(x,-dx/2).simplify()

ddujph= diff(uPoly,(x,2)).doit().subs(x,dx/2).simplify()
ddujmh= diff(uPoly,(x,2)).doit().subs(x,-dx/2).simplify()

Geq = u*h - b1*diff(h**3*diff(u,x), x)

#Simpsons 3/8 rule to integrate uh
GeqInt = dx/8*(hjmh*ujmh + 3*hjms*ujms + 3*hjps*ujps + hjph*ujph) - (b1*hjph**3*dujph - b1*hjmh**3*dujmh)
GeqInt =GeqInt.expand()
Geq =Geq.doit()

Geqjmh = Geq.subs(diff(u,(x,2)),ddujmh ).subs(diff(u,(x)),dujmh ).subs(u,ujmh ).subs(diff(h,(x)),Sdhjmh ).subs(h,hjmh )
Geqjmh =Geqjmh.expand()

Geqjph = Geq.subs(diff(u,(x,2)),ddujph ).subs(diff(u,(x)),dujph ).subs(u,ujph ).subs(diff(h,(x)),Sdhjph ).subs(h,hjph )
Geqjph =Geqjph.expand()

Geqj = Geq.subs(diff(u,(x,2)),dduj ).subs(diff(u,(x)),duj ).subs(u,uj ).subs(diff(h,(x)),Sdhj ).subs(h,Shj )
Geqj =Geqj.expand()


Eq1 =GeqInt/dx - Gaj
Eq1 = Eq1.expand()
Eq2 =Geqjmh - Gjmh
Eq3 = Geqjph - Gjph
Eq4 = Geqj - SGj

Eq1ujmh = collect(Eq1,ujmh)
Eq1ujms = collect(Eq1ujmh,ujms)
Eq1ujps = collect(Eq1ujms,ujps)
Eq1ujph = collect(Eq1ujps,ujph)


Eq1ujmhC = Eq1ujph.subs(Gaj,0).subs(ujmh,1).subs(ujms,0).subs(ujps,0).subs(ujph,0)
Eq1ujmsC = Eq1ujph.subs(Gaj,0).subs(ujmh,0).subs(ujms,1).subs(ujps,0).subs(ujph,0)
Eq1ujpsC = Eq1ujph.subs(Gaj,0).subs(ujmh,0).subs(ujms,0).subs(ujps,1).subs(ujph,0)
Eq1ujphC = Eq1ujph.subs(Gaj,0).subs(ujmh,0).subs(ujms,0).subs(ujps,0).subs(ujph,1)


a0,a1,a2,a3 =  symbols('a0 a1 a2 a3',real =True)
b0,b1,b2,b3 =  symbols('b0 b1 b2 b3',real =True)
c0,c1,c2,c3 =  symbols('c0 c1 c2 c3',real =True)
d0,d1,d2,d3 =  symbols('d0 d1 d2 d3',real =True)

Eq1R = Eq1ujph.subs(Eq1ujmhC,a0).subs(Eq1ujmsC,a1).subs(Eq1ujpsC,a2).subs(Eq1ujphC,a3)

Eq2ujmh = collect(Eq2,ujmh)
Eq2ujms = collect(Eq2ujmh,ujms)
Eq2ujps = collect(Eq2ujms,ujps)
Eq2ujph = collect(Eq2ujps,ujph)

Eq2ujmhC = Eq2ujph.subs(Gjmh,0).subs(ujmh,1).subs(ujms,0).subs(ujps,0).subs(ujph,0)
Eq2ujmsC = Eq2ujph.subs(Gjmh,0).subs(ujmh,0).subs(ujms,1).subs(ujps,0).subs(ujph,0)
Eq2ujpsC = Eq2ujph.subs(Gjmh,0).subs(ujmh,0).subs(ujms,0).subs(ujps,1).subs(ujph,0)
Eq2ujphC = Eq2ujph.subs(Gjmh,0).subs(ujmh,0).subs(ujms,0).subs(ujps,0).subs(ujph,1)

Eq2R = Eq2ujph.subs(Eq2ujmhC,b0).subs(Eq2ujmsC,b1).subs(Eq2ujpsC,b2).subs(Eq2ujphC,b3)

Eq3ujmh = collect(Eq3,ujmh)
Eq3ujms = collect(Eq3ujmh,ujms)
Eq3ujps = collect(Eq3ujms,ujps)
Eq3ujph = collect(Eq3ujps,ujph)

Eq3ujmhC = Eq3ujph.subs(Gjph,0).subs(ujmh,1).subs(ujms,0).subs(ujps,0).subs(ujph,0)
Eq3ujmsC = Eq3ujph.subs(Gjph,0).subs(ujmh,0).subs(ujms,1).subs(ujps,0).subs(ujph,0)
Eq3ujpsC = Eq3ujph.subs(Gjph,0).subs(ujmh,0).subs(ujms,0).subs(ujps,1).subs(ujph,0)
Eq3ujphC = Eq3ujph.subs(Gjph,0).subs(ujmh,0).subs(ujms,0).subs(ujps,0).subs(ujph,1)

Eq3R = Eq3ujph.subs(Eq3ujmhC,c0).subs(Eq3ujmsC,c1).subs(Eq3ujpsC,c2).subs(Eq3ujphC,c3)

Eq4ujmh = collect(Eq4,ujmh)
Eq4ujms = collect(Eq4ujmh,ujms)
Eq4ujps = collect(Eq4ujms,ujps)
Eq4ujph = collect(Eq4ujps,ujph)

Eq4ujmhC = Eq4ujph.subs(SGj,0).subs(ujmh,1).subs(ujms,0).subs(ujps,0).subs(ujph,0)
Eq4ujmsC = Eq4ujph.subs(SGj,0).subs(ujmh,0).subs(ujms,1).subs(ujps,0).subs(ujph,0)
Eq4ujpsC = Eq4ujph.subs(SGj,0).subs(ujmh,0).subs(ujms,0).subs(ujps,1).subs(ujph,0)
Eq4ujphC = Eq4ujph.subs(SGj,0).subs(ujmh,0).subs(ujms,0).subs(ujps,0).subs(ujph,1)

Eq4R = Eq4ujph.subs(Eq4ujmhC,d0).subs(Eq4ujmsC,d1).subs(Eq4ujpsC,d2).subs(Eq4ujphC,d3)

# #Most Expensive step
# uPoints = linsolve([Eq1R,Eq2R,Eq3R, Eq4R], (ujph,ujps,ujms,ujmh))

# uSolution = uPoints.args[0]
# ujphSol = uSolution[0]
# ujpsSol = uSolution[1]
# ujmsSol= uSolution[2]
# ujmhSol = uSolution[3]

# ujmhSol = ujmhSol.subs(a0,Eq1ujmhC).subs(a1,Eq1ujmsC).subs(a2,Eq1ujpsC).subs(a3,Eq1ujphC) \
#     .subs(b0,Eq2ujmhC).subs(b1,Eq2ujmsC).subs(b2,Eq2ujpsC).subs(b3,Eq2ujphC) \
#     .subs(c0,Eq3ujmhC).subs(c1,Eq3ujmsC).subs(c2,Eq3ujpsC).subs(c3,Eq3ujphC)\
#     .subs(d0,Eq4ujmhC).subs(d1,Eq4ujmsC).subs(d2,Eq4ujpsC).subs(d3,Eq4ujphC)



# ujmsSol=ujmsSol.subs(a0,Eq1ujmhC).subs(a1,Eq1ujmsC).subs(a2,Eq1ujpsC).subs(a3,Eq1ujphC) \
#     .subs(b0,Eq2ujmhC).subs(b1,Eq2ujmsC).subs(b2,Eq2ujpsC).subs(b3,Eq2ujphC) \
#     .subs(c0,Eq3ujmhC).subs(c1,Eq3ujmsC).subs(c2,Eq3ujpsC).subs(c3,Eq3ujphC)\
#     .subs(d0,Eq4ujmhC).subs(d1,Eq4ujmsC).subs(d2,Eq4ujpsC).subs(d3,Eq4ujphC)

# ujpsSol=ujpsSol.subs(a0,Eq1ujmhC).subs(a1,Eq1ujmsC).subs(a2,Eq1ujpsC).subs(a3,Eq1ujphC) \
#     .subs(b0,Eq2ujmhC).subs(b1,Eq2ujmsC).subs(b2,Eq2ujpsC).subs(b3,Eq2ujphC) \
#     .subs(c0,Eq3ujmhC).subs(c1,Eq3ujmsC).subs(c2,Eq3ujpsC).subs(c3,Eq3ujphC)\
#     .subs(d0,Eq4ujmhC).subs(d1,Eq4ujmsC).subs(d2,Eq4ujpsC).subs(d3,Eq4ujphC)

# ujpSolh=ujphSol.subs(a0,Eq1ujmhC).subs(a1,Eq1ujmsC).subs(a2,Eq1ujpsC).subs(a3,Eq1ujphC) \
#     .subs(b0,Eq2ujmhC).subs(b1,Eq2ujmsC).subs(b2,Eq2ujpsC).subs(b3,Eq2ujphC) \
#     .subs(c0,Eq3ujmhC).subs(c1,Eq3ujmsC).subs(c2,Eq3ujpsC).subs(c3,Eq3ujphC) \
#     .subs(d0,Eq4ujmhC).subs(d1,Eq4ujmsC).subs(d2,Eq4ujpsC).subs(d3,Eq4ujphC)

# ujmhSol = ujmhSol.simplify()   
# ujmsSol = ujmsSol.simplify()   
# ujpsSol = ujpsSol.simplify()   
# ujphSol = ujphSol.simplify()   
    
