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

def SympyToCode(Expr,u,h,xdepth):
    NewExpr = Expr
    
    for i in range(xdepth,-1,-1):
        if i == 0 :
            codestru = 'ui'
            codestrh = 'hi'
        elif i == 1 :
            codestru = 'dudx'
            codestrh = 'dhdx'
        else:
            codestru = 'd' + str(i) + 'udx' + str(i)
            codestrh = 'd' + str(i) + 'hdx' + str(i)
        NewExpr= NewExpr.subs(diff(u,(x,i)),codestru ).subs(diff(h,(x,i)),codestrh )
    return NewExpr

# symbols
x,t,ga = symbols('x t ga',real=True)

#Serre equations
h = Function('h')(x,t)
u = Function('u')(x,t)

tau = Function('tau')(x,t)

dhdt = diff(h,t)

fluxh = u*h
dfluxhdx = diff(fluxh,x)

LHSh = dhdt + dfluxhdx  

G = u*h - diff(h**3*diff(u,x),x)/3

dGdt = diff(G,t)
GReg = tau*(h*diff(h,(x,2)) - (diff(h,x)**2)/2)
dfluxGdx = diff(u*G + ga*h**2/2 - 2/3*h**3*diff(u,x)**2 - GReg,x)

LHSG = dGdt + dfluxGdx 


    
#Forced solutions
a0,a1,a2,a3,a4,a5,a6,a7 = symbols('a0,a1,a2,a3,a4,a5,a6,a7', positive = True, nonzero = True)

phi = x - a2*t
exp1 = exp(-(phi)**2 / (2*a3))
hforce = a0 + a1*exp1
dhfdx = simplify(diff(hforce,x).subs(exp1,'EXPPHI1')).subs(phi,'PHI')
d2hfdx2 = simplify(diff(hforce,(x,2)).subs(exp1,'EXPPHI1')).subs(phi,'PHI')
d3hfdx3 = simplify(diff(hforce,(x,3)).subs(exp1,'EXPPHI1')).subs(phi,'PHI')

uforce = a4*exp1
dufdx = simplify(diff(uforce,x).subs(exp1,'EXPPHI1')).subs(phi,'PHI')
d2ufdx2 = simplify(diff(uforce,(x,2)).subs(exp1,'EXPPHI1')).subs(phi,'PHI')
d3ufdx3 = simplify(diff(uforce,(x,3)).subs(exp1,'EXPPHI1')).subs(phi,'PHI')

tauforce = a6*(x - a5*t) + a7


#Condensed Form For Coding
Codeh = hforce.subs(exp1,'EXPPHI1').subs(phi,'PHI')
Codeu = uforce.subs(exp1,'EXPPHI1').subs(phi,'PHI')
CodeG = SympyToCode(G,u,h,3)

print('Code For Initial Conditions')
print('h(i) = '+ str(Codeh))
print('u(i) = '+ str(Codeu))
print('dhdx = '+ str(dhfdx))
print('dudx = '+ str(dufdx))
print('d2udx2 = '+ str(d2ufdx2))
print('G(i) = '+ str(CodeG))

print('\n\n\n')
print('Code For Forcing Terms')
Codedhdt = simplify((dhdt.subs(h,hforce).doit()).subs(exp1,'EXPPHI1') ).subs(phi,'PHI')
print('dhdt = ' + str(Codedhdt))
print()
CodedGdt = dGdt.subs(h,hforce).subs(u,uforce).subs(tau,tauforce).doit()
CodedGdt = simplify(CodedGdt.subs(uforce,'ui').subs(hforce,'hi') \
            .subs(tauforce,'tau').subs(exp1,'EXPPHI1').subs(phi,'PHI')).subs(a2*t - x,'(-PHI)')
CodedGdt = CodedGdt.expand()
print('dGdt = ' + str(CodedGdt))
print()


Codedfluxhdx = SympyToCode(dfluxhdx,u,h,3)
print('dfluxhdx = ' + str(Codedfluxhdx))
print()

CodedfluxGdx = SympyToCode(dfluxGdx,u,h,3)
print('dfluxGdx = ' + str(CodedfluxGdx))
print()
print('d2hdx2 = '+ str(d2hfdx2))
print('d3hdx3 = '+ str(d3hfdx3))
print('d3udx3 = '+ str(d3ufdx3))
print('\n\n\n')


