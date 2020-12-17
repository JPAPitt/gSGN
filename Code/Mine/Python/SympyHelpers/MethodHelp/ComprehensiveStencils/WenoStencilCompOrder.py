from sympy import *
from sympy.solvers.solveset import linsolve


Vimh,Viph,Vim3h,Vip3h,dx = symbols('V_{i-1/2} V_{i+1/2} V_{i-3/2} V_{i+3/2} dx')
x,xiph = symbols('x x_{i+1/2}')

#Poly0 at x^{-}_j+1/2
P0iph = Viph


#Poly1 at x^{-}_j+1/2
P1iphtoip3h = (Vip3h - Viph)/dx*(x - xiph) + Viph
P1imhtoiph = (Viph - Vimh)/dx*(x - xiph) + Viph

#Poly2 at x^{-}_j+1/2