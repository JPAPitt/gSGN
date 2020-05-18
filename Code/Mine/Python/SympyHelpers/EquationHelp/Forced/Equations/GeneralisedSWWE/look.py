#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 12:49:56 2020

@author: jp
"""
from numpy import *
import matplotlib.pyplot as plt

def minmod(a,b,c):
    d = 0
    if a > 0 and b> 0 and c>0:
        d = min(a,b,c)
    elif a < 0 and b < 0 and c < 0:
        d = max(a,b,c)
    
    return d

ga = 9.81
a0 = 1.0
a1 = 0.5
a2 = 5.0
a3 = 20.0
a4 = 0.3

dx = 0.1
t = 0
x = arange(-100,100,dx)

phi  = x - a2*t
expphi1 = exp(-phi**2 / (2*a3))
h = a0 + a1*expphi1
dh = -expphi1*phi*a1/a3
ddh = -a1*(-phi**2 + a3)*expphi1/a3**2 


#plt.plot(x,h,label='h')
#plt.plot(x,dh,label='dh')
#plt.plot(x,ddh,label='ddh')
#plt.legend()
#plt.show()

xe = x + 0.5*dx
phie  = xe - a2*t
expphi1e = exp(-phie**2 / (2*a3))
he = a0 + a1*expphi1e
dhe = -expphi1e*phie*a1/a3
ddhe = -a1*(-phie**2 + a3)*expphi1e/a3**2 

#dhir = (2*hbc(i) -3*hbc(i-1) + hbc(i-2))  /dx
dhapproxeCent= (h[1:] - h[:-1])/dx
dhxCent = xe[0:-1]
dhapproxeBack = (2*h[2:-1] -3*h[1:-2] + h[0:-3])/dx
dhxBack = xe[2:-1]
#dhip1l = (-2*hbc(i+1) + 3*hbc(i+2) - hbc(i+3)) /dx
dhapproxeforward = (-2*h[1:-3] +3*h[2:-2] - h[3:-1])/dx
dhxforward = xe[0:-4]


x_com = xe[2:-4]
dhe_com = dhe[2:-4]
dhcent_com = dhapproxeCent[2:-3]
dhback_com = dhapproxeBack [0:-3]
dhforw_com = dhapproxeforward[2:]

n = len(x_com)
dhminmod_com = zeros(n)
for i in range(n):
    
    dhminmod_com[i] = minmod(dhback_com[i],dhcent_com[i],dhforw_com[i] )

plt.plot(x_com,dhe_com - dhcent_com ,label = 'centapprox')
plt.plot(x_com,dhe_com -dhback_com ,label = 'backapprox')
plt.plot(x_com,dhe_com-dhforw_com,label = 'forwardapprox')
plt.plot(x_com,dhe_com-dhminmod_com,label = 'limderivapprox')
#plt.plot(xe,dhe,label = 'actual')
plt.legend()