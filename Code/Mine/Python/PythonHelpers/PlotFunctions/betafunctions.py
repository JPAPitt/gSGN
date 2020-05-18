#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 14:22:34 2020

@author: jp
"""
from scipy import *
from numpy import *
import matplotlib.pyplot as plt

def betafunc(a5,a6,a7,x,t):
    return a6*(x - a5*t) + a7


x0 = -50
x1 = 100
t0 = 0
t1 = 10
b0 = -2.0/3.0
b1 = 1
b2b1 =2
a5 = -0.9
phi0 = x0 - a5*t0
phi1 = x1 - a5*t1

b1a6 = (b1 - b0)/ (phi1 - phi0)
b1a7 = b0 - b1a6*phi0

b2a6 = (b2b1)/ (phi1 - phi0)
b2a7 = - b2a6*phi0


x = arange(x0,x1,0.1)
Beta1t0 = [betafunc(a5,b1a6,b1a7,xi,t0) for xi in x]
Beta1t1 = [betafunc(a5,b1a6,b1a7,xi,t1) for xi in x]

Beta2t0 = [betafunc(a5,b2a6,b2a7,xi,t0) for xi in x]
Beta2t1 = [betafunc(a5,b2a6,b2a7,xi,t1) for xi in x]
plt.plot(x,Beta1t0,label ='beta1 t = t0')
plt.plot(x,Beta1t1,label ='beta1 t = t1')

plt.plot(x,Beta2t0,label ='beta2 t = t0')
plt.plot(x,Beta2t1,label ='beta2 t = t1')

plt.plot(x,array(Beta2t1)/(2.0/3.0 + array(Beta1t1)),label ='ratio of beta2/beta1')
    
plt.legend(loc='upper left')
