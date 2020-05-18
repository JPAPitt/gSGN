#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  7 14:22:34 2020

@author: jp
"""
from scipy import *
from numpy import *
import matplotlib.pyplot as plt

def GenSWWESerreP(k,g,u0,h0,a,e):
    wp =u0*k + k*sqrt(g*h0)*sqrt(1 + (a-1)*e*(k**2 *h0**2 )/ (e*k**2 *h0**2 + 1)) 
    return wp

def GenSWWESerrem(k,g,u0,h0,a,e):
    wm =u0*k - k*sqrt(g*h0)*sqrt(1 + (a-1)*e*(k**2 *h0**2 )/ (e*k**2 *h0**2 + 1)) 
    return wm

def EulerP(u0,h0,g,k):
    wp = u0*k + sqrt(g*k*tanh(k*h0))
    return wp
    
def EulerM(u0,h0,g,k):
    wm = u0*k - sqrt(g*k*tanh(k*h0))
    return wm


u0 = 0
h0 = 1
g = 9.81

k = arange(0,10,0.1)
Eulp = [EulerP(u0,h0,g,ki) for ki in k]
plt.plot(k,Eulp,label ='Linear Wave Speed')
#list of alpha, epsilon
Comb = [(0,1.0/3.0),(1,0),(0.5,1.0/3.0),(0.1,1.0/3.0)]
for a,e in Comb:

    EqP = [GenSWWESerreP(ki,g,u0,h0,a,e) for ki in k]
    str1 = 'Eq Wave Speed |' + 'alpha = ' + str(a) + '|' + 'epsilon = ' + str(round(e,3))
    plt.plot(k,EqP,label=str1)
    
plt.legend(loc='upper left')
