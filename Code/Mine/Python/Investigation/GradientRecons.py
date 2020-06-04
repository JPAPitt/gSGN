#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 07:45:31 2020

@author: jp
"""
from scipy import *
from numpy import *
from numpy.linalg import norm
from matplotlib.pyplot import plot

def smooth(x,a0,a1,a2):
    y = a0 + a1*sin(a2*x)
    return y
    
def db(x,q0,q1):
    n = len(x)
    y = zeros(n)
    for i in range(0,n):
        if x[i] < 0:
            y[i] = q0
        else:
            y[i] = q1
    return y

def minmod(a):
    n = len(a)
    sgna = sign(a)
    if (1 in sgna) and (-1 not in sgna):
        return min(a)
    elif (1 not in sgna) and (-1 in sgna):
        return max(a)
    else:
        return 0
    
def minmod2(a):
    n = len(a)
    mina = min(a)
    maxa = max(a)
    minabsa = min(abs(a))
    if (mina > 0 and minabsa > 0):
        return min(a)
    elif (maxa < 0 and minabsa > 0):
        return max(a)
    else:
        return 0

def ReconGradM(x,q,dx):
    n = len(x)
    xn = []
    dqcn = []
    
    for i in range(2,n-3):
        xp = x[i] + 0.5*dx
        dqb = (2*q[i] - 3*q[i-1] + q[i-2])  /dx
        dqm = (q[i+1] - q[i]) /dx
        dqc = minmod((dqb,dqm))
        
        xn.append(xp)
        dqcn.append(dqc)
    return array(xn),array(dqcn)


def ReconGradP(x,q,dx):
    n = len(x)
    xn = []
    dqcn = []
    
    for i in range(2,n-3):
        xp = x[i] + 0.5*dx
        dqm = (q[i+1] - q[i]) /dx
        dqf = (-2*q[i+1] + 3*q[i+2] - q[i+3]) /dx
        dqc = minmod((dqf,dqm))
        
        xn.append(xp)
        dqcn.append(dqc)
    return array(xn),array(dqcn)


def Recon2GradM(x,q,dx):
    n = len(x)
    xn = []
    ddqcn = []
    
    for i in range(3,n-4):
        xp = x[i] + 0.5*dx
        
        ddqb = (5*q[i] - 13*q[i-1] + 11*q[i-2] - 3*q[i-3])/(2*dx**2)
        ddqbp = (3*q[i+1] -7*q[i] + 5*q[i-1] - q[i-2])/(2*dx**2)
        
        ddqm = (q[i+2]  - q[i+1] - q[i] + q[i-1] )/(2*dx**2)

        ddqc = minmod((ddqb,ddqbp,ddqm))
        
        xn.append(xp)        
        ddqcn.append(ddqc)
    return array(xn),array(ddqcn)        


def Recon2GradP(x,q,dx):
    n = len(x)
    xn = []
    ddqcn = []
    
    for i in range(3,n-4):
        xp = x[i] + 0.5*dx
                
        ddqm = (q[i+2]  - q[i+1] - q[i] + q[i-1] )/(2*dx**2)
        
        ddqfp = (3*q[i] - 7*q[i+1] + 5*q[i+2] - q[i+3]) /(2*dx**2)
        ddqf = (5*q[i+1] - 13*q[i+2] + 11*q[i+3] - 3*q[i+4]) /(2*dx**2)

        ddqc = minmod((ddqf,ddqfp,ddqm))
        
        xn.append(xp)
        ddqcn.append(ddqc)
    return array(xn),array(ddqcn)          

"""
q0 = 1
q1 = 0 
dbx = db(x,q0,q1)  

#Gradient with reconstruction - behaviour as expected, gradients 0 everywhere,
xn,db_dqbn,db_dqmn,db_dqfn,db_dqcn = ReconGrad(x,dbx,dx)

xn,db_ddqbn,db_ddqbpn,db_ddqmn,db_ddqfn,db_ddqfpn,db_ddqcn = Recon2Grad(x,dbx,dx)
"""    
        
a0 = 1
a1 = 0.5
a2 = pi/10

n = 5000

dxs = []
L1ms = []
L1ps = []
L12ms = []
L12ps = []

for i in range(5,6):
    n = 100*2**(i)
    x = linspace(-50,50,n)  
    dx = x[1] - x[0]
    
    sinx = smooth(x,a0,a1,a2)  
    
    xn,sin_dqcn_m = ReconGradM(x,sinx,dx)
    xn,sin_dqcn_p = ReconGradP(x,sinx,dx)
    dsinxn = a2*a1*cos(a2*xn)
    
    L1m = norm(dsinxn - sin_dqcn_m,ord=2)/ norm(dsinxn,ord=2)
    L1p = norm(dsinxn - sin_dqcn_p,ord=2)/ norm(dsinxn,ord=2)
    
    xn,sin_ddqcn_m = Recon2GradM(x,sinx,dx)
    xn,sin_ddqcn_p = Recon2GradP(x,sinx,dx)
    ddsinxn = -(a2**2)*a1*sin(a2*xn)
    
    L12m = norm(ddsinxn - sin_ddqcn_m,ord=2)/ norm(ddsinxn,ord=2)
    L12p = norm(ddsinxn - sin_ddqcn_p,ord=2)/ norm(ddsinxn,ord=2)
    
    dxs.append(dx)
    L1ms.append(L1m)
    L1ps.append(L1p)
    L12ms.append(L12m)
    L12ps.append(L12p)


