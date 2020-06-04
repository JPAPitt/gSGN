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
    y = a0 + (a1 - a0)/2.0 *(1.0 + tanh(x/a2))
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

def ReconGradM(x,q,dx,theta):
    n = len(x)
    xn = []
    dqcn = []
    
    for i in range(2,n-3):
        xp = x[i] + 0.5*dx
        dqb = (2*q[i] - 3*q[i-1] + q[i-2])  /dx
        dqm = (q[i+1] - q[i]) /dx
        dqc = minmod((theta*dqb,dqm))
        
        xn.append(xp)
        dqcn.append(dqc)
    return array(xn),array(dqcn)


def ReconGradP(x,q,dx,theta):
    n = len(x)
    xn = []
    dqcn = []
    
    for i in range(2,n-3):
        xp = x[i] + 0.5*dx
        dqm = (q[i+1] - q[i]) /dx
        dqf = (-2*q[i+1] + 3*q[i+2] - q[i+3]) /dx
        dqc = minmod((theta*dqf,dqm))
        
        xn.append(xp)
        dqcn.append(dqc)
    return array(xn),array(dqcn)


def Recon2GradM(x,q,dx,theta):
    n = len(x)
    xn = []
    ddqcn = []
    
    for i in range(3,n-4):
        xp = x[i] + 0.5*dx
        
        ddqb = (5*q[i] - 13*q[i-1] + 11*q[i-2] - 3*q[i-3])/(2*dx**2)
        ddqbp = (3*q[i+1] -7*q[i] + 5*q[i-1] - q[i-2])/(2*dx**2)
        
        ddqm = (q[i+2]  - q[i+1] - q[i] + q[i-1] )/(2*dx**2)

        ddqc = minmod((theta*ddqb,theta*ddqbp,ddqm))
        
        xn.append(xp)        
        ddqcn.append(ddqc)
    return array(xn),array(ddqcn)        


def Recon2GradP(x,q,dx,theta):
    n = len(x)
    xn = []
    ddqcn = []
    
    for i in range(3,n-4):
        xp = x[i] + 0.5*dx
                
        ddqm = (q[i+2]  - q[i+1] - q[i] + q[i-1] )/(2*dx**2)
        
        ddqfp = (3*q[i] - 7*q[i+1] + 5*q[i+2] - q[i+3]) /(2*dx**2)
        ddqf = (5*q[i+1] - 13*q[i+2] + 11*q[i+3] - 3*q[i+4]) /(2*dx**2)

        ddqc = minmod((theta*ddqf,theta*ddqfp,ddqm))
        
        xn.append(xp)
        ddqcn.append(ddqc)
    return array(xn),array(ddqcn)          

def sech(x):
    return 2.0/ (exp(x) + exp(-x))

"""
q0 = 1
q1 = 0 
dbx = db(x,q0,q1)  

#Gradient with reconstruction - behaviour as expected, gradients 0 everywhere,
xn,db_dqbn,db_dqmn,db_dqfn,db_dqcn = ReconGrad(x,dbx,dx)

xn,db_ddqbn,db_ddqbpn,db_ddqmn,db_ddqfn,db_ddqfpn,db_ddqcn = Recon2Grad(x,dbx,dx)
"""    
        
a0 = 2
a1 = 1
a2 = 0.1


dxs = []
L1ms = []
L1ps = []
L12ms = []
L12ps = []

n = 5000
x = linspace(-50,50,n)  
dx = x[1] - x[0]
theta = 0.5

sdbx = smooth(x,a0,a1,a2)  

xdn,sdbx_dqcn_m = ReconGradM(x,sdbx,dx,theta)
xdn,sdbx_dqcn_p = ReconGradP(x,sdbx,dx,theta)

dsdbx =  - (0.5*(a0 - a1)/ a2)* sech(xdn/a2)**2

xddn,sdbx_ddqcn_m = Recon2GradM(x,sdbx,dx,theta)
xddn,sdbx_ddqcn_p = Recon2GradP(x,sdbx,dx,theta)

ddsdbx =   ((a0 - a1)/ a2**2)*tanh(xddn/a2)*sech(xddn/a2)**2





