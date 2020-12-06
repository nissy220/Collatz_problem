# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 21:39:37 2020

@author: Nishiyama
"""

import numpy as np
import matplotlib.pyplot as plt

sigma=0.5*np.log(3)
v=0.5*np.log(0.75)/sigma
s=0.5*np.log(3)/sigma
alpha=4.5*s**2-sigma*(3*s-v)
beta=2*s**2

def collatz(n):
    count=0
    while (n > 1):
        if (n % 2 == 0):
            n = n // 2
        else:
            n = 3 * n + 1
        count+=1   
    return count

def calc_stopping_time(ns, ne):
    t_range=400
    hist=np.zeros((t_range, 1))
    for i in range(ns, ne):
        count=collatz(i)
        if count < t_range:
            hist[count]+=1
    return hist     

def calc_stopping_time_2d(ns, ne, step1, step2):
    
    t_range=300
    n_range=10**4
    hist=np.zeros((n_range // step2, t_range // step1))
    for i in range(1, n_range):
        count=collatz(i)
        if count < t_range:
            hist[i// step2, count // step1]+=1
    return hist   
#approximation of the erfc function
def erfcx(x):
    a=1.98
    b=1.135
    return (1-np.exp(-a*x))/(b * np.sqrt(np.pi)*x)
            
def integral_dist(n, T):
    C=np.sqrt(2)*sigma*s*(-s/np.sqrt(beta)-0.5*(3*s-v)/np.sqrt(alpha))
    D=-np.sqrt(2)*sigma*s*(s/np.sqrt(beta)-0.5*(3*s-v)/np.sqrt(alpha))
    
    t=(np.log(n)/sigma+2*s*T)/(3*s-v)
    x=np.sqrt(alpha*t)+T*np.sqrt(beta/t)
    y=-np.sqrt(alpha*t)+T*np.sqrt(beta/t)
    z=-alpha*t-beta*T**2/t
    g=(6*s**2-2*sigma*s)*T+z
      
    return (C*erfcx(x)+D*erfcx(y))*np.exp(g)   

def total_dist(ns, ne, T):
    return integral_dist(ne, T) - integral_dist(ns, T)

def calc_max_stopping_time(n):
    t=np.arange(1,10**4)    
    total_hist=total_dist(1, n, t ) 
    return np.argmin(np.abs(total_hist[30:]-1))    
    
def dist(n, T):
    x=np.log(n)/sigma
    t=(x+2*s*T)/(3*s-v)
    return 2*s/(3*s-v)/np.sqrt(2*np.pi*t**3)*x*np.exp(-0.5*(x+v*t)**2/t)
    
step=10
step2=100
t=np.arange(1,400, step)
#n=np.arange(1, 10**4, step2)
#T, N= np.meshgrid(t, n)
#
total_hist=total_dist(4*10**5, 5*10**5, t) 
#total_hist=total_dist(N, N+100, T) 
#total_hist[total_hist>1] = 1
#total_hist[total_hist<1] = 0
#hist_2d=calc_stopping_time_2d(1, 10**4, step, step2)
#plt.pcolormesh(N, T, total_hist, cmap='jet')          
#plt.pcolormesh(N, T, hist_2d, cmap='jet')              
#          
#hist=dist(N, T)
hist=calc_stopping_time(4*10**5, 5*10**5)
ave_range=5
v=np.ones(ave_range) / ave_range
hist=np.convolve(hist.reshape(hist.shape[0]),v, mode='same')
plt.plot(hist)
plt.plot(t, total_hist)



    
            