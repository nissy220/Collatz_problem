"""
Created on Sat Dec  5 21:39:37 2020

@author: T.Nishiyama
"""
import numpy as np
import matplotlib.pyplot as plt

#Definition of parameters
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

def calc_stopping_time(ns, ne, t_range):
    hist=np.zeros(t_range - 1)
    t_max=0
    for i in range(ns, ne):
        count=collatz(i)
        if count < t_range-1:
            hist[count]+=1
        if count > t_max:
            t_max=count
    return hist, t_max     

#Approximation of the erfcx function
def erfcx(x):
    a=1.98
    b=1.135
    return (1-np.exp(-a*x))/(b * np.sqrt(np.pi)*x)

#Function that approximates the stopping time distribution by the Brownian motion model           
def phi(n, T):
    C=np.sqrt(2)*sigma*s*(-s/np.sqrt(beta)-0.5*(3*s-v)/np.sqrt(alpha))
    D=-np.sqrt(2)*sigma*s*(s/np.sqrt(beta)-0.5*(3*s-v)/np.sqrt(alpha))
    
    t=(np.log(n)/sigma+2*s*T)/(3*s-v)
    x=np.sqrt(alpha*t)+T*np.sqrt(beta/t)
    y=-np.sqrt(alpha*t)+T*np.sqrt(beta/t)
    z=-alpha*t-beta*T**2/t
    g=(6*s**2-2*sigma*s)*T+z
      
    return (C*erfcx(x)+D*erfcx(y))*np.exp(g)   

#Estimate the stopping time distribution for the numbers from ns to ne
def estimate_distribution(ns, ne, t_range):
    T=np.arange(1, t_range)
    return phi(ne, T) - phi(ns, T)

#Estimate the max of stopping time for the numbers from 1 to n
def estimate_max_stopping_time(n, t_range): 
    total_hist=estimate_distribution(1, n, t_range ) 
    return np.argmin(np.abs(total_hist[30:]-1))    

#Inverse Gaussian distribution 
#def dist(n, T):
#    x=np.log(n)/sigma
#    t=(x+2*s*T)/(3*s-v)
#    return 2*s/(3*s-v)/np.sqrt(2*np.pi*t**3)*x*np.exp(-0.5*(x+v*t)**2/t)
#    
#

#range of the stopping time
t_range=500

ns=1
ne=10**6

#Calculate and estimate the stopping time distribution 
total_hist=estimate_distribution(ns, ne, t_range) 

hist, t_max=calc_stopping_time(ns, ne, t_range)
ave_range=5
filter=np.ones(ave_range) / ave_range
hist=np.convolve(hist,filter, mode='same')
plt.plot(np.arange(1, t_range), hist, label='Collatz sequence ')
plt.plot(np.arange(1, t_range),total_hist, label='Brownian motion model')
plt.xlabel('stopping time')
plt.ylabel('count')
plt.title('Stopping time distribution from' + "{:10.1e}".format(ns) + ' to' + "{:10.1e}".format(ne))
plt.legend()
plt.show()

t_range=3000

#Estimate the max of the stopping time
for i in range(4, 18):
    n=10**i
    pred_t_max=estimate_max_stopping_time(n, t_range)
#    print('max of stopping time =', t_max)
    print('less than'+ "{:10.1e}".format(n) + ', prediction of max of stopping time =', pred_t_max, 'steps')


    
            