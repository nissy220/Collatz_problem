"""
Created on Sat Dec  5 21:39:37 2020

@author: T.Nishiyama
"""
import numpy as np
import matplotlib.pyplot as plt

#Definitions of parameters
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
    hist=np.zeros(t_range)
    t_max=0
    for i in range(ns, ne+1):
        count=collatz(i)
        if count < t_range:
            hist[count]+=1
        if count > t_max:
            t_max=count
    return hist, t_max     

#Approximation of the scaled complementary error function
def erfcx(x):
    a=1.98
    b=1.135
    return (1-np.exp(-a*x))/(b * np.sqrt(np.pi)*x)

#Function that approximates the total stopping time distribution by the Brownian motion          
def log_phi(r, p, T):    
    x=np.log(r)+p*np.log(10)
    C=-s/np.sqrt(beta)-0.5*(3*s-v)/np.sqrt(alpha)
    D=-s/np.sqrt(beta)+0.5*(3*s-v)/np.sqrt(alpha)
    
    t=(x/sigma+2*s*T)/(3*s-v)
    X=np.sqrt(alpha*t)+T*np.sqrt(beta/t)
    Y=-np.sqrt(alpha*t)+T*np.sqrt(beta/t)
    Z=-alpha*t-beta*T**2/t
    g=(6*s**2-2*sigma*s)*T+Z
      
    return g + np.log(np.sqrt(2)*sigma*s) + np.log(C*erfcx(X)+D*erfcx(Y))   

#Estimate the total stopping time distribution for the numbers from ns=rs*10**ps to ne=re*10**pe
def estimate_distribution(rs, ps, re, pe, t_range):
    T=np.arange(t_range)
    ps = max(ps, 1)
    return np.exp(log_phi(re, pe, T)) - np.exp(log_phi(rs, ps, T))

#Estimate the max of total stopping time for the numbers from 1 to n=r*10**p
def estimate_max_stopping_time(r, p, t_range): 
    T=np.arange(t_range)
    log_count=log_phi(r, p, T)
    max_idx=np.argmax(log_count)
    min_idx=np.argmin(np.abs(log_count[max_idx:])) + max_idx 
    return min_idx    

#Inverse Gaussian distribution 
def estimate_density(r, p, t_range):
    T=np.arange(t_range)
    x=(np.log(r)+p*np.log(10))/sigma
    
    t=(x+2*s*T)/(3*s-v)
    return 2*s/(3*s-v)/np.sqrt(2*np.pi*t**3)*x*np.exp(-0.5*(x+v*t)**2/t)
    
def estimate_max_stopping_time2(r, p, t_range): 
    log_count=np.log(estimate_density(r, p, t_range)) + np.log(r) + p * np.log(10)
    max_idx=np.argmax(log_count)
    min_idx=np.argmin(np.abs(log_count[max_idx:])) + max_idx - 1
    return min_idx  

#range of the total stopping time
t_range=500
#range of numbers [rs*10^ps, re*10^pe]
rs=1
re=1
ps=0
pe=5

#Estimate the total stopping time distribution 
est_hist=estimate_distribution(rs, ps, re, pe, t_range) 
#est_hist2=estimate_density(re, pe, t_range) * re*10**pe

hist, t_max=calc_stopping_time(rs*10**ps, re*10**pe, t_range)
ave_range=3
#Smooth the histgram
filter=np.ones(ave_range) / ave_range
hist=np.convolve(hist,filter, mode='same')

plt.plot(np.arange(t_range), hist, label='Collatz sequences')
plt.plot(np.arange(t_range),est_hist, label='Brownian motion model')
#plt.plot(np.arange(t_range),est_hist2, label='inverse')

plt.xlabel('total stopping time')
plt.ylabel('counts')
plt.title('Total stopping time from' + "{:10.1e}".format(rs*10**ps) + ' to' + "{:10.1e}".format(re*10**pe))
plt.legend()
plt.show()

#Estimate the max of the total stopping time
t_range2=2500

#List of max of the total stopping time
collatz_t_max = np.array([0, 19, 118, 178, 261, 350, 524, 685, 949, 986, 1132, 1228,
                 1348, 1563, 1662, 1862, 1958, 2091, 2283])

##Estimate the max of the total stopping time
power_range=18
pred_t_max=np.zeros(power_range + 1)
for i in range(1, power_range + 1):
    pred_t_max[i]=estimate_max_stopping_time(1, i, t_range2)

plt.plot(np.arange(len(collatz_t_max)), collatz_t_max, 'o', label='Collatz sequences') 
plt.plot(np.arange(power_range + 1), pred_t_max, 'o', label='Brownian motion model')  
 
plt.xlabel('log_10(n)')   
plt.ylabel('Max of total stopping time')
plt.legend()
plt.show()

error=np.zeros(len(collatz_t_max))
error[1:]=np.abs(pred_t_max[1:len(collatz_t_max)]-collatz_t_max[1:]) / collatz_t_max[1:]
plt.plot(np.arange(2,len(collatz_t_max)), error[2:], 'o')             
plt.xlabel('log_10(n)')   
plt.ylabel('Error')            
    
         