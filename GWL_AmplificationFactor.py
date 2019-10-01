# PROCEDURE: Amplification Factor 

# OUTPUT: plots of amplification factor (mod) vs. frequency (w) for
#        both point mass and isothermal systems

#%% Import libraries 
import matplotlib.pyplot as plt
import numpy as np
from itertools import islice
import operator
import math
import astropy
import scipy.stats as stats 
import math
from scipy.special import gamma, factorial, hyp1f1
import cmath
from mpmath import *

#%% Set parameters 
w = np.linspace( 0.001, 100, 1000) #frequency w = 8pi * M * f
y = 0.1    # source position
x = (y + np.sqrt((y**2)+4)/2)    # image position
phi = (((x-y)**2)/2)- math.log(x)   # phase


#%% Set funtions
def AmplificationFactor(y):
    x = (y + np.sqrt((y**2)+4)/2)    # image position
    phi = (((x-y)**2)/2)- math.log(x)   # phase
    def Gamma(x): #Gamma function
        r1 = np.real(1 - ((1j/2)*x))
        r2 = np.imag(1 - ((1j/2)*x))
        r = float(r1) + (float(r2)*1j)
        l = gamma(r) 
        return l

    L = list(map(Gamma,w))

    def Expo(x): #Exponential function 
        e1 = np.real((np.pi*x/4)+(1j*x/2)*(math.log(x/2)-(2*phi)))
        e2 = np.imag((np.pi*x/4)+(1j*x/2)*(math.log(x/2)-(2*phi)))
        e = float(e1)+float(e2)*1j
        return np.exp(e)
    Exponential = list(map(Expo,w))


    def F1(x):
        mp.dps = 5; mp.pretty = True
        i = 0 
        a = ((1j/2)*x)
        b = ((1j/2)*x)*(y**2)
        F_real = np.zeros(len(x))
        F_imag = np.zeros(len(x))
        while i < len(x):
            temp=hyp1f1(a[i], 1, b[i])
            F_real[i]= float(np.real(temp))
            F_imag[i]= float(np.imag(temp))
            i = i + 1
        return F_real + F_imag*1j 


    Exp_L = np.multiply(Exponential,L) #Expo * Gamma
    F = F1(w)   #Hyper Geometrical function
    Exp_L_F = np.multiply(Exp_L, F) #Result (Amplification Factor)
    return np.absolute(Exp_L_F)

#%%
plt.plot(w, AmplificationFactor(0.1), 'w')
plt.grid(True)
plt.xscale('log')
plt.yscale('log')
plt.ylim([0.1, 10])


#%%
