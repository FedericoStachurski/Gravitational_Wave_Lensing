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
path="Plots/"

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
    return Exp_L_F

def PhaseFactor(y):
    mp.dps = 15; mp.pretty = True
    ratio = np.divide(AmplificationFactor(y), np.absolute(AmplificationFactor(y)))
    i = 0
    r_real = np.zeros(len(ratio))
    r_imag = np.zeros(len(ratio))
    while i < len(ratio):
        r_real[i] = float(np.real(ratio[i]))
        r_imag[i] = float(np.imag(ratio[i]))
        i += 1
    #print(r_real,type(r_real))
    theta = np.arctan(np.divide(r_imag,r_real))
    theta = (-1j)*np.log(np.exp(1j*theta))
    return theta


#%% Plot and save figure (Amplificatio factor PML)
fig = plt.figure()
plt.style.use('default')
ax = plt.subplot(111)
#ax.plot(w, np.absolute(AmplificationFactor(0.01)),'c',label='$y = 0.01$')
#ax.plot(w, np.absolute(AmplificationFactor(0.05)),'b' ,label='$y = 0.05$')
ax.plot(w, np.absolute(AmplificationFactor(0.1)), 'k',label='$y = 0.1$')
ax.plot(w, np.absolute(AmplificationFactor(0.3)), 'y',label='$y = 0.3$')
ax.plot(w, np.absolute(AmplificationFactor(1.0)), 'r',label='$y = 1.0$')
ax.plot(w, np.absolute(AmplificationFactor(3.0)), 'g',label='$y = 3.0$')

ax.legend()
matplotlib.pyplot.grid(True, which="both", alpha=0.2)
plt.title('Point Mass Lens')
plt.ylabel('$|F|$')
plt.xlabel('$w ( = 8\pi M_{z} f )$')
plt.xscale('log')
plt.yscale('log')
plt.ylim([0.1, 100])
plt.xlim([ 0.001, 100])

fig.savefig(path+'AmplFactor_PM.png', dpi=1200, bbox_inches = "tight")

#%% Plot and save figure (Phase PML)
fig = plt.figure()
plt.style.use('default')
ax = plt.subplot(111)
#ax.plot(w, PhaseFactor(0.01),'c',label='$y = 0.01$')
#ax.plot(w, PhaseFactor(0.05),'b' ,label='$y = 0.05$')
ax.plot(w, PhaseFactor(0.1), 'k',label='$y = 0.1$')
ax.plot(w, PhaseFactor(0.3), 'y',label='$y = 0.3$')
ax.plot(w, PhaseFactor(1.0), 'r',label='$y = 1.0$')
ax.plot(w, PhaseFactor(3.0), 'g',label='$y = 3.0$', alpha=0.7)

ax.legend()
matplotlib.pyplot.grid(True, which="both", alpha=0.2)
plt.title('Point Mass Lens')
plt.ylabel(r"$ {\theta}_{F} $")
plt.xlabel('$w ( = 8\pi M_{z} f )$')
plt.xscale('log')
plt.ylim([-1.6, 1.6])
plt.xlim([ 0.001, 100])

fig.savefig(path+'AmplFactorPhase_PM.png', dpi=1200, bbox_inches = "tight")


#%%
