#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 13:43:44 2024

Plotting Normalised densities vs age of the

"""

# First let's set up our packages
import numpy as np
from matplotlib import pyplot as plt
from scipy import integrate
import math

# And set some constants
c = 299792.458 # km/s (speed of light)

H0kmsmpc = 70.  # Hubble constant in km/s/Mpc
H0s = H0kmsmpc * 3.2408e-20 # H0 in inverse seconds is H0 in km/s/Mpc * (3.2408e-20 Mpc/km)
H0y = H0s * 3.154e7 * 1.e9 # H0 in inverse Giga years is H0 in inverse seconds * (3.154e7 seconds/year) * (1e9 years / Giga year)


# Start by making an array of scalefactors
astart = 1e-4
astop = 10
astep = 0.0001 # Make this finer to make the plot smoother
a_arr = np.arange(astart,astop,astep)


def adotinv(a,om,ol,oR):
    adot= a*( om*(a)**(-3) + (1 - ol - om - oR)*a**(-2) + ol + oR*a**(-4) )**(1/2)
    return 1.0/adot

#functions for normalised densities
def norm_om(a,om,oR,ol):
    omfrac= om*a**(-3)/(om*a**(-3) + oR*a**(-4) + ol)
    return omfrac

def norm_ol(a,om,oR,ol):
    olfrac=  ol/(om*a**(-3) + oR*a**(-4) + ol)
    return olfrac

#set densities
om= 0.3
oR= 1e-3
ol= 0.6
ok = 1 - om - oR - ol

fracom = norm_om(a_arr,om,oR,ol)
fracol = norm_ol(a_arr,om,oR,ol)

# Note that when you integrate something with more than one argument you pass it with args=(arg1,arg2) in the integrate function
# e.g. "integrate.quad(adotinv, lower_limit, uper_limit, args=(om,ol))""
t3_lookback_Gyr = np.array([integrate.quad(adotinv, 0.0001, a_end, args=(om,ol,oR))[0] for a_end in a_arr])/H0y


plt.plot(t3_lookback_Gyr,a_arr,label='$(\Omega_M,\Omega_\Lambda,\Omega_R)$=(%.2f,%.2f,%.2f)'%(om,ol,oR)) 
plt.xlabel('Age (Gyr)')
plt.ylabel('Scalefactor')
plt.legend(loc='lower right',frameon=False)
plt.show()

#normalised mass density vs age (i think)
plt.plot(t3_lookback_Gyr,fracom,label='$(\Omega_M,\Omega_\Lambda,\Omega_R)$=(%.2f,%.2f,%.2f)'%(om,ol,oR)) 
plt.xlabel('Age of the universe (Gyr)')
plt.xscale('log') 
plt.ylabel('Normalised mass density')
plt.legend(loc='lower right',frameon=False)
plt.show()

#age vs normalised mass density (swapped the axes)
plt.plot(fracom,t3_lookback_Gyr,label='$(\Omega_M,\Omega_\Lambda,\Omega_R)$=(%.2f,%.2f,%.2f)'%(om,ol,oR)) 
plt.xlabel('Normalised mass density')
plt.yscale('log') 
plt.ylabel('Age of the universe (Gyr)')
plt.legend(loc='lower right',frameon=False)
plt.show()


#normalised cosmo constant vs age (i think)
plt.plot(t3_lookback_Gyr,fracol,label='$(\Omega_M,\Omega_\Lambda,\Omega_R)$=(%.2f,%.2f,%.2f)'%(om,ol,oR)) 
plt.xlabel('Age of the universe (Gyr)')
plt.xscale('log') 
plt.ylabel('Normalised Cosmological constant')
plt.legend(loc='lower right',frameon=False)
plt.show()

R0 = c/(H0s*math.sqrt(ok))
print(R0)




"""
om_arr = np.arange(0.001,4,0.5)

for om in om_arr:
    t4_lookback_Gyr = np.array([integrate.quad(adotinv, 1, a_end, args=(om,ol,oR))[0] for a_end in a_arr])/H0y
    plt.plot(t4_lookback_Gyr,a_arr,label='$(\Omega_M,\Omega_\Lambda)$=(%.1f,%.1f)'%(om,ol))
    
plt.axvline(x=0,linestyle=':') # Plot some crosshairs 
plt.axhline(y=1,linestyle=':')
plt.xlabel('Lookback time (Gyr)')
plt.ylabel('Scalefactor')
#plt.legend(loc='lower right',frameon=False)
plt.show()

"""





