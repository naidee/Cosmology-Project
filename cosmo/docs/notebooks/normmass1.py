#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 16:54:15 2024

@author: sakhi
"""

# First let's set up our packages
import numpy as np
from matplotlib import pyplot as plt
from scipy import integrate

# And set some constants
c = 299792.458 # km/s (speed of light)

H0kmsmpc = 70.  # Hubble constant in km/s/Mpc
H0s = H0kmsmpc * 3.2408e-20 # H0 in inverse seconds is H0 in km/s/Mpc * (3.2408e-20 Mpc/km)
H0y = H0s * 3.154e7 * 1.e9 # H0 in inverse Giga years is H0 in inverse seconds * (3.154e7 seconds/year) * (1e9 years / Giga year)


# Start by making an array of scalefactors
astart = 1e-4
astop = 5
astep = 0.0001 # Make this finer to make the plot smoother
a_arr = np.arange(astart,astop,astep)

# First set up an array of times (initially all set to zero) into which we'll put our calculated times
t_Gyr = np.zeros(len(a_arr))  # len(a_arr) gives the length of the a_arr 


# First calculate the index corresponding to a=1.0.  (Find when |a-1.0| is minimum.  You could also do this by just redoing the integral from 0<a<1, but the way I've set it up above we know we have an a=1 in the array, so we can just find what we've already calculated.)
index_today = np.argmin(np.abs(a_arr - 1.0))


# First write a function that takes as input a, Omega_M (om), and Omega_Lambda (ol) and outputs 1/adot
def adotinv(a,om,ol):
    ok = 1.0-om-ol
    adot= a*( om*(a)**(-3) + ok*a**(-2) + ol )**(1/2)
    return 1.0/adot

#normalised density functions
def norm_om(a,om,oR,ol):
    omfrac= om*a**(-3)/(om*a**(-3) + oR*a**(-4) + ol)
    return omfrac

def norm_oR(a,om,oR,ol):
    oRfrac= oR*a**(-4)/(om*a**(-3) + oR*a**(-4) + ol)
    return oRfrac

def norm_ol(a,om,oR,ol):
    olfrac=  ol/(om*a**(-3) + oR*a**(-4) + ol)
    return olfrac


#set densities
#om= 0.3
#oR= 1e-3
#ol= 0.7

#trial densities
om = 0.175
oR = 0.05
ol = 0.8


fracom = norm_om(a_arr,om,oR,ol)
fracoR = norm_oR(a_arr,om,oR,ol)
fracol = norm_ol(a_arr,om,oR,ol)

#plt.semilogx
#plt.loglog
#plt.semilogy
plt.plot(a_arr, fracoR,'r',label="oR")
plt.plot(a_arr,fracol,'b',label="ol")
plt.plot(a_arr,fracom,'g',label="om")
plt.xlabel('Scalefactor')
plt.ylabel('Normalised density')
plt.xscale('log') 
plt.legend()
plt.show()

"""""
Plot in look back time:

for a in a_arr:
    t3_lookback_Gyr = np.array([integrate.quad(adotinv, 1, a_end, args=(om,ol))[0] for a_end in a_arr])/H0y
    fracom = norm_om(a_arr,om,oR,ol)
    fracoR = norm_oR(a_arr,om,oR,ol)
    fracol = norm_ol(a_arr,om,oR,ol)
    plt.plot(t3_lookback_Gyr,fracom,'g', t3_lookback_Gyr, fracoR, 'r', t3_lookback_Gyr, fracol, 'b')
    plt.xlabel('Lookback time (Gyr)')
    plt.ylabel('Normalised density')
    
"""""

"""
# Calculate for the universe we think we live in, with approximately matter density 0.3 and cosmological constant 0.7
om = 0.9
ol = 0.1

# Note that when you integrate something with more than one argument you pass it with args=(arg1,arg2) in the integrate function
# e.g. "integrate.quad(adotinv, lower_limit, uper_limit, args=(om,ol))""
t3_lookback_Gyr = np.array([integrate.quad(adotinv, 1, a_end, args=(om,ol))[0] for a_end in a_arr])/H0y

om_arr = np.arange(0.1,2.1,0.4)

for om in om_arr:
    t4_lookback_Gyr = np.array([integrate.quad(adotinv, 1, a_end, args=(om,ol))[0] for a_end in a_arr])/H0y
    plt.plot(t4_lookback_Gyr,a_arr,label='$(\Omega_M,\Omega_\Lambda)$=(%.1f,%.1f)'%(om,ol))
 
plt.axvline(x=0,linestyle=':') # Plot some crosshairs 
plt.axhline(y=1,linestyle=':')
plt.xlabel('Lookback time (Gyr)')
plt.ylabel('Scalefactor')
plt.legend(loc='lower right',frameon=False)
plt.show()
"""











