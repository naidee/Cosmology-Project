#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 16:53:10 2024

@author: sakhi
"""

# First let's set up our packages
import numpy as np
from matplotlib import pyplot as plt
from scipy import integrate

# And set some constants
c = 299792.458  # km/s (speed of light)

H0kmsmpc = 70.  # Hubble constant in km/s/Mpc
cH0mpc = c / H0kmsmpc  # c/H0 in Mpc (the km/s cancel out in the numerator and denominator)
cH0Glyr = cH0mpc * 3.262 / 1000  # c/H0 in billions of light years. There are 3.262 light year / parsec

# Write a function for the integrand, i.e. $1/E(z)$,
def Ezinv(z, om, ol):  # atm assumes or is zero
    Ez = np.sqrt((om * (1 + z) ** 3) + ((1.0 - om - ol) * (1 + z) ** 2) + (ol))
    return 1.0 / Ez

# Start by making an array of redshifts
zstart = 0.0
zstop = 1000
zstep = 0.1 # Make this finer to make the plot smoother
zarr = np.arange(zstart,zstop,zstep)
# print('zarr=',zarr)

with_a_emit = np.zeros(len(zarr))
# Now add your code to calculate distance vs redshift and then plot it.  
xarr = np.zeros(len(zarr))



# You may find this function useful to calculate the perpendicular comoving distance R0*Sk(X)
    # Corrects the comoving distance, xx, for curvature, ok.
    # Result is perpendicular comoving distance / (c/H0)  
    # i.e. it's the distance to use in angular diameter and luminosity distance
def Sk(xx, ok):
    if ok < 0.0:
        dk = np.sin(np.sqrt(-ok)*xx)/np.sqrt(-ok)
    elif ok > 0.0:
        dk = np.sinh(np.sqrt(ok)*xx)/np.sqrt(ok)
    else:
        dk = xx
    return dk
#"This whole function spits out: R_0Sk(\chi) divded by c/H_0"


#om_range = (0,1,0.1)
#ol_range = (0,0.4,0.1)

ol=0.7
om_range = np.array([0,0.3,0.6])

for i, z in enumerate(zarr):
    for OM in om_range:
                xarr[i] = integrate.quad(Ezinv,0,z,args=(OM,ol))[0]
                a_emit = (1+z)**-1
                with_a_emit[i] = a_emit*xarr[i]

    
# Sub in the required constants to get the comoving distance R_0*X
R0X = xarr*cH0Glyr # Distance in Glyr
D_prop = with_a_emit*cH0Glyr



 
for om in om_range:
    R0X = np.zeros(len(zarr))
    ok = 1.0-om-ol
    #CD = R0X                          # COMOVING distance
    DL = Sk(xarr,ok)*(1+zarr)*cH0Glyr # Luminosity distance
    DA = Sk(xarr,ok)*cH0Glyr/(1+zarr)  
    for i, z in enumerate(zarr):
        x = integrate.quad(Ezinv, 0, z, args=(om, ol))[0]
        R0X[i] = x * cH0Glyr  # Distance in Glyr   
    # Angular diameter distance
    #plt.plot(zarr,CD,label='Comoving Distance')
    label_str = f'$\Omega_k$={ok:.2f}'
    plt.plot(zarr,DL,label=f'{label_str} LD',linestyle='--')
    plt.plot(zarr,DA,label=f'{label_str} ADD')
    plt.plot(zarr, R0X, label=f'{label_str} CD',linestyle=':')
        
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small')
plt.xlim(0.001, 1000)
plt.xscale('log') 
plt.ylim(0.001, 1000)
plt.yscale('log')
plt.xlabel('Redshift')
plt.ylabel('$Distances$ (Glyr)')
plt.show()
