"""
Created on Sun April 16 15:42:39 2024

@author: MYMTEAM
"""


# For each graph if you want to download the graphs directly with names
# just uncomment the pdf or png lines (I didn't want a million graphs while debugging)

import numpy as np
import os 
from matplotlib import pyplot as plt
from scipy import integrate

### definitions for variables ###

def adotinv(a, om, ol, orad):
    ok = 1 - om - ol - orad
    adot = a * np.sqrt(orad * a**-4 + om * a**-3 + ok * a**-2 + ol)
    return 1.0/adot

def aadotinv(a, om, ol, orad):
    ok = 1 - om - ol - orad
    aadot = a**2 * np.sqrt(orad * a**-4 + om * a**-3 + ok * a**-2 + ol)
    return 1.0/aadot

def hubble(a, om, ol, orad):
    ok = 1 - om - ol - orad
    adotovera =  np.sqrt(orad * a**-4 + om * a**-3 + ok * a**-2 + ol)
    return 1 / adotovera

def Ez(z, om, ol, orad):
    ok = 1 - om - ol - orad
    Ez = np.sqrt(orad * (1 + z)**4 + om * (1 + z)**3 + ok * (1 + z)**2 + ol) 
    return 1 / Ez

def Sk(xx, om, ol, orad):
    ok = 1 - om - ol - orad
    if ok < 0:
        dk = np.sin(np.sqrt(-ok) * xx) / np.sqrt(-ok)
    elif ok > 0:
        dk = np.sinh(np.sqrt(ok) * xx) / np.sqrt(ok)
    else:
        dk = xx
    return dk

#plt.rcParams.update(plt.rcParamsDefault)   #only uncomment this if you want to restore the default font
plt.rcParams['font.family'] = 'Serif'


#def main():
dir_path = os.path.dirname(os.path.realpath(__file__))

### --- Constants --- ###
c = 299792.458 # km/s (speed of light)
H0kmsmpc = 70.  # Hubble constant in km/s/Mpc
H0s = H0kmsmpc * 3.2408e-20 # H0 in inverse seconds is H0 in km/s/Mpc * (3.2408e-20 Mpc/km)
H0y = H0s * 3.154e7 * 1.e9 # H0 in inverse Giga years is H0 in inverse seconds * (3.154e7 seconds/year) * (1e9 years / Giga year)
cH0mpc = c/H0kmsmpc   # c/H0 in Mpc  (the km/s cancel out in the numerator and denominator)
cH0Glyr = cH0mpc * 3.262 / 1000 #c/H0 in billions of light years.  There are 3.262 light year / parsec
    




### --- Scalefactor Graphs --- ### 

# Start by making an array of scalefactors
astart, astop, astep = 0, 2.1, 0.05
a_arr = np.arange(astart, astop, astep)

# Calculate for the universe we think we live in, with approximately matter density 0.3 and cosmological constant 0.7
om, ol, orad = 0.3, 0.7, 0.00005

ok = 1 - om - ol - orad
Rscale = a_arr * c / (np.sqrt(abs(ok)) * H0kmsmpc)
# Note that when you integrate something with more than one argument you pass it with args=(arg1,arg2) in the integrate function
# e.g. "integrate.quad(adotinv, lower_limit, uper_limit, args=(om,ol))""
t_lookback_Gyr = np.array([integrate.quad(adotinv, 1, upper_limit, args=(om,ol, orad))[0] for upper_limit in a_arr])/H0y

""" scalefactor vs Lookback time (Our Universe only) """
fig, ax = plt.subplots()
ax.plot(t_lookback_Gyr,a_arr,label='$(\Omega_M,\Omega_\Lambda)$=(%.2f,%.2f)'%(om,ol)) 
ax.axvline(x=0,linestyle=':'); plt.axhline(y=1,linestyle=':') # Plot some crosshairs
ax.plot(0,1,"ro") 
ax.set_xlabel('Lookback Time (Gyr)'); ax.set_ylabel('Scalefactor ($a$)')
ax.legend(loc='upper left',frameon=False)
ax2 = ax.twinx()        #makes a twin y axis so that we can compare a scale factor with R scale factor
ax2.plot(t_lookback_Gyr, Rscale, alpha=0); ax2.set_ylabel("Scalefactor $R$")   #plot an invisible curve so that there is the R scale factor axis too
ax2.ticklabel_format(axis='y', style='scientific', useMathText=True)
#fig.savefig(dir_path+'\\Part One Graphs\\Scalefactor-vs-LookbackTime for Our Universe.png', dpi=200, bbox_inches='tight', pad_inches = 0.01)
#fig.savefig(dir_path+'\\Part One Graphs\\Scalefactor-vs-LookbackTime for Our Universe.pdf', dpi=200, bbox_inches='tight', pad_inches = 0.01)
plt.show()
plt.close(fig)


"""" Scalefactor vs Lookback time (multiple matter densities)"""
om_arr = np.arange(0,2.1,0.4)
fig, ax = plt.subplots()
ax2 = ax.twinx()
for OM in om_arr:
    t_lookback_Gyr = np.array([integrate.quad(adotinv, 1, upper_limit, args=(OM, ol, orad))[0] for upper_limit in a_arr]) / H0y  #calculates lookback time for some om 
    ax.plot(t_lookback_Gyr, a_arr, label='$(\Omega_M,\Omega_\Lambda)$=(%.1f,%.1f)'%(OM,ol))     #plots the scale factor vs lookback time
    ax2.plot(t_lookback_Gyr, Rscale)    #as above, but for R scale factor too

ax.axvline(x=0,linestyle=':'); ax.axhline(y=1,linestyle=':') #Plot some crosshairs 
ax.set_xlabel('Lookback time (Gyr)'); ax.set_ylabel('Scalefactor ($a$)')  
ax.legend(loc='upper left',frameon=False)
ax2.set_ylabel("Scale Factor $R$"); ax2.ticklabel_format(axis='y', style='scientific', useMathText=True)
#fig.savefig(dir_path+'\\Part One Graphs\\Scalefactor-vs-LookbackTime.png', dpi=200, bbox_inches='tight', pad_inches = 0.01)
#fig.savefig(dir_path+'\\Part One Graphs\\Scalefactor-vs-LookbackTime.pdf', dpi=200, bbox_inches='tight', pad_inches = 0.01)
plt.show()
plt.close(fig)

""" scalefactor vs Lookback time (varying both densities) """
om_arr = np.array([0.3,1,5]) # played around with different values these seemed nicest/easiest to format
ol_0 = 0
fig, ax = plt.subplots()
ax2 = ax.twinx()
ax.plot(t_lookback_Gyr,a_arr,label='$(\Omega_M,\Omega_\Lambda)$=(%.2f,%.2f)'%(om,ol))
ax2.plot(t_lookback_Gyr, Rscale)    #as above, but for R scale factor too
for OM in om_arr:
    t_lookback_Gyr = np.array([integrate.quad(adotinv, 1, upper_limit, args=(OM, ol_0, orad))[0] for upper_limit in a_arr]) / H0y  #calculates lookback time for some om 
    ax.plot(t_lookback_Gyr, a_arr, label='$(\Omega_M,\Omega_\Lambda)$=(%.1f,%.1f)'%(OM,ol_0))     #plots the scale factor vs lookback time
    ax2.plot(t_lookback_Gyr, Rscale)    #as above, but for R scale factor too

ax.axvline(x=0,linestyle=':'); ax.axhline(y=1,linestyle=':') #Plot some crosshairs 
ax.set_xlabel('Lookback time (Gyr)'); ax.set_ylabel('Scalefactor ($a$)')  
ax.legend(loc='upper left',frameon=False)
ax2.set_ylabel("Scale Factor $R$"); ax2.ticklabel_format(axis='y', style='scientific', useMathText=True)
#fig.savefig(dir_path+'\\Part One Graphs\\Scalefactor-vs-LookbackTime2.png', dpi=200, bbox_inches='tight', pad_inches = 0.01)
#fig.savefig(dir_path+'\\Part One Graphs\\Scalefactor-vs-LookbackTime2.pdf', dpi=200, bbox_inches='tight', pad_inches = 0.01)
plt.show()
plt.close(fig)

"""" Redshift Vs Scalefactor """
fig, ax = plt.subplots()
redshift = -1 + (1 / a_arr[1:])         #this is the formula for z in terms of a
ax.plot(a_arr[1:], redshift)
ax.set_xlabel("Scalefactor ($a$)"); ax.set_ylabel("Redshift ($z$)")
#fig.savefig(dir_path+'\\Part One Graphs\\Redshift vs Scalefactor.png', dpi=200, bbox_inches='tight', pad_inches = 0.01)
#fig.savefig(dir_path+'\\Part One Graphs\\Redshift vs Scalefactor.pdf', dpi=200, bbox_inches='tight', pad_inches = 0.01)
plt.show()
plt.close(fig)


""" Lookback time vs Redshift """
amin, amax, astep = 0.1, 1, 0.05
a_arr = np.arange(amin, amax + astep, astep)
redshift = -1 + (1 / a_arr)         #calculates redshift from this new array
t_lookback_Gyr = np.array([integrate.quad(adotinv, 1, upper_limit, args=(om, ol, orad))[0] for upper_limit in a_arr]) / H0y
fig, ax = plt.subplots()
ax.plot(redshift, t_lookback_Gyr)
ax.set_xlabel("Redshift ($z$)"); ax.set_ylabel("Lookback Time (Gyr)")
ax.set_xlim(xmin=0); ax.set_ylim(ymax=0)
#fig.savefig(dir_path+'\\Part One Graphs\\Lookback Time vs Redshift.png', dpi=200, bbox_inches='tight', pad_inches = 0.01)
#fig.savefig(dir_path+'\\Part One Graphs\\Lookback Time vs Redshift.pdf', dpi=200, bbox_inches='tight', pad_inches = 0.01)
plt.show()
plt.close(fig)



""" Hubble Parameter vs Time """
amin, amax, astep, scalestep = 0.02, 2, 0.02, 0.2
a_arr = np.arange(amin, amax + astep, astep)
times = np.zeros(len(a_arr))
fig, ax = plt.subplots()
ax2 = ax.twiny()        #create alternate x-axis to sit at the top (represents scalefactor a)
scaletickslabs, scaletickslocs = [0], [0]       #initialize lists for the alternate (top) x-axis ticks
for i, a in enumerate(a_arr):
    times[i] = integrate.quad(adotinv, 0, a, args=(om, ol, orad))[0] / H0y      #calculates age of universe at scale factor a
    if a == 1:
        ax.axvline(x=times[i],linestyle=':')        #plots a vertical line at current scale factor (age of universe)
for a in np.arange(0, amax + scalestep, scalestep):
    if a == 0:
        pass        #we want to avoid a=0 or else encounter a math error
    else:
        time = integrate.quad(adotinv, 0, a, args=(om, ol, orad))[0] / H0y      #calculate at which time to place the scalefactor tick a
        scaletickslocs.append(time)     #location of this scalefactor tick with respect to the time axis
        scaletickslabs.append(round(a, 1))  #label of this scalefactor tick
hubble_arr = H0kmsmpc * (adotinv(a_arr, om, ol, orad) * a_arr)**-1    #inverse due to the adotinv function being the reciprocal
ax.plot(times, hubble_arr)
ax.set_xlabel("Age of Universe (Gyr)"); ax.set_ylabel("Hubble Parameter (km/s/Mpc)")
ax.set_yscale('log')
ax.axhline(y=H0kmsmpc,linestyle=':')
ax.set_xlim(xmin=min(a_arr))

mn, mx = ax.get_xlim()
ax2.set_xlim(mn, mx)    #set top x-axis to have same width as bottom x-axis
ax2.set_xticks(scaletickslocs); ax2.set_xticklabels(scaletickslabs)
ax2.set_xlabel("Scalefactor ($a$)")
#fig.savefig(dir_path+'\\Part One Graphs\\Hubble Parameter vs Time.png', dpi=200, bbox_inches='tight', pad_inches = 0.01)
#fig.savefig(dir_path+'\\Part One Graphs\\Hubble Parameter vs Time.pdf', dpi=200, bbox_inches='tight', pad_inches = 0.01)
plt.show()
plt.close(fig)




### --- Redshift Graphs --- ###

# Start by making an array of redshifts
zstart, zstop, zstep = 0, 4, 0.01
zarr = np.arange(zstart, zstop + zstep, zstep)

# Now add your code to calculate distance vs redshift and then plot it.  
xarr = np.zeros(len(zarr))
axarr = np.zeros(len(zarr))
for i, z in enumerate(zarr):
    xarr[i] = integrate.quad(Ez, 0, z, args=(om, ol, orad))[0] 
    axarr[i] = (1 / (1 + z)) * integrate.quad(Ez, 0, z, args=(om, ol, orad))[0]     #this is effectively a * R0x - distance to light source seen at redshift z at time of light emission
    
# Sub in the required constants to get the comoving distance R_0*X
R0X = xarr * cH0Glyr # Distance in Glyr
aR0X = axarr * cH0Glyr # Distance in Glyr

""" Comoving distance vs redshift """
fig, ax = plt.subplots()
ax.plot(zarr,R0X)
ax.set_xlabel('Redshift ($z$)'); ax.set_ylabel('$R_0\chi$ (Glyr)')
ax.set_xlim(xmin=0); ax.set_ylim(ymin=0)
#fig.savefig(dir_path+'\\Part One Graphs\\Comoving Distance vs Redshift.png', dpi=200, bbox_inches='tight', pad_inches = 0.01)
#fig.savefig(dir_path+'\\Part One Graphs\\Comoving Distance vs Redshift.pdf', dpi=200, bbox_inches='tight', pad_inches = 0.01)
plt.show()
plt.close(fig)