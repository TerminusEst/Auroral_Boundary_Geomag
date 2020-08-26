"""
Python to calculate E-field data using 1-minute INTERMAGNET data and a resistive
1-D Quebec ground profile. 

The INTERMAGNET data used here is taken from 3 days at the Eskdalemuir observatory.

1-D Z-tensor taken from http://www.digitalearthlab.com/tutorial/tutorial-1d-mt-forward/  

"""

import numpy as np
import math
import datetime
from matplotlib import pyplot as plt

############################################################################

def Z_Tensor_1D(resistivities, thicknesses, frequencies):
    """Calculate 1D Z-Tensor for given ground resistivity profile.

    Parameters
    -----------
    resistivities = array or list of resistivity values in Ohm.m

    thicknesses = array or list of thicknesses in m.
        **len(resistivities) must be len(thicknesses) + 1**

    frequencies = array or list of frequencies to get response of
    
    Returns
    -----------
    Z = complex array of Z tensor values
    
    Taken from:
    http://www.digitalearthlab.com/tutorial/tutorial-1d-mt-forward/  
    """
    if len(resistivities) != len(thicknesses) + 1:
        print("Length of inputs incorrect!")
        return 
    
    mu = 4*np.pi*1E-7; #Magnetic Permeability (H/m)
    n = len(resistivities);
    master_Z, master_absZ, master_phase = [], [], []

    for frequency in frequencies:   
        w =  2*np.pi*frequency;       
        impedances = list(range(n));
        #compute basement impedance
        impedances[n-1] = np.sqrt(w*mu*resistivities[n-1]*1j);
       
        for j in range(n-2,-1,-1):
            resistivity = resistivities[j];
            thickness = thicknesses[j];
      
            # 3. Compute apparent resistivity from top layer impedance
            #Step 2. Iterate from bottom layer to top(not the basement) 
            # Step 2.1 Calculate the intrinsic impedance of current layer
            dj = np.sqrt((w * mu * (1.0/resistivity))*1j);
            wj = dj * resistivity;
            # Step 2.2 Calculate Exponential factor from intrinsic impedance
            ej = np.exp(-2*thickness*dj);                     
        
            # Step 2.3 Calculate reflection coeficient using current layer
            #          intrinsic impedance and the below layer impedance
            belowImpedance = impedances[j + 1];
            rj = (wj - belowImpedance)/(wj + belowImpedance);
            re = rj*ej; 
            Zj = wj * ((1 - re)/(1 + re));
            impedances[j] = Zj;    
    
        # Step 3. Compute apparent resistivity from top layer impedance
        Z = impedances[0];
        phase = math.atan2(Z.imag, Z.real)
        master_Z.append(Z)
        master_absZ.append(abs(Z))
        master_phase.append(phase)
        #master_res.append((absZ * absZ)/(mu * w))
    return np.array(master_Z)
    
############################################################################
    
def E_Field_1D(bx, by, resistivities, thicknesses, timestep = 60., Z = None, calc_Z = True, pad = True):
    """Calculate horizontal E-field components given Bx, By, resistivities and thicknesses.
    
    Parameters
    -----------
    bx, by = array of Bx, By timeseries in nT

    resistivities = array or list of resistivity values in Ohm.m

    thicknesses = array or list of thicknesses in m.
        **len(resistivities) must be len(thicknesses) + 1**

    timestep = time between samples (default is 60. for minute sampling)
    
    Z = complex Z-tensor array. If not supplied, Z will be calculated from input
        resistivities and thicknesses
    
    Returns
    -----------
    ext, eyt = arrays of electric field components in mV/km
    """
    
    if pad == False:
        new_bx = bx
        new_by = by
    else:
        new_bx = np.concatenate((bx[:150], bx, bx[-150:][::-1]))
        new_by = np.concatenate((by[:150], by, by[-150:][::-1]))
    
    mu0 = 4*np.pi * 1e-7
    freq = np.fft.fftfreq(new_bx.size, d = timestep)
    freq[0] = 1e-100

    if calc_Z == True:  # if you need to calculate Z
        Z = Z_Tensor_1D(resistivities, thicknesses, freq)
        
    bx_fft = np.fft.fft(new_bx)
    by_fft = np.fft.fft(new_by)

    exw = Z * by_fft/mu0; 
    eyw = -1 * Z * bx_fft/mu0

    ext = 1e-3 * np.fft.ifft(exw).real
    eyt = 1e-3 * np.fft.ifft(eyw).real

    if pad == False:
        return ext, eyt
    else:
        return ext[150:-150], eyt[150:-150]
        
############################################################################    
# Read in Eskdalemuir 1 minute bx, by data (bastardised from INTERMAGNET data)    
bx, by = np.loadtxt("Data/ESK_B.txt", usecols = (2,3), unpack = True, skiprows = 1)
datestr, timestr = np.loadtxt("Data/ESK_B.txt", usecols = (0, 1), unpack = True, skiprows = 1, dtype = str)

# list of datetime objects
td = []
for d, t in zip(datestr, timestr):
    dt = datetime.datetime.strptime(d + " " + t, '%Y-%m-%d %H:%M:%S')
    td.append(dt)
    
# 1D Quebec Resistivity Profile:
Qres = np.array([20000., 200, 1000, 100, 3])
Qthick = 1000. * np.array([15, 10, 125, 200])
  
# calculate the 1-minute E-field components    
ex, ey = E_Field_1D(bx, by, Qres, Qthick, timestep = 60., Z = None, calc_Z = True, pad = True)


# Make Plot
plt.clf()

ax1 = plt.subplot(211)
plt.title("ESK, 57.8$^{\circ}$N", fontsize = 18)
plt.plot(td, bx, label = "$B_{X}$")
plt.plot(td, by, label = "$B_{Y}$")
plt.grid(True)

plt.ylabel("$B$ (nT)", fontsize = 18)
plt.legend(fontsize = 18)
plt.xlim([td[0], td[-1]])

#-----------------------------------------------------------------
ax2 = plt.subplot(212)

plt.plot(td, ey/1000., label = "$E_{Y}$")
plt.plot(td, ex/1000., label = "$E_{X}$")

plt.grid(True)

plt.ylabel("$E$ (V/km)", fontsize = 18)
plt.legend(fontsize = 18)
plt.xlim([td[0], td[-1]])

plt.xlabel("29-31 October 2003", fontsize = 18)
plt.show()












































































        
        
        
