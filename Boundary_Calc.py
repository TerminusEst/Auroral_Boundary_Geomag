"""
Script to calculate auroral equatorward boundary for a quiet and an active day
Appears in Blake et al. (2020) "Estimating Location of Auroral Equatorward Boundary 
using Historical and Simulated Surface Magnetic Field Data", JGR SP, DOI: XXXXXXXXXXXXX

This script calculates the boundaries in Fig. 3, and assumes magnetic latitudes
and maximum Eh have been previously calculated.

Uses the Python gcvspline library:  

"""

from gcvspline import SmoothedNSpline
import numpy as np
from matplotlib import pyplot as plt
############################################################################

def fit_spline_lat_e(lats, maxe, bootstrap_num = 500, lat_cutoff = 10, uselog = True):
    """Calculate Spline Bootstrap fit for input lats and max Eh. Calculates using
    random subselection bootstrap_num times, then calculates using all data
    
    inputs:
    -------------------------------------------------------------
    lats = magnetic latitude of sites in deg
    maxe = maximum calculated Eh at sites (usually V/km)
    bootstrap_num = number of times to calculate boundary with random 
                    subselection of points
    lat_cutoff = latitude South of which points are ignored for fitting
    uselog = whether to use log10 of maxe or not for fitting 
    
    outputs:
    -------------------------------------------------------------
    xthresh = all of the latitude threshold points calculated from bootstrap
    ythresh = all of the E values for each threshold point
    gradients = gradients at each point
    absfit = np.array([xx1, abs_yy1]) = the curve calculated with 100% of the points
    """

    X = lats
    Y = maxe
    ind = X >= lat_cutoff
    
    if uselog == True:
        x1, y1 = X[ind], np.log10(Y[ind])
    else:
        x1, y1 = X[ind], Y[ind]
        
    # This is the smoothing factor
    p = 400
    
    # Latitude points to calculate fit 
    xx1 = np.linspace(x1[0], x1[-1], 1000)

    # calculate the fit for 75% of the points, bootstrap_num times
    xthresh, ythresh, gradients = [], [], []
    while len(xthresh) < bootstrap_num:
        # Randomly select 75% of the points
        indices = np.random.choice(len(x1), int(len(x1)*(0.75)), replace=False)
        x1_new, y1_new = x1[indices], y1[indices]
        x1_new, y1_new = zip(*sorted(zip(x1_new, y1_new)))
        x1_new, y1_new = np.array(x1_new) + (np.arange(len(x1_new)) / 1e7), np.array(y1_new)

        # Calculate spline fit with this reduced set
        xx1 = np.linspace(x1[0], x1[-1], 1000)
        w1 = np.ones_like(y1_new)
        GCV_manual = SmoothedNSpline(x1_new, y1_new, w=w1, p=p)
        yy1 = GCV_manual(xx1)
        gradyy1 = np.gradient(yy1)
        xgradyy1 = xx1[gradyy1.argmax()]

        # If the fit misattributes it to a point below the lat_cutoff, ignore
        if xgradyy1 >= lat_cutoff:
            gradients.append(np.max(gradyy1))
            xthresh.append(xgradyy1)
            ythresh.append(yy1[gradyy1.argmax()])


    # Now calculate with all of the data
    w1 = np.ones_like(y1)
    GCV_manual2 = SmoothedNSpline(x1, y1, w=w1, p=p)
    abs_yy1 = GCV_manual2(xx1)
    abs_gradyy1 = np.gradient(abs_yy1)
    abs_xgradyy1 = xx1[abs_gradyy1.argmax()]
     
    return xthresh, ythresh, gradients, np.array([xx1, abs_yy1])

############################################################################

# Read in data, 1 = 2003-10-30, 2 = 2009-10-07
mlat1, maxe1 = np.loadtxt('Data/2003_10_30.txt', usecols = (0,1), unpack = True, skiprows = 1)
mlat2, maxe2 = np.loadtxt('Data/2009_10_07.txt', usecols = (0,1), unpack = True, skiprows = 1)

# Calculate the curves
xthresh1, ythresh1, gradients1, absfit1 = fit_spline_lat_e(mlat1, maxe1/1000., bootstrap_num = 500, lat_cutoff = 10, uselog = True)
xthresh2, ythresh2, gradients2, absfit2 = fit_spline_lat_e(mlat2, maxe2/1000., bootstrap_num = 500, lat_cutoff = 10, uselog = True)

# Get the boundaries (points where the gradient of the fit is greatest
absgrad1 = np.gradient(absfit1[1])
index1 = np.where(absgrad1 == np.max(absgrad1))
boundary1 = absfit1[0][index1]

# And for day 2:
absgrad2 = np.gradient(absfit2[1])
index2 = np.where(absgrad2 == np.max(absgrad2))
boundary2 = absfit2[0][index2]

############################################################################
# Now make the plot:

plt.clf()

ax1 = plt.subplot(111)

plt.title("Calculation of Auroral Boundary", fontsize = 18)

plt.plot(mlat1, np.log10(maxe1/1000.), 'ob', label = '30 Oct 2003')
plt.plot(absfit1[0], absfit1[1], 'b', lw = 3)#, label = 'fit')
plt.axvline(boundary1, lw = 2, color = 'blue', linestyle = 'dashed')

plt.plot(mlat2, np.log10(maxe2/1000.), 'or', label = '07 Oct 2009')
plt.plot(absfit2[0], absfit2[1], 'r', lw = 3)
plt.axvline(boundary2, lw = 2, color = 'red', linestyle = 'dashed')

plt.grid(True)
plt.ylabel("$log_{10}E$ (V/km)", fontsize = 18)

plt.xlim([0, 90])
plt.legend(loc = 'upper left', fontsize = 14)

plt.xlabel("Magnetic Latitude ($^{\circ}$N)", fontsize = 18)

plt.show()


















