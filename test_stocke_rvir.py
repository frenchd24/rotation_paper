from pylab import *

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


def stocke_rvir_interpolate(filename):
    """
    Interpolates the Stocke et al. (2013) Figure 1 to get Rvir as a function of L/L*
    
    Parameters
    ----------
    filename    : the directory to the csv file containing an extraction of the plot data
                  Two columns - x and y, both floats
    
    Returns
    -------
    interp      : interpolation object such that interp(L) returns the corresponding Rvir
    
    
    """
    # --- read in the data  and unpack
    log_L, log_rvir = np.loadtxt(filename, delimiter=',', usecols=(0,1), unpack=True)
    
    print('log_L[0]: ',log_L[0])
    print('log_rvir[0]: ',log_rvir[0])
    print()
    
    # --- get out of log space
    L = 10**log_L
    rvir = 10**log_rvir
    
    # --- interpolate it using a cubic spline
    interp = interp1d(L, rvir, kind='cubic')
    
    return interp
    
    
# --- interpolate the Stocke et al. (2013) Lstar vs Rvir relation
stocke_rvir_filename = '/Users/frenchd/Research/rotation_paper_data/Rvir_L_Lstar2.csv'
stocke_rvir = stocke_rvir_interpolate(stocke_rvir_filename)


fig = plt.figure(figsize=(7.7,5.7))
ax = fig.add_subplot(1,1,1)

x_fit = np.arange(0.00011, 80, 0.0001)
plt.plot(np.log10(x_fit), np.log10(stocke_rvir(x_fit)), 'r-', color='green', alpha=0.9, lw=2)


def lessOne(L):
    logR = 0.15 * np.log10(L) + 2.10
    return 10**logR

def middle(L):
    logR = 0.30*np.log10(L) + 2.25
    return 10**logR
    
def zeroPlus(L):
    logR = 0.33 * np.log10(L) + 2.25
    return 10**logR

x_fit = np.arange(0.00011, 0.1, 0.0001)
plt.plot(np.log10(x_fit), np.log10(lessOne(x_fit)), 'r-', color='red', alpha=0.7, lw=1.5)

x_fit = np.arange(0.1, 1, 0.0001)
plt.plot(np.log10(x_fit), np.log10(middle(x_fit)), 'r-', color='black', alpha=0.7, lw=1.5)

x_fit = np.arange(1.0, 80, 0.0001)
plt.plot(np.log10(x_fit), np.log10(zeroPlus(x_fit)), 'r-', color='blue', alpha=0.7, lw=1.5)

plt.ylabel(r'$\rm log(R_{vir}/kpc)$')
plt.xlabel(r'$\rm log(L/L^*)$')

# x-axis
majorLocator   = MultipleLocator(1.0)
majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
minorLocator   = MultipleLocator(0.1)
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)

# y axis
majorLocator   = MultipleLocator(0.5)
majorFormatter = FormatStrFormatter(r'$\rm %.1f$')
minorLocator   = MultipleLocator(0.1)
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)

plt.ylim(1.5, 3.5)
ax.grid(True)
fig.show()

# saveDirectory = '/Users/frenchd/Research/rotation_paper_data/'
# fig.savefig("{0}/stocke_2013_lstar.jpg".format(saveDirectory),dpi=300,bbox_inches='tight')



# print("stocke_rvir(1.1) = ",stocke_rvir(1.1))
# print("stocke_rvir(0.007) = ",stocke_rvir(0.007))
# print("stocke_rvir(0.16) = ",stocke_rvir(0.16))
# print("stocke_rvir(0.53) = ",stocke_rvir(0.53))
# print("stocke_rvir(1.0) = ",stocke_rvir(1.0))
# print("stocke_rvir(1.1) = ",stocke_rvir(1.1))
# print("stocke_rvir(1.2) = ",stocke_rvir(1.2))
# print("stocke_rvir(1.3) = ",stocke_rvir(1.3))
# print("stocke_rvir(1.4) = ",stocke_rvir(1.4))
# print("stocke_rvir(1.5) = ",stocke_rvir(1.5))
# print("stocke_rvir(0.21) = ",stocke_rvir(0.21))
print("stocke_rvir(0.0216986) = ",stocke_rvir(0.0216986))
print("stocke_rvir(2.09138) = ",stocke_rvir(2.09138))
print("stocke_rvir(0.21) = ",stocke_rvir(0.21))


