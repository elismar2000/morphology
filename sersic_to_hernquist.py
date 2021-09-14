from astropy.io import fits
import numpy as np
from scipy import integrate

#=========================
#Calculating bulge galaxy physical parameters
#=========================

Ie = 8.90029 #DN s^-1
re = 1.85068 #pixel
n  = 2.02058
bulge_radius = 5.0 #estimated (by eye) total radius of the bulge

#d  = 38 #Mpc (Theureau et al. 2007)
# image = '/home/elismar/Documentos/Fisica/IC/imfit-1.7.1/ngc2992-93/images/NICMOS/n4sb08040/n4sb08040_mos.fits'
# hdr = fits.getheader(image)
# PHOTFLAM = hdr['PHOTFLAM']
#
# Ie *= PHOTFLAM * 4 * np.pi * d**2 #converting Ie to units of luminosity, i.e., erg s^-1 angstrom^-1
#
# mpc_to_cm = 3.086e+24
# d *= mpc_to_cm #converting distance d from Mpc to cm
#
# pixel_scale = 0.2
# re *= pixel_scale * d #this is re in cm
#
# bulge_radius = 5 * pixel_scale * d #approximate total radius os the bulge

#===========================
#Applying Abel's integral
#===========================
import math

b = 2*n - (1/3)
result = []
epsilon = 1e-5
for r in np.arange(0.1, bulge_radius, bulge_radius/30):
    def func(R):
        '''
        Function to be integrated by the Romberg method
        '''
        dI_dR = -(b * np.exp(-b * (-1 + (R/re)**(1/n))) * Ie * (R/re)**(-1 + (1/n))) / (n * re)
        return (-1/np.pi) * dI_dR * (1/np.sqrt(R**2 - r**2))

    i = integrate.romberg(func, r + epsilon, 1e+80, show=True)
    result.append(i)

result = np.asarray(result)
a = np.arange(0, len(result), 1)

import matplotlib.pyplot as plt
plt.plot(a, result)
plt.show()
