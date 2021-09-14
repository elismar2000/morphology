import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits

#=============================
#Setting arguments
#=============================

bulge_model = '/home/elismar/Documentos/Fisica/IC/imfit-1.7.1/ngc2992/images/2MASS/bulge_model.fits'
#bulge_model = '/home/elismar/Documentos/Fisica/IC/imfit-1.7.1/ngc2992/images/2MASS/NGC2993_2MASS_h-band.fits'
magzp = 20.4107 #h-band
solar_abs_mag = 3.32 #h-band
d  = 38e+6 #pc
gama = 1.0
ps = 1.0 #is it really the pixel scale for 2MASS?
re = 3.73

#=============================
#Calculating bulge mass
#=============================

bulge_model = fits.getdata(bulge_model)

tot_flux = np.sum(bulge_model)

app_mag = magzp - 2.5 * np.log10(tot_flux)
abs_mag = app_mag - 2.5 * np.log10((d/10)**2)
lum = 10**(-0.4 * (abs_mag - solar_abs_mag))

mass = gama * lum

print('The total mass of the bulge in solar masses is {:.2e}'.format(mass))

#==============================
#Calculating parameter a
#==============================

ps *= (1/3600) * np.pi/180 #from arcsecs to radians

re *= ps * d

a = re / 1.8153

print('The value of parameter a in parsecs is {:.2f}'.format(a))
