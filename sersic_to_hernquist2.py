import matplotlib.pyplot as plt
import numpy as np

#third part
import abel
from astropy.io import fits
from lmfit import Model, Parameters
from photutils.isophote import EllipseGeometry, Ellipse

#=======================
#Performing the inverse Abel transformation over the bulge 2D image
#=======================

#This must be a 2D model image of the bulge to be transformed
image = '/home/elismar/Documentos/Fisica/IC/imfit-1.7.1/ngc2992/images/2MASS/bulge_model.fits'
image = fits.getdata(image)

#The center pixel coordinates of the bulge
x0 = 98.5
y0 = 98.5
bulge_3D = abel.Transform(image, direction='inverse', method='three_point', center=(y0, x0)).transform

#=======================
#Fitting ellipses to the 3D image in order to get its brigthness profile
#=======================

# geometry = EllipseGeometry(x0=98, y0=98, sma=15, eps=0.19, pa=(np.pi/180)*226)
geometry = EllipseGeometry(x0=98, y0=98, sma=15, eps=0.0, pa=0)

ellipse = Ellipse(bulge_3D, geometry)
isolist = ellipse.fit_image(fix_pa=True, fix_eps=True)

#========================
#Converting to mass profile
#========================

#Considering 2MASS H-band of course
magzp = 20.4107
distance = 38e+6 #pc
solar_abs_mag = 3.32
gamma = 1.0

def count2mass(image, magzp, d, solar_abs_mag, gamma):
    app_mag = magzp - (2.5 * np.log10(image))
    abs_mag = app_mag - (2.5 * np.log10((d/10)**2))
    lum = 10**(-0.4*(abs_mag - solar_abs_mag))
    mass = gamma * lum
    return mass

mass = count2mass(isolist.intens, magzp, distance, solar_abs_mag, gamma)
mass = mass[1:]

#=========================
#Converting sma from px to pc
#=========================

ps = 1.0
ps *= (1/3600) * np.pi/180 #from arcsecs to radians

sma = isolist.sma[1:]
sma *= ps * distance

#=========================
#Fitting Hernquist profile
#=========================

def Hernquist(r, M, a):
    return (M*a) / (2 * np.pi * r * (r+a)**3)

model = Model(Hernquist)
params = Parameters()
params.add('M', value=1e+10)
params.add('a', value=1e+3)
result = model.fit(mass, params, r=sma)

print(result.fit_report())

plt.plot(sma, mass, 'X', color='yellow')
plt.plot(sma, Hernquist(sma, result.params['M'].value, result.params['a'].value), '--', color='blue')
plt.xlabel(r'$r\ [pc]$')
plt.ylabel(r'$\rho(r)\ [\frac{M_\odot}{pc^3}]$')

plt.show()
