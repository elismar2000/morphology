import matplotlib.pyplot as plt
import numpy as np

from lmfit import Model, Parameters

#=======================
#Setting parameters
#=======================

ps = 1.0
d  = 38 #Mpc
magzp = 20.4107 #h-band
solar_abs_mag = 3.32 #h-band
Ie = 49.24 #countrate
re = 3.73 #pixels
n = 1.41
gama = 1.0

#===========================
#Generating 1D Sersic profile
#===========================

#Translating re from pixels to parsecs
ps *= (1/3600) * np.pi/180 #from arcsecs to radians
d *= 1e+6 #from Mpc to pc

re *= ps * d

#Translating Ie from countrate to luminosity in solar units
Ie_app_Mh = magzp - 2.5 * np.log10(Ie)
Ie_abs_Mh = Ie_app_Mh - 2.5 * np.log10((d/10)**2)
Ie_lum = 10**(-0.4 * (Ie_abs_Mh - solar_abs_mag))

bn = n - (1/3)

R = np.arange(re/100, 3.1*re, re/100)
S = Ie_lum * np.exp(-bn * ((R/re)**(1/n) - 1))

#============================
#Fitting Hernquist-based model
#============================

def I(R, a, M):
    X = []
    for r in R:
        s = r/a
        if (s<=1) and (s>=0):
            X.append((1/np.sqrt(1 - s**2)) * np.arccosh(1/s))
        elif (s>=1):
            X.append((1/np.sqrt(s**2 - 1)) * np.arccos(1/s))
    X = np.array(X)
    return (M / (2*np.pi * a**2 * gama * (1 - s**2)**2)) * ((2 + s**2) * X - 3)

model = Model(I)
params = Parameters()
params.add('M', value=1e+10)
params.add('a', value=1e+3)
result = model.fit(S, params, R=R)

print(result.fit_report())

plt.plot(R, S, '--X', color='yellow', label='Sersic Profile')
plt.plot(R, I(R, result.params['a'].value, result.params['M'].value), '--', color='blue', label='Hernquist-based model')
plt.xlabel(r'$r\ [pc]$')
plt.ylabel(r'$\rho(r)\ [\frac{M_\odot}{pc^3}]$')

plt.legend()
plt.show()
