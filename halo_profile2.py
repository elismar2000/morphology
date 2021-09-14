import astropy.units as u
from lmfit import Model, Parameters
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cons

#========================================

light_curve = np.genfromtxt('ngc2992_vel_curve.dat')

radii_neg = np.sort(light_curve[:, 8][light_curve[:, 8] < 0])[::-1] * -1
radii = radii_neg

v_obs = np.array([light_curve[:, 9][light_curve[:, 8] == -r][0] for r in radii_neg])

v_obs -= v_obs[0]
v_obs = abs(v_obs) * u.km / u.s

#deleting the first few due to resolution issues
radii = radii[3:-3]
v_obs = v_obs[3:-3]

#========================================

def M_r_bulge(r, M_bulge, a):
    return M_bulge * r**2 / (r + a)**2

def M_r_disk(r, M_disk, Rd):
    return M_disk * ( (-r/Rd)*np.exp(-r/Rd) - np.exp(-r/Rd) + 1 )

M_bulge = 1.4e+10  #Msun
a = 0.38  #kpc
M_disk = 5.0e+10 #Msun
Rd = 2.1 #kpc

distance = 39e+3 #Kpc

radii = np.deg2rad(radii / 3600) * distance #Kpc

mass_lum = np.array([M_r_bulge(r, M_bulge, a) + M_r_disk(r, M_disk, Rd) for r in radii])

radii = radii * u.kpc
G = cons.G * u.m**3 / (u.kg * u.s**2)
mass_lum *= u.Msun
v_lum = np.sqrt(G * mass_lum / radii).to(u.km / u.s)

#========================================

M_dark = (radii / G) * (v_obs**2 - v_lum**2)
M_dark = M_dark.to(u.Msun)

fig, axs = plt.subplots(1, 2)
axs[0].plot(radii, v_obs, marker='s', color='orangered', label='Observed rotation curve')
axs[0].plot(radii, v_lum, 'o', color='gold', label='Keplerian rotation curve from luminous matter')
axs[0].set_ylabel('km/s', fontsize=20)
axs[0].set_xlabel('kpc', fontsize=20)

axs[1].plot(radii, M_dark, marker='x', color='slateblue', label='Calculated value of halo cumulative density profile')
axs[1].set_ylabel(r'$M_{\odot}$', fontsize=20)
axs[1].set_xlabel('kpc', fontsize=20)

#========================================

def Hernquist(r, M, a):
    return M * r**2/(r + a)**2

model = Model(Hernquist)
params = Parameters()
params.add('M', value=1e+12, min=1e+9, max=1e+13)
params.add('a', value=100, min=1, max=50)
result = model.fit(M_dark.value, params, r=radii.value)

print(result.fit_report())

axs[1].plot(radii, Hernquist(radii.value, result.params['M'].value, result.params['a'].value), linestyle='--',
    color='indigo', label='Line of best fit')

axs[0].legend()
axs[1].legend()
plt.show()
