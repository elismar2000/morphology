#!/usr/bin/env python
# coding: utf-8

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

def stellar_mass_halo(m_halo):

    E_o = 0.0351
    Gamma = 0.608
    Beta = -1.376
    Mo = 11.590

    Mr = 10**m_halo/10**Mo

    E_n = 2*E_o*(1./ (Mr**Beta + Mr**Gamma))

    return E_n

M_halo_range = np.linspace(10,14.5,100)
E_range = stellar_mass_halo(M_halo_range)

M_gal = np.log10(10**(M_halo_range)*E_range)
mhalo_func = interp1d(M_gal, M_halo_range,
                      bounds_error=False,
                      fill_value=(np.nan, np.nan))


fig, axs = plt.subplots(1, 2)
axs[0].plot(M_halo_range, np.log10(E_range), label='Moster et al. 2013')
axs[0].set_xlabel('log(M_halo)')
axs[0].set_ylabel('log (M_star/M_halo)')

axs[1].plot(M_gal, M_halo_range, label='Moster et al. 2013')
axs[1].set_ylabel('log(M_halo)')
axs[1].set_xlabel('log (M_star)')
axs[1].axvline(np.log10(2.88e+10), label='NGC2992', color='purple')
axs[1].axvline(np.log10(2.57e+10), label='NGC2993', color='orange')
axs[1].axhline(np.log10(55.0e+10), color='purple')
axs[1].axhline(np.log10(20.5e+10), color='orange')

plt.legend()
plt.show()
