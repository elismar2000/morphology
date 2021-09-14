import numpy as np
import matplotlib.pyplot as plt


def M_r_bulge(r, M_bulge, a):
    return M_bulge * r**2 / (r + a)**2

def M_r_disk(r, M_disk, Rd):
    return M_disk * ( (-r/Rd)*np.exp(-r/Rd) - np.exp(-r/Rd) + 1 )

M_bulge = 1.4e+10  #Msun
a = 0.38  #kpc
M_disk = 5.0e+10 #Msun
Rd = 2.1 #kpc

half_mass = (M_bulge + M_disk) / 2

radii = np.arange(0 + Rd/1000, 10*Rd, Rd/1000)

cumulative_stellar_mass = M_r_bulge(radii, M_bulge, a) + M_r_disk(radii, M_disk, Rd)

plt.plot(radii, cumulative_stellar_mass, label='bulge + disk stellar masses')
plt.axhline(half_mass, color='red', label='Half stellar mass')
plt.legend()
plt.show()
