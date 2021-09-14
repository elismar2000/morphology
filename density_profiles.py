import numpy as np
import matplotlib.pyplot as plt

#z = 0
def exp(r, M, Rd, z0):
    return M/(4* np.pi * Rd**2 * z0) * np.exp(-r/Rd)

def hernquist(r, M, a):
    return (M / (2 * np.pi)) * (a/r) * 1/(r+a)**3

M_star = 5e+10
M_gas = 0.5e+10
M_bulge = 1.4e+10
M_halo = 0.55e+12

Rd_star = 2.0
Rd_gas = 2.0

z0_star = 0.5
z0_gas = 0.5 * 0.2

a_bulge = 0.38
a_halo = 14.92

r = np.linspace(a_halo/200, a_halo, 1000)


rho_star = exp(r, M_star, Rd_star, z0_star)
rho_gas = exp(r, M_gas, Rd_gas, z0_gas)
rho_bulge = hernquist(r, M_bulge, a_bulge)
rho_halo = hernquist(r, M_halo, a_halo)

#plt.plot(r, rho_star, label='Stellar disk')
#plt.plot(r, rho_gas, label='Gas disk')
plt.plot(r, rho_bulge, label='Bulge')
plt.plot(r, rho_halo, label='Halo')
plt.xlabel('r [Kpc]')
plt.ylabel(r'$\rho [M_{\odot}]$')
plt.legend()
plt.show()
