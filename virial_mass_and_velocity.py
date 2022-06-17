#Esse código calcula alguns dos parâmetros necessários pra gerar as galáxias no Dice
#A princípio, os únicos inputs que ele precisa são o redshift e os parâmetros dos componentes da galáxia
#Aqui estou assumindo halo/bojo --> Hernquist, discos --> perfil exponencial

from astropy.cosmology import FlatLambdaCDM
from astropy.constants import G
import astropy.units as u

import matplotlib.pyplot as plt
import numpy as np


#===========================
#Definindo a densidade crítica do Universo no tempo de formação das galáxias
#===========================
#Redshift das galáxias
z = 0.00771

#Definindo a cosmologia e o parâmetro de Hubble no tempo de formação das galáxias
#os wmap carregam os valores dos parametros, esse aqui viram do ned
cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)

H = cosmo.H0 * cosmo.efunc(z=z)

#Densidade critica
rho_c = 3 * H**2 / (8 * np.pi * G)
rho_c = rho_c.to(u.Msun / u.kpc**3)


#===========================
#Definindo o perfil cumulativo de massa das galaxias
#===========================

def cumulative_hernquist(r, M, a):
    return M * r**2 / (r + a)**2

def cumulative_exponential(r, M, Rd, z0):
    return M * ( (-r/Rd)*np.exp(-r/Rd) - np.exp(-r/Rd) + 1 ) * np.tanh(r/z0)


dr = 0.5
radii = np.arange(120, 180+dr, dr)

#Parâmetros dos perfis das galáxias
#Massas em 10^10 M_sun, comprimentos em Kpc
#Esses valores aí são os que eu tava usando antes do ajuste que o Jose fez pro halo
galaxy = 'NGC 2993'
print(galaxy)

if galaxy == 'NGC 2992':
    M_halo = 55.0
    a_halo = 15.0

    M_disk = 2.3
    Rd_disk = 2.1
    z0_disk = 0.42

    M_bulge = 0.58
    a_bulge = 0.38

    M_gas = 0.16
    Rd_gas = 4.2
    z0_gas = 0.2 * z0_disk


if galaxy == 'NGC 2993':
    M_halo = 20.5
    a_halo = 5.58

    M_disk = 1.46
    Rd_disk = 0.802
    z0_disk = 0.266

    M_bulge = 0.11
    a_bulge = 0.21

    M_gas = 0.46
    Rd_gas = 1.604
    z0_gas = 0.2 * z0_disk


mass_halo_r = cumulative_hernquist(radii, M_halo, a_halo)
mass_disk_r = cumulative_exponential(radii, M_disk, Rd_disk, z0_disk)
mass_bulge_r = cumulative_hernquist(radii, M_bulge, a_bulge)
mass_gas_r = cumulative_exponential(radii, M_gas, Rd_gas, z0_gas)

masses_from_models = mass_halo_r + mass_disk_r + mass_bulge_r + mass_gas_r

volume = (4/3) * np.pi * radii**3

density = masses_from_models / volume


#===========================
#Plotando os resultados, com uma linha horizontal pra demarcar 200*rho_c
#===========================
#Perceber que as massas no Galstep são em 1e+10M_sun, mas aqui estamos deixando tudo em M_sun

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(radii, density)
ax.set_ylabel(r'$\rho(r)\ [\frac{10^{10} M_{\odot}}{Kpc^3}]$')
ax.set_xlabel(r'$r\ [Kpc^3]$')
ax.axhline(200 * rho_c.value / 1e+10, color='firebrick')
plt.show()


#===========================
#Determinando o r_200 e calculando M200 e v200
#===========================

r200_index = np.abs(density - (200 * rho_c.value / 1e+10)).argmin()
r200 = radii[r200_index] * u.kpc

M200 = masses_from_models[r200_index] * 1e+10 * u.Msun
v200 = np.sqrt(G * M200 / r200).to(u.km / u.s)

print('r_200 = {:.3f}'.format(r200))
print('M_200 = {:e}'.format(M200))
print('v_200 = {:.3f}'.format(v200))


#===========================
#Determinando frações de massa dos componentes
#===========================

M200_halo = mass_halo_r[r200_index] * 1e+10
M200_gas = mass_gas_r[r200_index] * 1e+10
M200_disk = mass_disk_r[r200_index] * 1e+10
M200_bulge = mass_bulge_r[r200_index] * 1e+10

print('mass_frac_halo = {:.3f}'.format(M200_halo/M200.value))
print('mass_frac_gas = {:.3f}'.format(M200_gas/M200.value))
print('mass_frac_disk = {:.3f}'.format(M200_disk/M200.value))
print('mass_frac_bulge = {:.3f}'.format(M200_bulge/M200.value))
