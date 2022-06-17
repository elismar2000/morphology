import numpy as np
import matplotlib.pyplot as plt

from pygadgetreader import *



#Analytical models for cumulative Hernquist and disk profiles
def cumulative_hernquist(r, M, a):
    return M * r**2 / (r + a)**2

def cumulative_exponential(r, M, Rd, z0):
    return M * ( (-r/Rd)*np.exp(-r/Rd) - np.exp(-r/Rd) + 1 ) * np.tanh(r/z0)


#Reading positions and masses of NGC2992 initial condition
ngc2992 = '/home/elismar/Documentos/Fisica/IC/simulations_ICs/galstep/galstep/ngc2992.ic'

pos_halo = readsnap(ngc2992, 'pos', 'dm')
pos_disk = readsnap(ngc2992, 'pos', 'disk')
pos_gas = readsnap(ngc2992, 'pos', 'gas')
pos_bulge = readsnap(ngc2992, 'pos', 'bulge')
pos = np.concatenate((pos_halo, pos_disk, pos_gas, pos_bulge))

mass_halo = readsnap(ngc2992, 'mass', 'dm')
mass_disk = readsnap(ngc2992, 'mass', 'disk')
mass_gas = readsnap(ngc2992, 'mass', 'gas')
mass_bulge = readsnap(ngc2992, 'mass', 'bulge')
mass = np.concatenate((mass_halo, mass_disk, mass_gas, mass_bulge))


#Calculating distances of all particles to the center of the galaxy ([0, 0, 0])
def _distances(position, point):
    '''
    Evaluates distances of particles from a given point
    '''
    distances = np.sqrt(np.sum(np.square(position[i] - point[i]) for i in range(3)))
    return distances

dist = np.array([_distances(pos[i], np.array([0, 0, 0])) for i in range(len(pos))])


#We need to loop over masses inside certain radii to make an array of cumulative masses
dr = 0.1
radii = np.arange(dr, 30+dr, dr)


#Getting the mass inside each radii using the analytical models of the galaxy
#Masses in 1e+10solar mass, lengths in Kpc
M_halo = 55.0
a_halo = 15.0

M_disk = 2.3
Rd_disk = 2.1
z0_disk = 0.42

M_bulge = 0.58
a_bulge = 0.38

M_gas = 0.16
Rd_gas = 2.1
z0_gas = 0.2 * z0_disk

mass_halo_r = cumulative_hernquist(radii, M_halo, a_halo)
mass_disk_r = cumulative_exponential(radii, M_disk, Rd_disk, z0_disk)
mass_bulge_r = cumulative_hernquist(radii, M_bulge, a_bulge)
mass_gas_r = cumulative_exponential(radii, M_gas, Rd_gas, z0_gas)

masses_from_models = mass_halo_r + mass_disk_r + mass_bulge_r + mass_gas_r


#Loop to get the mass inside each radii using the initial snapshot of the galaxy
masses_from_snapshot = []
for r in radii:
    mask = dist < r
    mass_r = np.sum(mass[mask])
    masses_from_snapshot.append(mass_r)


#Plotting both analytical and 'numerical' cumulative masses for NGC2992
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(radii, masses_from_models, color='indigo', label='Cumulative mass from models')
ax.plot(radii, masses_from_snapshot, color='yellowgreen', label='Cumulative mass from snapshot')
plt.legend()
plt.show()
