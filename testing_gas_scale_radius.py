import numpy as np
import matplotlib.pyplot as plt

from pygadgetreader import *


snapshot = '/home/elismar/Documentos/Fisica/IC/dice/ngc2992-93/ngc2993.g2'

pos = readsnap(snapshot, 'pos', 'gas')

mass = readsnap(snapshot, 'mass', 'gas')


M = np.sum(mass)
x = pos[:, 0]
y = pos[:, 1]
z = pos[:, 2]

x_cm = np.sum(x*mass)/M
y_cm = np.sum(y*mass)/M
z_cm = np.sum(z*mass)/M

cm = np.array([x_cm, y_cm, z_cm])

plt.plot(x, y, ',')
plt.plot(x_cm, y_cm, 'rX')
plt.show()


def _distances(position, point):
    '''
    Evaluates distances of particles from a given point
    '''
    distances = np.sqrt(np.sum(np.square(position[i] - point[i]) for i in range(3)))
    return distances


dists = np.array([_distances(pos[i], cm) for i in range(len(pos))])

h_disk = 0.8

mask_in = dists < h_disk
mask_out = dists > h_disk

mass_in = np.sum(mass[mask_in])
mass_out = np.sum(mass[mask_out])

print('Massa interna a h_disk = ', mass_in)
print('Massa externa a h_disk = ', mass_out)
