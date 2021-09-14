from lmfit import Model, Parameters
import numpy as np
import matplotlib.pyplot as plt

from counts2mass import abs_mag2mass
from px2kpc import px2kpc

#========================================

I0 = 64.75
r_out = px2kpc(8.0, ps=1.0, distance=38)
alpha = 1.0
beta = 1.47

def ferrer(r, I0, r_out, alpha, beta):
    I = I0 * (1 - (r/r_out)**(2 - beta))**alpha
    mass = abs_mag2mass(I)
    return mass

def hernquist(r, M, a):
    return (M / (2 * np.pi)) * (a/r) * (1/(r + a)**3)

#========================================

dr = r_out/1000
r = np.arange(dr, r_out, dr)

mass = ferrer(r, I0, r_out, alpha, beta)

model = Model(hernquist)
params = Parameters()
params.add('M', value=1e+10)
params.add('a', value=0.3)
result = model.fit(mass, params, r=r)

print(result.fit_report())

plt.plot(r, ferrer(r, I0, r_out, alpha, beta), 'bx', label='Ferrer')
plt.plot(r, hernquist(r, result.params['M'].value, result.params['a'].value), 'k-', label='Hernquist')
plt.legend()
plt.show()
