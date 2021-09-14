import matplotlib.pyplot as plt
import numpy as np

from lmfit import Model, Parameters

#================================

data = np.genfromtxt('disk_profile.dat')

data[:, 0] = abs(data[:, 0] - data[np.argmax(data[:, 1]), 0])

lower = 6
upper = 20
mask = (data[:, 0] >= lower) & (data[:, 0] <= upper)
folded = data[mask]

i = (upper - lower) + 1
mean = (folded[0:i, 1][::-1] + folded[i:, 1]) / 2
a = folded[i:, 0]

fig, axs = plt.subplots(2, 1, figsize=(40, 40))

axs[0].plot(data[:, 0], data[:, 1], '--o', color='magenta')
axs[0].plot(folded[:, 0], folded[:, 1], 'o', color='blue')
axs[0].plot(a, mean, 'X', color='yellow')

def exp(a, I0, h):
    return I0 * np.exp(-a/h)

model = Model(exp)
params = Parameters()
params.add('I0', value=11.0)
params.add('h', value=7.0)
result = model.fit(mean, params, a=a)

print(result.fit_report())

x = np.arange(0, 21, 0.5)
axs[1].plot(a, mean, 'X', color='yellow')
axs[1].plot(x, exp(x, result.params['I0'].value, result.params['h'].value), '--', color='blue')

plt.show()
