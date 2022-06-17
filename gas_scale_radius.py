import matplotlib.pyplot as plt
import numpy as np

h_disk = 0.9

def _f1(h_gas):
    return -h_disk/h_gas

def _f2(h_gas):
    return np.log(h_gas / (2*(h_gas + h_disk)))

r = np.linspace(0, 5*h_disk, 1000)

f1 = _f1(r)
f2 = _f2(r)

plt.plot(r, f1, '--', label='Function 1')
plt.plot(r, f2, '--', label='Function 2')
plt.axvline(h_disk, color='red')
plt.legend()
plt.show()
