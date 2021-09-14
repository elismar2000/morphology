import matplotlib.pyplot as plt
import numpy as np

#===============================
#Exponential
#===============================

def exp(a, I0, h):
    return I0 * np.exp(-a/h)

I0 = 3.64052
h = 13.5896

#================================
#Sérsic
#================================

def sersic(a, n, Ie, re):
    bn = 2*n - (1/3)
    return Ie * np.exp(-bn * ((a/re)**(1/n) - 1) )

n = 5.08936
Ie = 5.5241
re = 3.90792

#================================
#Plotting
#================================

a = np.arange(1, 20, 0.5)
plt.plot(a, exp(a, I0, h), '--', color='magenta', label='Exponential')
plt.plot(a, sersic(a, n, Ie, re), '--', color='orange', label='Sérsic')
plt.legend()
plt.xlabel('a')
plt.ylabel('brightness')
plt.show()
