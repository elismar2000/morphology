from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

from photutils.datasets import make_noise_image
from photutils.isophote import EllipseGeometry, Ellipse
from photutils import EllipticalAperture

#===================================
#Reading fits file containing data to be fitted
#===================================

bulge = fits.getdata('/home/elismar/Documentos/Fisica/IC/imfit-1.7.1/ngc2992-93/images/NICMOS/n4sb08040/bulge.fits')

#===================================
#Fitting elliptical isophotes to bulge
#===================================

geometry = EllipseGeometry(x0=50, y0=65, sma=3, eps=0.1, pa=0)

ellipse = Ellipse(bulge, geometry)
isolist = ellipse.fit_image()

#===================================
#Plotting results
#===================================

fig = plt.figure()
axs = fig.add_subplot(111)
norm = plt.Normalize(-2, 10)
axs.imshow(bulge, origin='lower', norm=norm)
axs.set_title('Data')

smas = np.linspace(1, 10, 10)
for sma in smas:
    iso = isolist.get_closest(sma)
    x, y, = iso.sampled_coordinates()
    axs.plot(x, y, color='white')

plt.show()
