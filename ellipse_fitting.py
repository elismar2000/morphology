from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

from photutils.datasets import make_noise_image
from photutils.isophote import EllipseGeometry, Ellipse
from photutils import EllipticalAperture

#================================

fits_file = '/home/elismar/Documentos/Fisica/IC/imfit-1.7.1/ngc2992-93/images/GMOS/ngc2992/ngc2992_model3.fits'

data = fits.getdata(fits_file)

#================================

geometry = EllipseGeometry(x0=237, y0=168, sma=87, eps=0.3, pa=-58.*np.pi/180)

aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,
                           geometry.sma*(1 - geometry.eps), geometry.pa)

norm = plt.Normalize(0, 5000)
plt.imshow(data, origin='lower', norm=norm)
aper.plot(color='white')
plt.show()

ellipse = Ellipse(data, geometry)
isolist = ellipse.fit_image()
