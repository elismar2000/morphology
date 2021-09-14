from astropy.modeling.models import Gaussian2D
import matplotlib.pyplot as plt
import numpy as np

from photutils.datasets import make_noise_image
from photutils.isophote import EllipseGeometry, Ellipse
from photutils import EllipticalAperture

g = Gaussian2D(100., 75, 75, 20, 12, theta=40.*np.pi/180.)
ny = nx = 150
y, x = np.mgrid[0:ny, 0:nx]

noise = make_noise_image((ny, nx), distribution='gaussian', mean=0.,
                         stddev=2., random_state=12345)
data = g(x, y) + noise

#==========================================

geometry = EllipseGeometry(x0=75, y0=75, sma=20, eps=0.5,
                           pa=20.*np.pi/180)

aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma,
                           geometry.sma*(1 - geometry.eps), geometry.pa)

plt.imshow(data, origin='lower')
aper.plot(color='white')
plt.show()

ellipse = Ellipse(data, geometry)
isolist = ellipse.fit_image()
