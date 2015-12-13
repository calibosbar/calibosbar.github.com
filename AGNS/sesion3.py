# coding: utf-8
get_ipython().magic(u'cd Desktop')
get_ipython().magic(u'cd AGNS')
hdus = fits.open('all_AGN_V.fits.gz')
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
hdus = fits.open('all_AGN_V.fits.gz')
img = hdus[0].data
plt.imshow(img, origin = 'lower')
plt.figure()
profile = img.sum(axis=1)
plt.plot(profile)
get_ipython().magic(u'save sesion3 1-13')
