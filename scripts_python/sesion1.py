# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
hdus = fits.open('UGC11680NED01.p_e.rad_SFH_lum.fits.gz')
get_ipython().magic(u'pinfo hdus')
img = hdus[0].data
plt.imshow(img)
plt.clf()
plt.imshow(img, origin = 'lower')
img.shape
img.min
img.min()
plt.figure()
plt.plot(img[:,25])
plt.figure()
plt.plot(img[:,30])
plt.figure()
plt.plot(img[:,35])
plt.plot(img[:,35])
plt.figure()
plt.plot(img[:,40])
plt.plot(img[:,37])
plt.plot(img[:,36])
plt.plot(img[:,35])
plt.figure()
plt.plot(img[:,35])
plt.figure()
profile = img.sum(axis=1)
plt.figure()
plt.plot(profile)
spectrum =img.sum(axis=0)
plt.figure()
plt.plot(spectrum)
time =img.sum(axis=2)
time =img.sum(axis=2)
spectrum =img.sum(axis=0)
plt.plot(spectrum)
plt.plot(spectrum)
plt.plot(profile)
plt.plot(profile)
get_ipython().magic(u'save sesion1 1-42')