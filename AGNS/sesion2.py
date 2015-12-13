# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits

image_list = [ 'AGN_1.fits.gz','AGN_2.fits.gz','AGN_3.fits.gz','AGN_4.fits.gz','AGN_5.fits.gz','AGN_6.fits.gz','AGN_7.fits.gz','AGN_8.fits.gz','AGN_9.fits.gz','AGN_10.fits.gz','AGN_11.fits.gz','AGN_12.fits.gz','AGN_13.fits.gz','AGN_14.fits.gz','AGN_15.fits.gz','AGN_16.fits.gz','AGN_17.fits.gz','AGN_18.fits.gz','AGN_19.fits.gz','AGN_20.fits.gz','AGN_21.fits.gz','AGN_22.fits.gz','AGN_23.fits.gz','AGN_24.fits.gz']
image_concat =[]
for image in image_list:
    image_concat.append(fits.getdata(image))
    
final_image = np.zeros(shape=image_concat[0].shape)
for image in image_concat:
    final_image += image
    
image_hist = plt.hist(final_image.flat, 1000)
image_hist = plt.hist(final_image.flat, 30)
image_hist = plt.hist(final_image.flat, 1)
image_hist = plt.hist(final_image.flat, 100000)
image_hist = plt.hist(final_image.flat)
plt.imshow()
plt.imshow(final_image)
plt.imshow(final_image)
outfile = 'all_AGN_V.fits.gz'
hdu = fits.PrimaryHDU(SFH_radial_V)
hdu = fits.PrimaryHDU(final_image)
hdu.writeto(outfile,clobber = True)
get_ipython().magic(u'save sesion2 1-32')
