
source actiate py34
source deactiate

import numpy as np
from astropy.io import fits
import pylab

image_list = ['/home/jeffrey/Documents/theo/MCG-02-08-014.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/theo/UGC03973.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/theo/UGC03995.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/theo/UGC09711.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/Documents/theo/UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz',]
image_concat = []

for image in image_list:
    image_concat.append(fits.getdata(image))

img_concat_lin= (np.power(10, image_concat))/5
where_are_NaNs = np.isnan(img_concat_lin)
img_concat_lin[where_are_NaNs] = 0

final_image = np.zeros(shape=img_concat_lin[0].shape)


for image in img_concat_lin:
    final_image += image

outfile = 'theos1.p_e.rad_SFH_lim_Mass.fits.gz'

hdu = fits.PrimaryHDU(final_image)
hdu.writeto(outfile, clobber=True)
