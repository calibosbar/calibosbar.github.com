
source actiate py34
source deactiate

import numpy as np
from astropy.io import fits
import pylab

image_list = ['/home/jeffrey/tesis/mass_11_2/ARP220.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/IC0307.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/IC0485.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/IC0674.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/IC4566.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC0036.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC0169.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC0180.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC0192.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC0426.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC0477.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC0508.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC0787.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC1070.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC1324.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC2410.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC2449.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC2565.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC2572.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC2639.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC2916.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC4003.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC5406.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC5533.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC5616.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC5635.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC5675.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC5720.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC5772.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC5888.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC6060.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC6154.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC6394.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC6478.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC6497.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC6941.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC7364.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC7466.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/NGC7591.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC00036.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC01368.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC02099.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC02367.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC03038.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC03973.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC03995.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC04132.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC04262.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC04455.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC05108.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC05111.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC08107.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC08234.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC08781.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC09401.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC09537.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC10205.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC10337.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC10710.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC10811.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC11717.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC12185.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC12250.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/UGC12348.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/mass_11_2_chi.p_e.rad_SFH_lum_Mass.fits.gz',
'/home/jeffrey/tesis/mass_11_2/MCG-02-03-015.p_e.rad_SFH_lum_Mass.fits.gz']


image_concat = []

for image in image_list:
    image_concat.append(fits.getdata(image))

img_concat_lin=(np.power(10, image_concat))/66
where_are_NaNs = np.isnan(img_concat_lin)
img_concat_lin[where_are_NaNs] = 0


final_image = np.zeros(shape=img_concat_lin[0].shape)


for image in img_concat_lin:
    final_image += image


outfile = 'all_mass_11_2_chi.p_e.rad_SFH_lum_Mass.fits.gz'

hdu = fits.PrimaryHDU(final_image)
hdu.writeto(outfile, clobber=True)
