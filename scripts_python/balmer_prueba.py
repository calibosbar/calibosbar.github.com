

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.axes as ax
import seaborn as sns
from astropy.io import fits
import pylab

#plt.xkcd()
hdus = fits.open('map_Ha.UGC11680.fits')
img = hdus[0].data

fig=plt.figure()
ax1 = fig.add_subplot(111)
img_masked= img
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0

plt.imshow(img_masked, origin = 'lower',aspect='auto') 
plt.pcolor(img_masked,norm=LogNorm()) 
plt.set_cmap('seismic')
cbar=plt.colorbar() 

hdus1 = fits.open('map_Hb.UGC11680.fits')
img1 = hdus1[0].data

fig=plt.figure()
ax1 = fig.add_subplot(111)
img_masked1= img1
where_are_NaNs = isnan(img_masked1)
img_masked1[where_are_NaNs] = 0

plt.imshow(img_masked1, origin = 'lower',aspect='auto') 
plt.pcolor(img_masked1,norm=LogNorm()) 
plt.set_cmap('seismic')
cbar=plt.colorbar() 


with np.errstate(divide='ignore', invalid='ignore'):
    ratio = np.true_divide(img_masked,img_masked1)
    ratio[ratio == np.inf] = 0
    ratio= np.nan_to_num(ratio)

 ratio_norm=ratio/2.76
  balmer = 2.5*(ma.log10(ratio_norm))

plt.imshow(balmer, origin = 'lower',aspect='auto') 
plt.pcolor(balmer,norm=LogNorm()) 
plt.set_cmap('jet')
cbar=plt.colorbar() 


ratio_norm=ratio/2.76
balmer = 2.5*(ma.log(ratio_norm))
print balmer.filled(0)
outfile = 'balmer.fits.gz'

hdu = fits.PrimaryHDU(balmer)
hdu.writeto(outfile, clobber=True)

hdus2 = fits.open('ratio.fits')
img2 = hdus2[0].data

fig=plt.figure()
ax1 = fig.add_subplot(111)
img_masked2= img2
where_are_NaNs = isnan(img_masked2)
img_masked2[where_are_NaNs] = 0

plt.imshow(img_masked2, origin = 'lower',aspect='auto') 
plt.pcolor(img_masked2,norm=LogNorm()) 
plt.set_cmap('seismic')
cbar=plt.colorbar() 









