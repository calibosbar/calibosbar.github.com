import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.axes as ax
import seaborn as sns
from astropy.io import fits
import pylab

#plt.xkcd()
hdus = fits.open('UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
fig=plt.figure()
ax1 = fig.add_subplot(111)
img_mask= np.power(10, img)
#diff= np.diff(img_mask, axis=0)
#diff1= diff/0.37837838
img_masked= diff1
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0
ticksy=np.linspace(0,2,9)
ticksx=np.linspace(0,13,9)
ax1.set_xticklabels(ticksx)
ax1.set_yticklabels(ticksy)
plt.imshow(img_masked, origin = 'lower',aspect='auto') 
plt.pcolor(img_masked, norm=LogNorm())
plt.set_cmap('seismic')
cbar=plt.colorbar() #ticks=[0,1,2,3,4,5]
#cbar.ax.set_yticklabels(np.array([0.01,0.1,1,10,100]))
cbar.ax.set_ylabel('Log $\Sigma_{*}$ $[M_{sun}/pc^{2}]$', rotation=270, fontsize=10, verticalalignment='top')
plt.ylabel('$R/R_e$')
plt.xlabel('Look Back in time (Gigayears)')
plt.title('UGC11680NED01') 

outfile = 'UGC11680NED01_SFH_SFR_.fits.gz'

hdu = fits.PrimaryHDU(img_masked)
hdu.writeto(outfile, clobber=True)


