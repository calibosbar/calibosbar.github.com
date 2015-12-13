

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
img_masked= img_mask
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
cbar.ax.set_ylabel('Log $\Sigma_{*}$ $[M_{\odot} kpc^{-2}]$', fontsize=10)
plt.ylabel('$R/R_e$')
plt.xlabel('Look Back in time (Gigayears)')
plt.title('UGC11680NED01') 


hdus2 = fits.open('all_agns_SFH_Mass.fits.gz')
img2 = hdus2[0].data

fig=plt.figure()
ax1 = fig.add_subplot(111)
img_masked2= img2
where_are_NaNs = isnan(img_masked2)
img_masked2[where_are_NaNs] = 0
ticksy=np.linspace(0,2,9)
ticksx=np.linspace(0,13,9)
ax1.set_xticklabels(ticksx)
ax1.set_yticklabels(ticksy)
plt.imshow(img_masked2, origin = 'lower',aspect='auto') 
plt.pcolor(img_masked2, norm=LogNorm())
plt.set_cmap('seismic')
cbar=plt.colorbar() #ticks=[0,1,2,3,4,5]
#cbar.ax.set_yticklabels(np.array([0.01,0.1,1,10,100]))
cbar.ax.set_ylabel('Log $\Sigma_{*}$ $[M_{\odot} kpc^{-2}]$', fontsize=10)
plt.ylabel('$R/R_e$')
plt.xlabel('Look Back in time (Gigayears)')
plt.title('All AGNs') 

sum_chi= ((img_masked-img_masked2)**2)

norm=0.5*(img_masked+img_masked2)

sum_tot= sum_chi/norm


chi_j= np.sum(sum_tot, axis=1)

chi_ugc11680= np.sum(chi_j)

img_masked.shape

norm_df= 36*38

chi_red_ugc11680= chi_ugc11680/norm_df





chiugc= img_masked[:,range(32,39)]
#np.savetxt("chiugc.csv", chiugc, delimiter=",")

chiagn= img_masked2[:,range(32,39)]
#np.savetxt("chiagn.csv", chiagn, delimiter=",")

chiugc2= img_masked[:,range(23,31)]
#np.savetxt("chiugc2.csv", chiugc2, delimiter=",")

chiagn2= img_masked2[:,range(23,31)]
#np.savetxt("chiagn2.csv", chiagn2, delimiter=",")

chiugc3= img_masked[:,range(15,22)]
#np.savetxt("chiugc3.csv", chiugc3, delimiter=",")

chiagn3= img_masked2[:,range(15,22)]
#np.savetxt("chiagn3.csv", chiagn3, delimiter=",")


chiugc4= img_masked[:,range(7,14)]
#np.savetxt("chiugc4.csv", chiugc4, delimiter=",")

chiagn4= img_masked2[:,range(7,14)]
#np.savetxt("chiagn4.csv", chiagn4, delimiter=",")

chiugc5= img_masked[:,range(0,6)]
#np.savetxt("chiugc5.csv", chiugc5, delimiter=",")

chiagn5= img_masked2[:,range(0,6)]
#np.savetxt("chiagn5.csv", chiagn5, delimiter=",")
