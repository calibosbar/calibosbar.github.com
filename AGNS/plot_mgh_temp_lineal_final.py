


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as ax
from astropy.io import fits
import pylab
import seaborn as sn

#plt.xkcd()
hdus = fits.open('/home/jeffrey/Documents/dataproducts/NGC5394.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
img_mask= np.power(10, img)
img_masked= img_mask
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0

fig=plt.figure()
x =np.logspace(0, 14, 39)
ax1 = fig.add_subplot(111)
#x_lin= np.power(10, x)
ax1.set_ylabel('$M(t)/M_{0}$')
ax1.set_xlabel('Look back in time (Gyrs)')
ticksx=np.array([0,3,7,10,12,14])
ax1.set_xticklabels(ticksx)
#plt.ylimt(0,4)
#plt.yscale('log')
#plt.xscale('log')
ax1.set_title('NGC5394')


#t0= img_masked[0,:]
#t0[:] = t0[::-1].cumsum()
#t0[:] = t0[::-1]
#s0= t0/np.max(t0)
#plt.plot(x,s0, label= '0.0 R/Re' )
#plt.legend(loc=3)

t1= img_masked[3,:]
t1[:] = t1[::-1].cumsum()
t1[:] = t1[::-1]
#s1=t1
s1= t1/np.max(t1)
plt.plot(x ,s1, label= '0.0 < R/Re < 0.5' )
plt.legend(loc=3)

t20= img_masked[20,:]
t20[:] = t20[::-1].cumsum()
t20[:] = t20[::-1]
s20= t20/np.max(t20)
plt.plot(x,s20, label= '0.5 < R/Re < 1.0' )
plt.legend(loc=3)


t31= img_masked[31,:]
t31[:] = t31[::-1].cumsum()
t31[:] = t31[::-1]
s31= t31/np.max(t31)
plt.plot(x,s31, label= '1 < R/Re < 2' )
plt.legend(loc=3)






