
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as ax
from astropy.io import fits
from matplotlib.colors import LogNorm
import pylab

#plt.xkcd()
hdus = fits.open('NGC6394.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
img_mask= np.power(10, img)
img_masked= img_mask
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0


fig=plt.figure()
x =np.logspace(0, 14, 39)
ax1= fig.add_subplot(111)
ax2= ax1.twiny()
#z= 0.026
ax1.set_ylabel('$M(t)/M_0$ ')
ax1.set_xlabel('Look back in time (Gigayears)')
ax1.set_xticks(np.array([0,2,4,6,8,10,12,14]))
ax1.set_xticklabels(np.array([0, 2 , 4, 6, 8, 10,12,14]))
ax2.set_xticks(np.array([0,2,4,6,8,10,12,14]))
#ax2.set_xticklabels(np.array([0.1, 0.25 , 0.5, 1.0, 2.0, 3.5,4.0,4.2]))
[l.set_rotation(45) for l in ax2.get_xticklabels()]
ax2.tick_params(axis='x', labelbottom='off')
#ax1.axhline(y=0.8,xmin=0,xmax=3,c="black",linewidth=1,zorder=0)
#ax2.set_xlabel('Redshift')



#plt.ylimt(0,4)
#plt.yscale('log')
#plt.xscale('log')
plt.title(' MGH $1<R/R_e<2$  ')

t26= img_masked[26,:]
t27= img_masked[27,:]
t28= img_masked[28,:]
t29= img_masked[29,:]
t30= img_masked[30,:]
t31= img_masked[31,:]
t32= img_masked[32,:]
t33= img_masked[33,:]
t34= img_masked[34,:]
t35= img_masked[35,:]


ts3= t26+t27+t28+t29+t30+t31+t32+t33+t34+t35
ts3[:] = ts3[::-1].cumsum()
ts3[:] = ts3[::-1]
s3= ts3/np.max(ts3)
plt.plot(x,s3, label= 'NGC6394 ')
plt.legend(loc=1)

hdus = fits.open('UGC11680NED02.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
img_mask= np.power(10, img)
img_masked= img_mask
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0


t26= img_masked[26,:]
t27= img_masked[27,:]
t28= img_masked[28,:]
t29= img_masked[29,:]
t30= img_masked[30,:]
t31= img_masked[31,:]
t32= img_masked[32,:]
t33= img_masked[33,:]
t34= img_masked[34,:]
t35= img_masked[35,:]


ts3= t26+t27+t28+t29+t30+t31+t32+t33+t34+t35
ts3[:] = ts3[::-1].cumsum()
ts3[:] = ts3[::-1]
s3= ts3/np.max(ts3)
plt.plot(x,s3, label= 'UGC11680NED02 ')
plt.legend(loc=1)


hdus = fits.open('NGC3160.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
img_mask= np.power(10, img)
img_masked= img_mask
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0


t26= img_masked[26,:]
t27= img_masked[27,:]
t28= img_masked[28,:]
t29= img_masked[29,:]
t30= img_masked[30,:]
t31= img_masked[31,:]
t32= img_masked[32,:]
t33= img_masked[33,:]
t34= img_masked[34,:]
t35= img_masked[35,:]


ts3= t26+t27+t28+t29+t30+t31+t32+t33+t34+t35
ts3[:] = ts3[::-1].cumsum()
ts3[:] = ts3[::-1]
s3= ts3/np.max(ts3)
plt.plot(x,s3, label= 'NGC3160 ')
plt.legend(loc=1)

hdus = fits.open('NGC5635.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
img_mask= np.power(10, img)
img_masked= img_mask
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0


t26= img_masked[26,:]
t27= img_masked[27,:]
t28= img_masked[28,:]
t29= img_masked[29,:]
t30= img_masked[30,:]
t31= img_masked[31,:]
t32= img_masked[32,:]
t33= img_masked[33,:]
t34= img_masked[34,:]
t35= img_masked[35,:]


ts3= t26+t27+t28+t29+t30+t31+t32+t33+t34+t35
ts3[:] = ts3[::-1].cumsum()
ts3[:] = ts3[::-1]
s3= ts3/np.max(ts3)
plt.plot(x,s3, label= 'NGC5635 ')
plt.legend(loc=1)

hdus = fits.open('UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
img_mask= np.power(10, img)
img_masked= img_mask
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0


t26= img_masked[26,:]
t27= img_masked[27,:]
t28= img_masked[28,:]
t29= img_masked[29,:]
t30= img_masked[30,:]
t31= img_masked[31,:]
t32= img_masked[32,:]
t33= img_masked[33,:]
t34= img_masked[34,:]
t35= img_masked[35,:]


ts3= t26+t27+t28+t29+t30+t31+t32+t33+t34+t35
ts3[:] = ts3[::-1].cumsum()
ts3[:] = ts3[::-1]
s3= ts3/np.max(ts3)
plt.plot(x,s3, label= 'UGC11680NED01 ')
plt.legend(loc=1)




