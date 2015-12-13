

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
plt.title(' MGH $0.5<R/R_e<1$  ')

t8= img_masked[8,:]
t9= img_masked[9,:]
t10= img_masked[10,:]
t11= img_masked[11,:]
t12= img_masked[12,:]
t13= img_masked[13,:]
t14= img_masked[14,:]
t15= img_masked[15,:]
t16= img_masked[16,:]
t17= img_masked[17,:]
t18= img_masked[18,:]
t19= img_masked[19,:]
t20= img_masked[20,:]
t21= img_masked[21,:]
t22= img_masked[22,:]
t23= img_masked[23,:]
t24= img_masked[24,:]
t25= img_masked[25,:]


ts2= t8+t9+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20+t21+t22+t23+t24+t25
ts2[:] = ts2[::-1].cumsum()
ts2[:] = ts2[::-1]
s2= ts2/np.max(ts2)
plt.plot(x,s2, label= 'NGC6394 ')
plt.legend(loc=1)


hdus = fits.open('UGC11680NED02.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
img_mask= np.power(10, img)
img_masked= img_mask
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0


t8= img_masked[8,:]
t9= img_masked[9,:]
t10= img_masked[10,:]
t11= img_masked[11,:]
t12= img_masked[12,:]
t13= img_masked[13,:]
t14= img_masked[14,:]
t15= img_masked[15,:]
t16= img_masked[16,:]
t17= img_masked[17,:]
t18= img_masked[18,:]
t19= img_masked[19,:]
t20= img_masked[20,:]
t21= img_masked[21,:]
t22= img_masked[22,:]
t23= img_masked[23,:]
t24= img_masked[24,:]
t25= img_masked[25,:]


ts2= t8+t9+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20+t21+t22+t23+t24+t25
ts2[:] = ts2[::-1].cumsum()
ts2[:] = ts2[::-1]
s2= ts2/np.max(ts2)
plt.plot(x,s2, label= 'UGC11680NED02')
plt.legend(loc=1)

hdus = fits.open('NGC3160.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
img_mask= np.power(10, img)
img_masked= img_mask
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0


t8= img_masked[8,:]
t9= img_masked[9,:]
t10= img_masked[10,:]
t11= img_masked[11,:]
t12= img_masked[12,:]
t13= img_masked[13,:]
t14= img_masked[14,:]
t15= img_masked[15,:]
t16= img_masked[16,:]
t17= img_masked[17,:]
t18= img_masked[18,:]
t19= img_masked[19,:]
t20= img_masked[20,:]
t21= img_masked[21,:]
t22= img_masked[22,:]
t23= img_masked[23,:]
t24= img_masked[24,:]
t25= img_masked[25,:]


ts2= t8+t9+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20+t21+t22+t23+t24+t25
ts2[:] = ts2[::-1].cumsum()
ts2[:] = ts2[::-1]
s2= ts2/np.max(ts2)
plt.plot(x,s2, label= 'NGC3160 ')
plt.legend(loc=1)

hdus = fits.open('NGC5635.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
img_mask= np.power(10, img)
img_masked= img_mask
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0


t8= img_masked[8,:]
t9= img_masked[9,:]
t10= img_masked[10,:]
t11= img_masked[11,:]
t12= img_masked[12,:]
t13= img_masked[13,:]
t14= img_masked[14,:]
t15= img_masked[15,:]
t16= img_masked[16,:]
t17= img_masked[17,:]
t18= img_masked[18,:]
t19= img_masked[19,:]
t20= img_masked[20,:]
t21= img_masked[21,:]
t22= img_masked[22,:]
t23= img_masked[23,:]
t24= img_masked[24,:]
t25= img_masked[25,:]


ts2= t8+t9+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20+t21+t22+t23+t24+t25
ts2[:] = ts2[::-1].cumsum()
ts2[:] = ts2[::-1]
s2= ts2/np.max(ts2)
plt.plot(x,s2, label= 'NGC5635 ')
plt.legend(loc=1)

hdus = fits.open('UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
img_mask= np.power(10, img)
img_masked= img_mask
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0


t8= img_masked[8,:]
t9= img_masked[9,:]
t10= img_masked[10,:]
t11= img_masked[11,:]
t12= img_masked[12,:]
t13= img_masked[13,:]
t14= img_masked[14,:]
t15= img_masked[15,:]
t16= img_masked[16,:]
t17= img_masked[17,:]
t18= img_masked[18,:]
t19= img_masked[19,:]
t20= img_masked[20,:]
t21= img_masked[21,:]
t22= img_masked[22,:]
t23= img_masked[23,:]
t24= img_masked[24,:]
t25= img_masked[25,:]


ts2= t8+t9+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20+t21+t22+t23+t24+t25
ts2[:] = ts2[::-1].cumsum()
ts2[:] = ts2[::-1]
s2= ts2/np.max(ts2)
plt.plot(x,s2, label= 'UGC11680NED01')
plt.legend(loc=1)







