

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as ax
from astropy.io import fits
from matplotlib.colors import LogNorm
import pylab
import seaborn as sn

#plt.xkcd()
hdus = fits.open('UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
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
ax1.set_ylabel('$M(t)/M_{T}$ ')
ax1.set_xlabel('Look back in time (Gyrs)')
ax1.set_xticks(np.array([0,2,4,6,8,10,12,14]))
ax1.set_xticklabels(np.array([0, 2 , 4, 6, 8, 10,12,14]))
ax2.set_xticks(np.array([0,2,4,6,8,10,12,14]))
ax2.set_xticklabels(np.linspace(0,6,14))
[l.set_rotation(45) for l in ax2.get_xticklabels()]
ax2.tick_params(axis='x', labelbottom='off')
#ax1.axhline(y=0.8,xmin=0,xmax=3,c="black",linewidth=1,zorder=0)
#ax2.set_xlabel('Redshift')



#plt.ylimt(0,4)
#plt.yscale('log')
#plt.xscale('log')
plt.title('UGC11680NED01 ')


t0= img_masked[0,:]
t1= img_masked[1,:]
t2= img_masked[2,:]
t3= img_masked[3,:]
t4= img_masked[4,:]
t5= img_masked[5,:]
t6= img_masked[6,:]
t7= img_masked[7,:]



ts1= t0+t1+t2+t3+t4+t5+t6+t7
ts1[:] = ts1[::-1].cumsum()
ts1[:] = ts1[::-1]
s1= ts1/np.max(ts1)
plt.plot(x,s1, label= '0.0 < R/Re < 0.5 ')
plt.legend(loc=3)

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
plt.plot(x,s2, label= '0.5 < R/Re < 1 ')
plt.legend(loc=3)



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
plt.plot(x,s3, label= '1 < R/Re < 2 ')
plt.legend(loc=3)












