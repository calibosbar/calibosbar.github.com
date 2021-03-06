

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.axes as ax
import seaborn as sns
from astropy.io import fits
import pylab

#plt.xkcd()
hdus = fits.open('UGC11680NED01.p_e.rad_SFH_lum.fits.gz')
img = hdus[0].data

fig=plt.figure()
ax1 = fig.add_subplot(111)
#img_mask= np.power(10, img)
img_masked= img
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0
ticksy=np.linspace(0,2,9)
ticksx=np.linspace(0,14,9)
ax1.set_xticklabels(ticksx)
ax1.set_yticklabels(ticksy)
plt.imshow(img_masked, origin = 'lower',aspect='auto') 
plt.pcolor(img_masked,norm=LogNorm()) 
plt.set_cmap('bwr')
cbar=plt.colorbar() #ticks=[0,1,2,3,4,5]
#cbar.ax.set_yticklabels(np.array([0.01,0.1,1,10,100]))
cbar.ax.set_ylabel('Log $L_*$ $(L_{sun}kpc^{2})$', rotation=270, fontsize=10, verticalalignment='top')
plt.ylabel('$R/R_e$')
plt.xlabel('Look Back in time (Gigayears)')
plt.title('UGC11680NED01') 

fig=plt.figure()
x =np.linspace(0, 10, 39)
ax1= fig.add_subplot(111)
ax2= ax1.twiny()
ax1.set_ylabel('Light Fraction ')
ax1.set_xlabel('Log(Age/year)')
ax1.set_xticks(np.array([0,2,4,6,8,10]))
ax1.set_xticklabels(np.array([0, 2 , 4, 6, 8, 10]))
ax2.set_xticks(np.array([0,2,4,6,8,10]))
ax2.set_xticklabels(np.array([0.01, 0.25 , 0.5, 1.0, 2.0, 3.5,4.0]))
[l.set_rotation(45) for l in ax2.get_xticklabels()]
ax2.tick_params(axis='x', labelbottom='off')
ax2.set_xlabel('UGC11680NED01 (log) Redshift')
#plt.title('UGC11680NED01 MGH Temporal ')


t0= img_masked[0,:]
t1= img_masked[1,:]
t2= img_masked[2,:]
t3= img_masked[3,:]
t4= img_masked[4,:]
t5= img_masked[5,:]
t6= img_masked[6,:]
t7= img_masked[7,:]

ts1= (t0+t1+t2+t3+t4+t5+t6+t7)/8

plt.plot(x,ts1, label= '0.0 < R/Re < 0.5 ')
plt.legend(loc=1)

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

ts2= (t8+t9+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20+t21+t22+t23+t24+t25)/17
plt.plot(x,ts2, label= '0.5 < R/Re < 1 ')
plt.legend(loc=1)

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

ts3= (t26+t27+t28+t29+t30+t31+t32+t33+t34+t35)/10
plt.plot(x,ts3, label= '1 < R/Re < 2 ')
plt.legend(loc=1)





