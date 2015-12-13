



import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as ax
from astropy.io import fits
from matplotlib.colors import LogNorm
import pylab

plt.xkcd()
hdus = fits.open('all_fo_mass_8_SFH_Mass.fits.gz')
img = hdus[0].data
#img_mask= np.power(10, img)
img_masked= img
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0


fig=plt.figure()
x =np.logspace(0, 14, 39)
ax1= fig.add_subplot(111)
#ax2= ax1.twiny()
#z= 0.026
ax1.set_ylabel('$M(t)/M_0$ ')
ax1.set_xlabel('Look back in time (Gigayears)')
ax1.set_xticks(np.array([0,2,4,6,8,10,12,14]))
ax1.set_xticklabels(np.array([0, 2 , 4, 6, 8, 10,12,14]))
#ax2.set_xticks(np.array([0,2,4,6,8,10,12,14]))
#ax2.set_xticklabels(np.array([0.1, 0.25 , 0.5, 1.0, 2.0, 3.5,4.0,4.2]))
[l.set_rotation(45) for l in ax2.get_xticklabels()]
#ax2.tick_params(axis='x', labelbottom='off')
#ax1.axhline(y=0.8,xmin=0,xmax=3,c="black",linewidth=1,zorder=0)
#ax2.set_xlabel('Redshift')



#plt.ylimt(0,4)
#plt.yscale('log')
#plt.xscale('log')
plt.title(' MGH $0<R/R_e<0.5$  ')


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
plt.plot(x,s1, label= 'Face On Mass 8 ')
plt.legend(loc=1)

hdus = fits.open('all_fo_mass_9_SFH_Mass.fits.gz')
img = hdus[0].data
#img_mask= np.power(10, img)
img_masked= img
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0

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
plt.plot(x,s1, label= 'Face On Mass 9')
plt.legend(loc=1)

hdus = fits.open('all_fo_mass_10_SFH_Mass.fits.gz')
img = hdus[0].data
#img_mask= np.power(10, img)
img_masked= img
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0

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
plt.plot(x,s1, label= 'Face On Mass 10')
plt.legend(loc=1)

hdus = fits.open('all_fo_mass_11_SFH_Mass.fits.gz')
img = hdus[0].data
#img_mask= np.power(10, img)
img_masked= img
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0

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
plt.plot(x,s1, label= 'Face On Mass 11')
plt.legend(loc=1)



hdus = fits.open('UGC11680NED01_Mass_lin.fits.gz')
img = hdus[0].data
#img_mask= np.power(10, img)
img_masked= img
where_are_NaNs = isnan(img_masked)
img_masked[where_are_NaNs] = 0

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
plt.plot(x,s1, label= 'UGC11680NED01 ')
plt.legend(loc=1)




