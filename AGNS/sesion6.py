# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
hdus = fits.open('UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
plt.imshow(img)
img = hdus[0].data
plt.imshow(img)
plt.clf()
plt.imshow(img, origin = 'lower')
img.shape
img.min()
t38= img[:,38]
plt.figure()
plt.plot(t38)

mass= img[:,33].cumsum()
plt.plot(mass)
mass= img[:,32].cumsum()
plt.plot(mass)
mass= img[:,31].cumsum()
plt.plot(mass)
mass= img[:,30].cumsum()
plt.plot(mass)

mass= img[:,:]
plt.plot(mass)
mass_sum= img[11,:].cumsum()
mass_sum_norm= img /mass_sum[:, np.newaxis]
mass_sum= img.sum(axis=0)
mass_sum_norm = img/mass_sum[:, np.newaxis]
mass_sum= img.sum(axis=1)
mass_sum_norm = img/mass_sum[:, np.newaxis]
plt.plot(mass_sum_norm)
norm= mass_sum_norm[11,:].cumsum()
plt.plot(norm)
norm= mass_sum_norm[:,0].cumsum()
plt.plot(norm)
norm= mass_sum_norm[:,1].cumsum()
plt.plot(norm)
norm= mass_sum_norm[:,2].cumsum()
plt.plot(norm)
norm= mass_sum_norm[:,3].cumsum()
plt.plot(norm)
norm= mass_sum_norm[:,3].cumsum()
mass_sum= img.sum(axis=0)
mass_sum_norm = img/mass_sum[np.newaxis,:]
plt.plot(mass_sum_norm)
norm= mass_sum_norm[:,0].cumsum()
plt.plot(norm)
norm= mass_sum_norm[:,1].cumsum()
plt.plot(norm)
norm= mass_sum_norm[:,2].cumsum()
plt.plot(norm)

norm= mass_sum_norm[:,38].cumsum()
norm1= mass_sum_norm[:,38].cumsum()
plt.plot(norm1 label= '10.1 log (Age/yr)')
plt.plot(norm1, label= '10.1 log (Age/yr)')
norm2= mass_sum_norm[:,35].cumsum()
plt.plot(norm2, label= '9.1 log (Age/yr)')
norm3= mass_sum_norm[:,32].cumsum()
plt.plot(norm3, label= '8.8 log (Age/yr)')
norm4= mass_sum_norm[:,31].cumsum()
plt.plot(norm3, label= '8.3 log (Age/yr)')
norm4= mass_sum_norm[:,29].cumsum()
plt.plot(norm4, label= '7.8 log (Age/yr)')
plt.legend()
norm= mass_sum_norm[:,:].cumsum()
plt.plot(norm)
norm= mass_sum_norm[1,:].cumsum()
plt.plot(norm)

norm= mass_sum_norm[0,:].cumsum()
plt.plot(norm)
sum_time = img.sum(axis=1)
sum_time_norm = img / sum_time[:, np.newaxis]
plt.plot(sum_time_norm)
time_norm=sum_time_norm[0,:].cumsum()
plt.plot(time_norm)
time_norm=sum_time_norm[:,0].cumsum()
plt.plot(time_norm)
time_norm=sum_time_norm[:,35].cumsum()
plt.plot(time_norm)

time_norm1=sum_time_norm[8,:].cumsum()
plt.plot(time_norm1,label='0.6 R/Re' )
time_norm2=sum_time_norm[9,:].cumsum()
plt.plot(time_norm2,label='1.1 R/Re' )
time_norm3=sum_time_norm[10,:].cumsum()
plt.plot(time_norm3,label='2.1 R/Re' )
plt.legend()
time=img[0,:].cumsum
plt.plot(time)

time1=img[8,:].cumsum()
plt.plot(time1, label='0.6 R/Re')
time2=img[9,:].cumsum()
plt.plot(time2, label='1.1 R/Re')
time2=img[10,:].cumsum()
plt.plot(time3, label='1.6 R/Re')
time3=img[10,:].cumsum()
plt.plot(time3, label='1.6 R/Re')
plot.legend()
plt.legend()
get_ipython().magic(u'save sesion6 1-405')
