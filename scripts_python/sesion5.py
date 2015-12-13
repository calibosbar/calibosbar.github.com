# coding: utf-8
get_ipython().magic(u'ls ')
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
hdus = fits.open('UGC11680NED01.p_e.rad_SFH_lum.fits.gz')
img = hdus[0].data
plt.imshow(img)
img = hdus[0].data
plt.imshow(img)
plt.clf()
plt.imshow(img, origin = 'lower')
img.shape
img.min()
lum= img[:,0].cumsum()
plt.figure()
plt.plot(lum)
lum= img[:,1].cumsum()
plt.plot(lum)
lum= img[:,2].cumsum()
plt.plot(lum)
lum= img[:,3].cumsum()
plt.plot(lum)
lum= img[:,4].cumsum()
plt.plot(lum)
lum= img[:,5].cumsum()
plt.plot(lum)
lum= img[:,6].cumsum()
plt.plot(lum)
lum= img[:,7].cumsum()
plt.plot(lum)
lum= img[0,:].cumsum()
plt.plot(lum)
lum= img[0,:].cumsum()
plt.plot(lum)
lum= img[1,:].cumsum()
plt.plot(lum)
lum= img[2,:].cumsum()
plt.plot(lum)
lum= img[3,:].cumsum()
plt.plot(lum)
lum= img[4,:].cumsum()
plt.plot(lum)
lum= img[0,:]
plt.plot(lum)
lum= img[:,0]
plt.plot(lum)
lum= img[:,1]
plt.plot(lum)
plt.plot(img[:, 4:10])
plt.plot(img[:, 20:30])
lum= img[:,39]
lum= img[:,38]
plt.plot(lum)
lum= img[:,38].cumsum()
plt.plot(lum)
lum1= img.sum(axis=1)
plt.plot(lum1)
plt.plot(lum1)
lum2= img.sum(axis=0)
plt.plot(lum2)
plt.plot(lum2)
lum= img[:,38]
plt.plot(lum)
plt.imshow(img)
plt.imshow(img, origin = 'lower')
lum= img[36,:]
lum= img[35,:]
plt.plot(lum)
plt.figure()
lum= img[:,38]
plt.plot(lum)
plt.figure()
lum= img[:,38].cumsum()
plt.plot(lum)
lum= img[:,38].sum()
plt.plot(lum)
plt.figure()
lum= img[:,38].sum()
plt.plot(lum)
lum= img[:,38]
plt.plot(lum)
plt.figure()
lum= img[35,:]
plt.plot(lum)
lum= img[35,:].cumsum()
plt.plot(lum)
lum= img[:,38].cumsum()
plt.figure()
plt.plot(lum)
lum= img[34,:].cumsum()
plt.plot(lum)
lum= img[33,:].cumsum()
lum= img[32,:].cumsum()
plt.plot(lum)
lum= img[31,:].cumsum()
plt.plot(lum)
lum= img[30,:].cumsum()
plt.plot(lum)
lum= img[29,:].cumsum()
plt.plot(lum)
lum= img[28,:].cumsum()
plt.plot(lum)
lum= img[27,:].cumsum()
plt.plot(lum)
lum= img[26,:].cumsum()
plt.plot(lum)
lum= img[25,:].cumsum()
plt.plot(lum)
plt.plot(lum)
lum= img[24,:].cumsum()
plt.plot(lum)
lum= img[23,:].cumsum()
plt.plot(lum)
lum= img[22,:].cumsum()
plt.plot(lum)
lum= img[21,:].cumsum()
plt.plot(lum)
lum= img[20,:].cumsum()
plt.plot(lum)
lum= img[20,:].cumsum()
lum
len(lum)
lum.shape
lum1= img[20,:]
lum1
lum1.shape
lum= img[20,:].cumsum()
lum
plt.plot(lum)
hdus = fits.open('UGC11680NED01.p_e.rad_SFH_lum_V.fits.gz')
img = hdus[0].data
plt.imshow(img)
img = hdus[0].data
plt.imshow(img)
plt.clf()
plt.imshow(img, origin = 'lower')
img.shape
img.min()
lum= img[20,:].cumsum()
plt.plot(lum)
plt.plot(lum)
lum= img[20,:]
plt.plot(lum)
lum1= img[:,20]
plt.figure()
plt.plot(lum1)
lum1= img[:,20].cumsum()
plt.plot(lum1)
lum1= img[:,39]
lum1= img[:,38]
plt.plot(lum1)
lum1= img[:,37]
plt.plot(lum1)
lum1= img[:,36]
plt.plot(lum1)
lum1= img[:,35]
plt.plot(lum1)
lum1= img[:,34]
plt.plot(lum1)
plt.plot(lum1)
lum1= img[:,33]
plt.plot(lum1)
lum1= img[:,32]
plt.plot(lum1)
lum1= img[:,31]
plt.plot(lum1)
lum1= img[:,30]
plt.plot(lum1)
lum1= img[:,0]
plt.plot(lum1)
lum1= img[:,0]
plot_lum= plt.plot(lum1)
plot_lum= plt.plot(lum1)
import matplotlibticker as ticker
import matplotlib.ticker as ticker
ticks= plot_lum.get_xticks()*0.05
scale = 0.05
ticks =ticker.FuncFormatter(lambda x, pos '{0:g}'.format(x*scale))
ticks =ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x*scale))
plot_lum.xaxis.set_major_formatter(ticks)
plot_lum.set_major_formatter(ticks)
lum1= img[:,1]
plt.plot(lum1)
lum1= img[:,2]
plt.plot(lum1)
lum1= img[:,3]
plt.plot(lum1)
lum1= img[:,4]
plt.plot(lum1)
lum1= img[:,5]
plt.plot(lum1)
plt.plot(lum1)
lum1= img[:,10]
plt.plot(lum1)
plt.plot(lum1 , extent=[6,10])
lum1= img[:,10]
plt.plot(lum1)
plt.plot.xkcd(lum1)
plt.xkcd(lum1)
plt.xkcd.plot(lum1)
plt.plot(lum1)
plt.plot(lum1)
lum1= img[:,10]
plt.plot(lum1)
lum1= img[:,10]
plt.plot(lum1)
lum1= img[:,10]
plt.plot(lum1)
lum1= img[:,10]
plt.plot(lum1)
lum1= img[:,10, 1000]
lum1= img[:,10]
plt.plot(lum1)
plt.plot(lum1)
lum1= img[:,10]
plt.plot(lum1)
lum1= img[:,11]
plt.plot(lum1)
lum1= img[:,39]
lum1= img[:,38]
plt.plot(lum1)
lum1= img[:,19]
plt.plot(lum1)
lum1= img[:,1]
plt.plot(lum1)
lum1= img[:,38]
plt.plot(lum1)
lum1= img[:,19]
plt.plot(lum1)
lum1= img[:,0]
plt.plot(lum1)
lum1= img[:,38]
plt.plot(lum1 label= $t= 10.1$)
plt.plot(lum1, label= $t= 10.1$)
plt.plot(lum1, label= 10.1 Gyrs)
plt.plot(lum1, label='t=10.1 Gyrs')
plt.plot(lum1, label='10.1 Gyrs')
lum1= img[:,19]
plt.plot(lum1, label='5.1 Gyrs')
lum1= img[:,19]
plt.plot(lum1, label='5.1 Gyrs')
plt.plot(lum1, label='5.1 Gyrs')
plt.legend()
lum38= img[:,38]
lum19= img[:,19]
lum19= img[:,25]
lum0= img[:,0]
lum25= img[:,25]
lum19= img[:,19]
plt.plot(lum38, label='10.1 Gyrs')
plt.plot(lum25, label='9 Gyrs')
plt.plot(lum19, label='8 Gyrs')
plt.plot(lum0, label='6.5 Gyrs')
plt.legend()
plt.plot(lum0, label='6.5 Gyrs')
plt.plot(lum0, label='6.5 Gyrs')
plt.plot(lum19, label='8 Gyrs')
plt.plot(lum25, label='9 Gyrs')
plt.plot(lum38, label='10.1 Gyrs')
plt.plot(lum38, label='10.1 Log(age/yr)')
plt.plot(lum25, label='9 Log(age/yr)')
plt.plot(lum19, label='8 Log(age/yr)')
plt.plot(lum0, label='Log(age/yr)')
plt.legend()
get_ipython().magic(u'save sesion5 1-265')
