# coding: utf-8
plt.imshow(img)
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
t0= img[:,0] 
plt.plot(t0)
t1= img[:,1]
s1= t1 + t0
plt.plot(s1)
plt.plot(t1)
t2= img[:,2]
plt.plot(t2)
plt.plot(t3)
t3= img[:,3]
plt.plot(t3)
t4= img[:,4]
plt.plot(t4)
t5= img[:,5]
plt.plot(t5)
t6= img[:,6]
plt.plot(t6)
t7= img[:,7]
plt.plot(t7)
t8= img[:,8]
plt.plot(t8)
t9= img[:,9]
plt.plot(t9)
t10= img[:,10]
plt.plot(t10)
t11= img[:,11]
plt.plot(t11)
t12= img[:,12]
plt.plot(t12)
s12 = t12+ t0
plt.plot(s12)
plt.plot(t0)
t13= img[:,13]
s13 = t13+ s12
plt.plot(s13)
t14= img[:,14]
s14= t14+s13
plt.plot(s14)
t15= img[:,15]
s15= t15+s14
plt.plot(s15)
t16= img[:,16]
s16= t16+s15
plt.plot(s16)
get_ipython().magic(u'save sesion 8 1-57')
