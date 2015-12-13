

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits



plt.xkcd()
hdus = fits.open('IC0540.p_e.rad_SFH_lum_Mass.fits.gz')
img = hdus[0].data
plt.imshow(img)
img = hdus[0].data
plt.imshow(img)
plt.clf()
plt.imshow(img, origin = 'lower')
plt.colorbar()
img.shape
img.min()

plt.figure()
img_mask= np.ma.masked_less_equal(img, -5.0)
plt.imshow(img_mask, origin = 'lower')
plt.colorbar()
plt.title('IC0540 Lum Mass Map') 
img_mask.min()
img_masked= img_mask*0.282051282


plt.figure()
plt.ylabel('Log(M/Mo)')
plt.xlabel('Log(Age/year)')
plt.yscale('log')
plt.title('IC0540 Mass Grow History Temporal')
x = linspace (0, 11, 39)


t0= img_masked[0,:]
plt.plot(x,t0)
t1= img_masked[1,:]
s1= t1+ t0
plt.plot(x,s1)
t2= img_masked[2,:]
s2= t2+ s1
plt.plot(x,s2)
t3= img_masked[3,:]
s3= t3+ s2
plt.plot(x,s3)
t4= img_masked[4,:]
s4= t4+ s3
plt.plot(x,s4)
t5= img_masked[5,:]
s5= t5+ s4
plt.plot(x,s5)
t6= img_masked[6,:]
s6= t6+ s5
plt.plot(x,s6)
t7= img_masked[7,:]
s7= t7+ s6
plt.plot(x,s7)
t8= img_masked[8,:]
s8= t8+ s7
plt.plot(x,s8)
t9= img_masked[9,:]
s9= t9+ s8
plt.plot(x,s9)
t10= img_masked[10,:]
s10= t10+ s9
plt.plot(x,s10)


