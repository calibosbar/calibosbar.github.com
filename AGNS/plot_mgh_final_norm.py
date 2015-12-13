

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pylab



#plt.xkcd()
hdus = fits.open('UGC11680NED01.p_e.rad_SFH_lum_Mass.fits.gz')
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

img_mask= np.power(10, img)
#img_mask= np.ma.masked_less_equal(img, -5.0)
img_masked= img_mask #*0.055555556
plt.imshow(img_masked, origin = 'lower')
plt.colorbar()
plt.ylabel('r/re (pixels)')
plt.xlabel('Age/yrs (pixels) ')
plt.title('UGC11680NED01 Lum Mass Map') 
img_masked.min()




plt.figure()
x =np.linspace(0, 2, 36)
plt.ylabel('Log Density Stellar Mass')
plt.xlabel('r/re')
#plt.yscale('log')
#plt.xscale('log')
plt.title('UGC11680NED01 Mass Grow History Radial')

t38= img_masked[:,38]
s38 = t38/np.linalg.norm(t38)
plt.plot(x,s38)


t37= img_masked[:,37]
t37n= t37+s38
s37= t37n/np.linalg.norm(t37n)
plt.plot(x,s37)

t36= img_masked[:,36]
t36n= t36+s37
s36= t36n/np.linalg.norm(t36n)
plt.plot(x,s36)

t35= img_masked[:,35]
t35n= t35+s36
s35= t35n/np.linalg.norm(t35n)
plt.plot(x,s35)

t34= img_masked[:,34]
t34n= t34+s35
s34= t34n/np.linalg.norm(t34n)
plt.plot(x,s34)


t33= img_masked[:,33]
t33n= t33+s33
s33= t33n/np.linalg.norm(t33n)
plt.plot(x,s33)


t32= img_masked[:,32]
s32= t32+s
plt.plot(x,s32)
t31= img_masked[:,31]
s31= t31+s32
plt.plot(x,s31)
t30= img_masked[:,30]
s30= t30+s31
plt.plot(x,s30)
t29= img_masked[:,29]
s29= t29+s30
plt.plot(x,s29)
t28= img_masked[:,28]
s28= t28+s29
plt.plot(x,s28)
t27= img_masked[:,27]
s27= t27+s28
plt.plot(x,s27)
t26= img_masked[:,26]
s26= t26+s27
plt.plot(x,s26)
t25= img_masked[:,25]
s25= t25+s26
plt.plot(x,s25)
t24= img_masked[:,24]
s24= t24+s25
plt.plot(x,s24)
t23= img_masked[:,23]
s23= t23+s24
plt.plot(x,s23)
t22= img_masked[:,22]
s22= t22+s23
plt.plot(x,s22)
t21= img_masked[:,21]
s21= t21+s22
plt.plot(x,s21)
t20= img_masked[:,20]
s20= t20+s21
plt.plot(x,s20)
t19= img_masked[:,19]
s19= t19+s20
plt.plot(x,s19)
t18= img_masked[:,18]
s18= t18+s19
plt.plot(x,s18)
t17= img_masked[:,17]
s17= t17+s18
plt.plot(x,s17)
t16= img_masked[:,16]
s16= t16+s17
plt.plot(x,s16)
t15= img_masked[:,15]
s15= t15+s16
plt.plot(x,s15)
t14= img_masked[:,14]
s14= t14+s15
plt.plot(x,s14)
t13= img_masked[:,13]
s13= t13+s14
plt.plot(x,s13)
t12= img_masked[:,12]
s12= t12+s13
plt.plot(x,s12)
t11= img_masked[:,11]
s11= t11+s12
plt.plot(x,s11)
t10= img_masked[:,10]
s10= t10+s11
plt.plot(x,s10)
t9= img_masked[:,9]
s9= t9+s10
plt.plot(x,s9)
t8= img_masked[:,8]
s8= t8+s9
plt.plot(x,s8)
t7= img_masked[:,7]
s7= t7+s8
plt.plot(x,s7)
t6= img_masked[:,6]
s6= t6+s7
plt.plot(x,s6)
t5= img_masked[:,5]
s5= t5+s6
plt.plot(x,s5)
t4= img_masked[:,4]
s4= t4+s5
plt.plot(x,s4)
t3= img_masked[:,3]
s3= t3+s4
plt.plot(x,s3)
t2= img_masked[:,2]
s2= t2+s3
plt.plot(x,s2)
t1= img_masked[:,1]
s1= t1+s2
plt.plot(x,s1)
t0= img_masked[:,0]
s0= t0+s1
plt.plot(x,s0)




