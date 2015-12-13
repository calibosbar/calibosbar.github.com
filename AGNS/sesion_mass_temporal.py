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

plt.figure()
img_masked= np.ma.masked_less_equal(img, -5.0)
plt.imshow(img_masked, origin = 'lower')
img_masked.min()


plt.figure()
t34= img_masked[34,:]
plt.plot(t34)
t33= img_masked[33,:]
s33= t33+t34
plt.plot(s33)
t32= img_masked[32,:]
s32= t32+s33
plt.plot(s32)
t31= img_masked[31,:]
s31= t31+s32
plt.plot(s31)
t30= img_masked[30,:]
s30= t30+s31
plt.plot(s30)
t29= img_masked[29,:]
s29= t29+s30
plt.plot(s29)
t28= img_masked[28,:]
s28= t28+s29
plt.plot(s28)
t27= img_masked[27,:]
s27= t27+s28
plt.plot(s27)
t26= img_masked[26,:]
s26= t26+s27
plt.plot(s26)
t25= img_masked[25,:]
s25= t25+s26
plt.plot(s25)
t24= img_masked[24,:]
s24= t24+s25
plt.plot(s24)
t23= img_masked[23,:]
s23= t23+s24
plt.plot(s23)
t22= img_masked[22,:]
s22= t22+s23
plt.plot(s22)
t21= img_masked[21,:]
s21= t21+s22
plt.plot(s21)
t20= img_masked[20,:]
s20= t20+s21
plt.plot(s20)
t19= img_masked[19,:]
s19= t19+s20
plt.plot(s19)
t18= img_masked[18,:]
s18= t18+s19
plt.plot(s18)
t17= img_masked[17,:]
s17= t17+s18
plt.plot(s17)
t16= img_masked[16,:]
s16= t16+s17
plt.plot(s16)
t15= img_masked[15,:]
s15= t15+s16
plt.plot(s15)
t14= img_masked[14,:]
s14= t14+s15
plt.plot(s14)
t13= img_masked[13,:]
s13= t13+s14
plt.plot(s13)
t12= img_masked[12,:]
s12= t12+s13
plt.plot(s12)
t11= img_masked[11,:]
s11= t11+s12
plt.plot(s11)
t10= img_masked[10,:]
s10= t10+s11
plt.plot(s10)
t9= img_masked[9,:]
s9= t9+s10
plt.plot(s9)
t8= img_masked[8,:]
s8= t8+s9
plt.plot(s8)
t7= img_masked[7,:]
s7= t7+s8
plt.plot(s7)
t6= img_masked[6,:]
s6= t6+s7
plt.plot(s6)
t5= img_masked[5,:]
s5= t5+s6
plt.plot(s5)
t4= img_masked[4,:]
s4= t4+s5
plt.plot(s4)
t3= img_masked[3,:]
s3= t3+s4
plt.plot(s3)
t2= img_masked[2,:]
s2= t2+s3
plt.plot(s2)
t1= img_masked[1,:]
s1= t1+s2
plt.plot(s1)
t0= img_masked[0,:]
s0= t0+s1
plt.plot(s0)
